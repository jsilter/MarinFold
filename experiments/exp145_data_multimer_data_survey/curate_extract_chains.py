"""Split downloaded dimers into one structure file per distinct chain.

PINDER clusters at the level of the individual chain, not the complex: it runs
Foldseek chain-against-chain, builds a graph over chains, finds communities, and
only then labels a complex by the pair of communities its two chains fall into.
Our fetch produced whole two-chain mmCIFs, so this is the step in between.

The unit of identity here is the UniProt accession. A homodimer's two chains are
the same accession and contribute one entry; a heterodimer contributes two. Over
the 881,010 val-eligible complexes that is 858,360 distinct chains.

Each accession is emitted once, from the lexicographically smallest complex that
contains it, so a rerun with the same inputs produces byte-identical output and a
partial run can be resumed by simply skipping accessions already on disk.
"""

import argparse
import csv
import gzip
import io
import json
import os
import tarfile
import time

import pyarrow.fs
from pathlib import Path

# Zero-based positions in an _atom_site row. The AFDB writes these columns in a
# fixed order (verified against the release), so a positional read is safe and is
# roughly twenty times faster than a real mmCIF parser at this volume. Field 6 is
# label_asym_id; do not use field 7, which is label_entity_id and is legitimately
# "1" for both chains of a homodimer.
ASYM_FIELD = 6

ASSIGNMENT = "dimer_split_assignment_id30_cov50.csv.gz"


def choose_sources(path: Path, split: str) -> dict[str, tuple[str, str]]:
    """Map each accession to the (complex, chain) it will be extracted from.

    Returns ``{accession: (model_entity_id, "A" | "B")}``. Chain A carries
    ``uniprot_1`` and chain B carries ``uniprot_2``; for a homodimer the two are
    equal and only chain A is ever requested.
    """
    opener = gzip.open if path.suffix == ".gz" else open
    sources: dict[str, tuple[str, str]] = {}
    with opener(path, "rt") as handle:
        for row in csv.DictReader(handle):
            if row["split"] != split:
                continue
            model = row["model_entity_id"]
            for accession, chain in ((row["uniprot_1"], "A"), (row["uniprot_2"], "B")):
                prior = sources.get(accession)
                if prior is None or model < prior[0]:
                    sources[accession] = (model, chain)
    return sources


def split_chains(text: str, wanted: set[str]) -> dict[str, str]:
    """Pull the requested chains out of one dimer mmCIF.

    Emits a minimal file per chain: the ``data_`` header, the ``_atom_site`` loop
    header, and that chain's ATOM records. Foldseek reads only the coordinates and
    the residue names, so carrying the entity and reference blocks across would
    cost bytes without changing any alignment.
    """
    lines = text.splitlines()
    header = next(line for line in lines if line.startswith("data_"))
    loop = [line for line in lines if line.startswith("_atom_site.")]
    if not loop:
        raise ValueError(f"{header}: no _atom_site loop")

    rows: dict[str, list[str]] = {chain: [] for chain in wanted}
    for line in lines:
        if not line.startswith("ATOM"):
            continue
        chain = line.split()[ASYM_FIELD]
        if chain in rows:
            rows[chain].append(line)

    out: dict[str, str] = {}
    for chain, atoms in rows.items():
        if not atoms:
            raise ValueError(f"{header}: chain {chain} has no ATOM records")
        out[chain] = "\n".join([header, "loop_", *loop, *atoms, "#", ""])
    return out


def iter_tar(stream, name: str):
    """Yield (model_entity_id, mmCIF text) from one open tar stream.

    Mode "r|" is the streaming reader: it never seeks, so the same code path works
    on a local file and on an S3 object being read once from the front.
    """
    with tarfile.open(fileobj=stream, mode="r|") as tar:
        for member in tar:
            if not member.name.endswith(".cif.gz"):
                continue
            handle = tar.extractfile(member)
            if handle is None:
                raise ValueError(f"{name}: {member.name} is not a regular file")
            yield member.name.removesuffix(".cif.gz"), gzip.decompress(handle.read()).decode()


def iter_members(source: str, shard_stride: int = 1):
    """Yield (model_entity_id, mmCIF text) from a local path or an s3:// prefix.

    ``shard_stride`` takes every Nth tar rather than the first N. Shards are
    contiguous bands of a sorted ID list, and the survey found AFDB ID windows are
    organism-blocked, so consecutive shards would sample one slice of taxonomy.
    Striding spreads the sample across the whole range.
    """
    if source.startswith("s3://"):
        # pyarrow, never boto3: boto3 wedges ephemeral EC2 workers about a minute
        # in. resolve_s3_region avoids a cross-region read of every shard.
        path = source[len("s3://"):].rstrip("/")
        bucket = path.split("/")[0]
        fs = pyarrow.fs.S3FileSystem(region=pyarrow.fs.resolve_s3_region(bucket))
        listing = fs.get_file_info(pyarrow.fs.FileSelector(path, recursive=False))
        tars = sorted(f.path for f in listing if f.path.endswith(".tar"))[::shard_stride]
        print(f"{len(tars)} tar shards from {source} (stride {shard_stride})")
        started = time.time()
        for index, name in enumerate(tars):
            # A run this long is unobservable without a progress line: the first
            # attempt spent 8 hours in here with nothing in the log but heartbeats.
            print(f"[{time.time() - started:8.0f}s] shard {index + 1}/{len(tars)} "
                  f"{name.rsplit('/', 1)[-1]}", flush=True)
            # Read each tar as one object rather than streaming it. tarfile pulls
            # 10240-byte records, and open_input_stream turns every one of those
            # into its own round trip: measured 0.08 MB/s against 12 MB/s for
            # whole-object reads of the same tar, a 47x difference that showed up
            # in production as an extraction pass on track to take 74 hours. The
            # shards are ~250 MB each, so holding one in memory costs nothing.
            with fs.open_input_file(name) as handle:
                blob = handle.read()
            yield from iter_tar(io.BytesIO(blob), name)
        return

    root = Path(source)
    tars = sorted(root.glob("*.tar"))[::shard_stride]
    if not tars:
        for path in sorted(root.glob("*.cif.gz")):
            yield path.name.removesuffix(".cif.gz"), gzip.decompress(path.read_bytes()).decode()
        return
    for tar_path in tars:
        with tar_path.open("rb") as stream:
            yield from iter_tar(stream, tar_path.name)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True,
                        help="directory of .cif.gz files, a directory of .tar shards, "
                             "or an s3:// prefix holding the shards")
    parser.add_argument("--shard-stride", type=int, default=1,
                        help="take every Nth shard; shards are contiguous ID bands and "
                             "AFDB ID windows are organism-blocked, so a subsample must "
                             "stride rather than truncate")
    parser.add_argument("--assignment", type=Path,
                        default=Path("data/curation") / ASSIGNMENT)
    parser.add_argument("--split", default="val")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--limit", type=int, default=0,
                        help="stop after this many chains (0 = no limit)")
    args = parser.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)
    print(f"reading {args.assignment}")
    sources = choose_sources(args.assignment, args.split)
    print(f"{len(sources):,} distinct chains in split {args.split!r}")

    # Invert: complex -> the chains this complex is the chosen source for. Most
    # complexes source one or two chains; many source none, because an earlier
    # complex already covered both accessions.
    by_model: dict[str, dict[str, str]] = {}
    for accession, (model, chain) in sources.items():
        by_model.setdefault(model, {})[chain] = accession

    written = skipped = 0
    started = time.time()
    for model, chains in iter_members(args.source, args.shard_stride):
        wanted = by_model.get(model)
        if not wanted:
            continue
        need = {chain: acc for chain, acc in wanted.items()
                if not (args.out / f"{acc}.cif.gz").exists()}
        skipped += len(wanted) - len(need)
        if not need:
            continue
        for chain, text in split_chains(chains, set(need)).items():
            # Level 6, not the default 9: measured on this data, level 9 costs
            # 25.8 ms per structure against 6.1 and produces a marginally larger
            # file. Over 880,248 structures that is 5 hours against 1.5. mtime=0
            # keeps the output byte-identical across reruns.
            (args.out / f"{need[chain]}.cif.gz").write_bytes(
                gzip.compress(text.encode(), compresslevel=6, mtime=0))
            written += 1
            if written % 50_000 == 0:
                print(f"  {written:,} chains written, {skipped:,} already on disk", flush=True)
        if args.limit and written + skipped >= args.limit:
            break

    elapsed = time.time() - started
    total = sum(p.stat().st_size for p in args.out.glob("*.cif.gz"))
    summary = {
        "split": args.split,
        "source": args.source,
        "shard_stride": args.shard_stride,
        "n_chains_in_split": len(sources),
        "n_written": written,
        "n_already_present": skipped,
        "bytes_on_disk": total,
        "mean_kb": round(total / max(written + skipped, 1) / 1024, 1),
        "elapsed_seconds": round(elapsed, 1),
    }
    print(json.dumps(summary, indent=2))
    (args.out.parent / f"{args.out.name}_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
