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


def iter_members(source: Path):
    """Yield (model_entity_id, mmCIF text) from a directory or a directory of tars."""
    if source.is_dir() and not any(source.glob("*.tar")):
        for path in sorted(source.glob("*.cif.gz")):
            yield path.name.removesuffix(".cif.gz"), gzip.decompress(path.read_bytes()).decode()
        return
    for tar_path in sorted(source.glob("*.tar")):
        with tarfile.open(tar_path) as tar:
            for member in tar:
                if not member.name.endswith(".cif.gz"):
                    continue
                handle = tar.extractfile(member)
                assert handle is not None
                payload = gzip.decompress(handle.read()).decode()
                yield member.name.removesuffix(".cif.gz"), payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True,
                        help="directory of .cif.gz files, or a directory of .tar shards")
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
    for model, chains in iter_members(args.source):
        wanted = by_model.get(model)
        if not wanted:
            continue
        need = {chain: acc for chain, acc in wanted.items()
                if not (args.out / f"{acc}.cif.gz").exists()}
        skipped += len(wanted) - len(need)
        if not need:
            continue
        for chain, text in split_chains(chains, set(need)).items():
            (args.out / f"{need[chain]}.cif.gz").write_bytes(
                gzip.compress(text.encode(), mtime=0))
            written += 1
        if args.limit and written + skipped >= args.limit:
            break

    elapsed = time.time() - started
    total = sum(p.stat().st_size for p in args.out.glob("*.cif.gz"))
    summary = {
        "split": args.split,
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
