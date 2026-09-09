"""Measure how Foldseek chain-against-chain search scales toward PINDER's 50k.

Phase 3 of CURATION_PLAN.md. PINDER partitions chains into sub-databases of
50,000 and runs every ordered pair of sub-databases as an independent job, so the
cost of the whole clustering is (cost of one pair) x (number of pairs). The only
Foldseek numbers this experiment had came from 45 structures, and nothing about
the 50k scale follows from them.

This runs an all-against-all on the same chain set at several sizes and reports
wall time and peak RSS for each, which is what turns "one pair costs X" into a
fan-out width. Search cost is expected to grow super-linearly (each added chain
is both an extra query and an extra target), so the deliverable is the measured
exponent, not an assumption about it.

Parameters are PINDER's, from the table in CURATION_PLAN.md: -s 11.0, -e 0.05,
--max-seqs 1000, --alignment-type 2, lddt scoring.
"""

import argparse
import json
import os
import re
import resource
import shutil
import subprocess
import time
from pathlib import Path

import foldseek_env

# PINDER's FoldseekConfig. -e 0.05 is fifty times looser than Foldseek's own
# default and is deliberate: the graph wants distant hits, because a similarity
# too weak to cluster on is still strong enough to leak between splits.
SENSITIVITY = "11.0"
EVALUE = "0.05"
MAX_SEQS = "1000"
ALIGNMENT_TYPE = "2"

# lddt is the score the graph is built on. qcov/tcov feed the 0.5 coverage
# filter and alnlen feeds the 10-residue minimum, so all three are collected
# here rather than recomputed later.
FORMAT_OUTPUT = "query,target,fident,alnlen,qstart,qend,tstart,tend,evalue,bits,lddt,qcov,tcov"


def peak_rss_gb() -> float:
    """Peak RSS of all child processes so far, in GB.

    ru_maxrss is a high-water mark that never decreases, so this is the peak
    across every Foldseek invocation to date, not the peak of the last one. The
    caller takes differences only where a stage is known to dominate.
    """
    return resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / (1024 ** 2)


def assert_all_chains_ingested(db: Path, stdout: str, expected: int) -> None:
    """Fail if Foldseek's database holds fewer chains than were handed to it.

    Foldseek reports what it dropped and then exits 0 whether that number is zero
    or not, so a truncated database silently produces a clustering of whatever it
    managed to read. Two independent checks: the "Ignore N out of M" line, and the
    line count of the .lookup file, which is one row per stored entry and is the
    database's own account of itself rather than a message about it.
    """
    match = re.search(r"Ignore (\d+) out of (\d+)\.", stdout)
    if match is None:
        raise RuntimeError("createdb printed no 'Ignore N out of M' line, cannot verify")
    ignored, seen = int(match.group(1)), int(match.group(2))
    lookup = db.with_suffix(".lookup")
    stored = sum(1 for _ in lookup.open()) if lookup.exists() else -1
    if seen != expected or ignored != 0 or stored != expected:
        raise RuntimeError(
            f"createdb was given {expected} chain files but saw {seen}, ignored "
            f"{ignored} (too short / not protein), and stored {stored} entries in "
            f"{lookup.name}. A truncated database still exits 0."
        )


def run(binary: str, args: list[str], log: Path) -> tuple[float, str]:
    """Run one Foldseek subcommand, returning (seconds, stdout).

    Foldseek writes its progress and its ingestion counts to stdout, not stderr.
    """
    started = time.time()
    proc = subprocess.run([binary, *args], capture_output=True, text=True)
    elapsed = time.time() - started
    log.write_text(proc.stdout + "\n----- stderr -----\n" + proc.stderr)
    if proc.returncode != 0:
        raise RuntimeError(f"foldseek {args[0]} exited {proc.returncode}, see {log}")
    return elapsed, proc.stdout


def measure(binary: str, chains: list[Path], n: int, work: Path, threads: int) -> dict:
    """Build a database of the first n chains and search it against itself."""
    stage = work / f"n{n}_t{threads}"
    if stage.exists():
        shutil.rmtree(stage)
    inputs = stage / "input"
    inputs.mkdir(parents=True)
    for path in chains[:n]:
        # Hard links, so a 50k-chain stage costs no extra disk and no copy time.
        (inputs / path.name).hardlink_to(path)

    db = stage / "db"
    before = peak_rss_gb()
    db_seconds, stdout = run(binary, ["createdb", str(inputs), str(db)], stage / "createdb.log")
    assert_all_chains_ingested(db, stdout, n)
    db_rss = max(peak_rss_gb() - before, 0.0)

    hits = stage / "hits.tsv"
    search_seconds, _ = run(binary, [
        "easy-search", str(inputs), str(inputs), str(hits), str(stage / "tmp"),
        "-s", SENSITIVITY, "-e", EVALUE, "--max-seqs", MAX_SEQS,
        "--alignment-type", ALIGNMENT_TYPE, "--format-output", FORMAT_OUTPUT,
        "--threads", str(threads),
    ], stage / "search.log")

    n_hits = sum(1 for _ in hits.open())
    db_bytes = sum(p.stat().st_size for p in stage.glob("db*") if p.is_file())
    row = {
        "n_chains": n,
        "threads": threads,
        "createdb_seconds": round(db_seconds, 1),
        "search_seconds": round(search_seconds, 1),
        "total_seconds": round(db_seconds + search_seconds, 1),
        "peak_rss_gb": round(peak_rss_gb(), 2),
        "createdb_rss_gb": round(db_rss, 2),
        "db_bytes": db_bytes,
        "n_hits": n_hits,
        "hits_per_chain": round(n_hits / n, 1),
        "comparisons": n * n,
        "microseconds_per_comparison": round(search_seconds * 1e6 / (n * n), 3),
    }
    print(json.dumps(row))
    # The hit table for the largest size feeds phase 4; the tmp trees do not.
    shutil.rmtree(stage / "tmp", ignore_errors=True)
    return row


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chains", type=Path, required=True,
                        help="directory of per-chain .cif.gz written by curate_extract_chains.py")
    parser.add_argument("--sizes", type=int, nargs="+", default=[500, 1000, 2000, 4000])
    parser.add_argument("--threads", type=int, nargs="+", default=[8],
                        help="thread counts to measure; more than one sweeps parallel "
                             "efficiency, which is what lets a measurement on 8 cores "
                             "predict a 32-vCPU instance")
    parser.add_argument("--tag", default="",
                        help="suffix for the output file, to keep sweeps side by side")
    parser.add_argument("--work", type=Path, default=Path("data/curation/foldseek_scaling"))
    parser.add_argument("--out", type=Path, default=Path("data/curation"))
    args = parser.parse_args()

    binary = foldseek_env.ensure_foldseek()
    version = foldseek_env.foldseek_version()
    chains = sorted(args.chains.glob("*.cif.gz"))
    print(f"foldseek {version}, {len(chains):,} chains available, "
          f"threads {args.threads}, sizes {args.sizes}")

    sizes = [n for n in args.sizes if n <= len(chains)]
    if not sizes:
        raise SystemExit(f"no requested size fits {len(chains)} chains")
    args.work.mkdir(parents=True, exist_ok=True)

    # Descending thread order so the slowest configuration runs last and partial
    # output is still the useful half of the curve.
    rows = [measure(binary, chains, n, args.work, threads)
            for n in sizes for threads in sorted(args.threads, reverse=True)]

    summary = {"foldseek_version": version, "cpu_count": os.cpu_count(),
               "n_chains_available": len(chains), "rows": rows}
    name = f"foldseek_scaling{args.tag}.json"
    (args.out / name).write_text(json.dumps(summary, indent=2) + "\n")
    print(f"\nwrote {args.out / name}")


if __name__ == "__main__":
    main()
