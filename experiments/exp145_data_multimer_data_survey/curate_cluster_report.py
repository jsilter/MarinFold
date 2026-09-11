"""Summarize a Foldseek chain clustering into small committable artifacts.

Input is the two-column `easy-cluster` TSV (representative, member), gzipped,
as published by the phase 3 run to
s3://marinfold-exp91-usw2/MarinFold/exp145-afdb-dimers/cluster/chain_cluster.tsv.gz.

The TSV itself is 857,599 rows and does not belong in git. This writes three
things that do: the full size distribution (one row per distinct cluster size),
the largest clusters by name, and a summary JSON holding the numbers phase 4 and
phase 5 actually plan against.
"""

import argparse
import csv
import gzip
import json
from collections import Counter
from pathlib import Path


def read_clusters(path: Path) -> dict[str, int]:
    """Map each representative to its member count, from the easy-cluster TSV."""
    sizes: Counter[str] = Counter()
    with gzip.open(path, "rt") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if not row:
                continue
            sizes[row[0]] += 1
    return dict(sizes)


def coverage_points(sizes: list[int], fractions: tuple[float, ...]) -> dict[str, int]:
    """How many of the largest clusters it takes to cover each fraction of chains.

    Phase 5 draws one complex per cluster pair, so the shape of the head decides
    how much of the corpus a small draw can represent.
    """
    total = sum(sizes)
    points: dict[str, int] = {}
    running = 0
    ordered = sorted(sizes, reverse=True)
    targets = sorted(fractions)
    index = 0
    for rank, size in enumerate(ordered, start=1):
        running += size
        while index < len(targets) and running >= targets[index] * total:
            points[f"clusters_to_cover_{targets[index]:.2f}"] = rank
            index += 1
    return points


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tsv", type=Path, required=True,
                        help="gzipped easy-cluster TSV (representative, member)")
    parser.add_argument("--out", type=Path, default=Path("data/curation"))
    parser.add_argument("--top", type=int, default=25,
                        help="how many of the largest clusters to name")
    args = parser.parse_args()

    sizes = read_clusters(args.tsv)
    counts = sorted(sizes.values(), reverse=True)
    n_chains = sum(counts)
    histogram = Counter(counts)

    dist_path = args.out / "chain_cluster_sizes.csv"
    with open(dist_path, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["cluster_size", "n_clusters", "n_chains", "cumulative_chains"])
        running = 0
        for size in sorted(histogram, reverse=True):
            chains = size * histogram[size]
            running += chains
            writer.writerow([size, histogram[size], chains, running])

    top_path = args.out / "chain_cluster_largest.csv"
    ordered = sorted(sizes.items(), key=lambda kv: (-kv[1], kv[0]))
    with open(top_path, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["rank", "representative", "size", "share_of_chains"])
        for rank, (rep, size) in enumerate(ordered[:args.top], start=1):
            writer.writerow([rank, rep, size, round(size / n_chains, 6)])

    midpoint = sorted(counts)
    summary = {
        "n_chains": n_chains,
        "n_clusters": len(counts),
        "representative_fraction": round(len(counts) / n_chains, 6),
        "largest_cluster": counts[0],
        "median_cluster_size": midpoint[len(midpoint) // 2],
        "mean_cluster_size": round(n_chains / len(counts), 2),
        "n_singletons": histogram[1],
        "singleton_share_of_clusters": round(histogram[1] / len(counts), 4),
        "singleton_share_of_chains": round(histogram[1] / n_chains, 6),
        "n_clusters_ge_100": sum(1 for size in counts if size >= 100),
        "chains_in_clusters_ge_100": sum(size for size in counts if size >= 100),
        **coverage_points(counts, (0.50, 0.80, 0.90, 0.99)),
    }
    summary_path = args.out / "chain_cluster_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    print(f"wrote {dist_path}, {top_path}, {summary_path}")


if __name__ == "__main__":
    main()
