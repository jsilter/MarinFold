# Copyright The MarinFold Authors
# SPDX-License-Identifier: Apache-2.0

"""B1 layer ablation: R-precision lost when one attention or MLP output is mean-ablated.

Reads ``b1_layer_ablation.py``'s ``scores/*.csv`` and ``timings/*.csv``
(``--input-dir``), concatenates them into ``data/b1_scores.csv`` and
``data/b1_timings.csv``, and writes ``data/b1_deltas.csv``.

Per condition, separation band (all, short, medium, long) and split (all
proteins, each eval set, each length tercile) the quantity is

    delta = R-precision(ablated) - R-precision(clean)

on the seed-averaged P(contact) matrix, the production readout, averaged
over proteins. Ablated and clean see the same proteins and document seeds,
so the comparison is paired.

Two CIs (``N_BOOT``, seeded):

- ``lo``/``hi``: bootstrap over proteins on the seed-averaged delta. This is
  the sampling error from which proteins are in the eval set.
- ``lo_seed``/``hi_seed``: two-level bootstrap, proteins and then the
  4 document seeds within each protein, on per-seed deltas. It adds the
  error from which document realizations were drawn. Its centre is the mean
  per-seed delta, which differs from ``delta`` because R-precision of an
  averaged matrix is not the average of per-seed R-precisions.

No threshold is applied; the deltas are the raw effect sizes.

    uv run python analyze_b1.py --input-dir out/b1
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

DATA = Path(__file__).resolve().parent / "data"
BANDS = ["all", "short", "medium", "long"]
N_BOOT = 10_000
SEED = 341


def load(input_dir: Path, sub: str) -> pd.DataFrame:
    frames = [pd.read_csv(path, dtype={"seed": str}) for path in sorted((input_dir / sub).glob("*.csv"))]
    return pd.concat(frames, ignore_index=True)


def splits(proteins: pd.DataFrame) -> dict[str, np.ndarray]:
    """Boolean masks over ``proteins`` (one row per stem) for each split."""
    terciles = pd.qcut(proteins.n_residues, 3, labels=False)
    edges = proteins.groupby(terciles).n_residues.agg(["min", "max"])
    out = {"all": np.ones(len(proteins), dtype=bool)}
    for eval_set in sorted(proteins.eval_set.unique()):
        out[eval_set] = (proteins.eval_set == eval_set).to_numpy()
    for k, (low, high) in edges.iterrows():
        out[f"length_T{k + 1}_{low}-{high}"] = (terciles == k).to_numpy()
    return out


def bootstrap(per_protein: np.ndarray, rng: np.random.Generator) -> tuple[float, float]:
    """95% CI of the mean over proteins (axis 0)."""
    idx = rng.integers(0, len(per_protein), size=(N_BOOT, len(per_protein)))
    means = per_protein[idx].mean(axis=1)
    return float(np.quantile(means, 0.025)), float(np.quantile(means, 0.975))


def bootstrap_seeds(per_seed: np.ndarray, rng: np.random.Generator) -> tuple[float, float]:
    """95% CI of the mean, resampling proteins (axis 0) then seeds (axis 1) within each."""
    n, s = per_seed.shape
    proteins = rng.integers(0, n, size=(N_BOOT, n))
    seeds = rng.integers(0, s, size=(N_BOOT, n, s))
    means = per_seed[proteins[:, :, None], seeds].mean(axis=(1, 2))
    return float(np.quantile(means, 0.025)), float(np.quantile(means, 0.975))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input-dir", type=Path, required=True)
    args = parser.parse_args()

    scores = load(args.input_dir, "scores")
    scores.to_csv(DATA / "b1_scores.csv", index=False)
    load(args.input_dir, "timings").to_csv(DATA / "b1_timings.csv", index=False)

    mean = scores[scores.seed == "mean"]
    per_seed = scores[scores.seed != "mean"]
    proteins = mean[mean.condition == "clean"][["stem", "eval_set", "n_residues"]].sort_values("stem").reset_index(drop=True)
    stems = proteins.stem.to_numpy()
    seeds = sorted(per_seed.seed.unique(), key=int)
    masks = splits(proteins)

    def matrix(frame: pd.DataFrame, condition: str, band: str) -> np.ndarray:
        rows = frame[frame.condition == condition]
        if frame is mean:
            return rows.set_index("stem").loc[stems, f"rp_{band}"].to_numpy()
        return rows.pivot(index="stem", columns="seed", values=f"rp_{band}").loc[stems, seeds].to_numpy()

    rng = np.random.default_rng(SEED)
    rows = []
    for condition in sorted(c for c in mean.condition.unique() if c != "clean"):
        layer, site = int(condition[1:3]), condition.split("_")[1]
        for band in BANDS:
            delta = matrix(mean, condition, band) - matrix(mean, "clean", band)
            delta_seed = matrix(per_seed, condition, band) - matrix(per_seed, "clean", band)
            clean = matrix(mean, "clean", band)
            for split, mask in masks.items():
                lo, hi = bootstrap(delta[mask], rng)
                lo_seed, hi_seed = bootstrap_seeds(delta_seed[mask], rng)
                rows.append({
                    "condition": condition, "layer": layer, "site": site, "band": band, "split": split,
                    "n_proteins": int(mask.sum()), "clean_rp": clean[mask].mean(), "delta": delta[mask].mean(),
                    "lo": lo, "hi": hi, "delta_seed": delta_seed[mask].mean(), "lo_seed": lo_seed, "hi_seed": hi_seed,
                })
    out = pd.DataFrame(rows)
    out.to_csv(DATA / "b1_deltas.csv", index=False, float_format="%.5g")

    table = out[(out.split == "all")].pivot_table(index=["layer"], columns=["site", "band"], values="delta")
    with pd.option_context("display.width", 200, "display.float_format", "{:+.3f}".format):
        print(table)
    return 0


if __name__ == "__main__":
    sys.exit(main())
