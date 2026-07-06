# exp91 distillation run — status (2026-07-06)

Status of the ESM Atlas (ESMFold2 Atlas, 1.1B predicted monomers) distillation
run that produces a quality-filtered, novel, de-duplicated subset to augment the
afdb-24M training set. For the *why* (characterization, funnel design) see
[`README.md`](README.md); for the *how* (runbook) see [`PIPELINE.md`](PIPELINE.md).

## TL;DR

The full funnel **completed successfully** on 2026-07-06. We now have, for
**66,759,963 representative proteins**, their **sequences + selection metadata**
in S3. We do **not** yet have their **3D structures** — decoding those is the
remaining `materialize.py` step (multi-TB, see "What remains").

## The run

- Instance: `i-073c8a0b3fa1a592f` (r7i.24xlarge, 96 vCPU / 768 GB, us-west-2),
  launched 2026-07-05 19:50 UTC, pipeline `_DONE` 2026-07-06 10:57 UTC.
- Funnel wall-time ~15.1h; run cost ~$95. Self-terminates after the final upload.
- Third attempt: the first two runs OOM-wedged (see "Failures & fixes"). This one
  used the readahead-capped scan + chunked novelty + per-chunk S3 checkpointing.

### Final funnel (1,095,530,880 Atlas rows in)

| Stage | Result | Time |
|---|---|---|
| scan | 383,988,806 survivors (pLDDT ≥ 0.70, pTM ≥ 0.50, len 60-1000) | 21 min |
| novelty | −220,802,158 (57.5% redundant vs afdb-24M @ 40% id, `-s 4.0`) | ~13.0h (77 chunks) |
| leakage | −41,517 (vs eval set; negligible, as predicted) → kept **163,144,153** | 42 min |
| cluster | 163M → **66,759,963 clusters** (mmseqs linclust @ 40% id) | 34 min |
| select | **66,759,963 reps** (1/cluster: longest, then lowest pLDDT std) | 23 min |

## What we have now (`s3://marinfold-exp91-usw2/exp91/out/`)

Durable — no bucket lifecycle/expiry policy, versioning off; persists until
explicitly deleted. ~$4/month in S3 Standard. **The structures are NOT here** (they
remain as encoded `structure_blob` in the Atlas Lance dataset).

| Object | Size | Contents |
|---|---|---|
| `selected_manifest.csv` | 8.64 GB | **The deliverable.** 66,759,963 rows; columns `cluster_id, protein_hash, seq_len, mean_plddt, ptm, plddt_std, cluster_size`. Authoritative list of selected reps + quality scores. Joins to the FASTA on `protein_hash`. |
| `clu_rep_seq.fasta` | 18.2 GB | Amino-acid **sequence** of every rep, keyed by `protein_hash`. This is the training sequence set. |
| `clu_cluster.tsv` | 10.8 GB | Full cluster membership (rep → all members), if the non-rep neighbors are ever wanted. |
| `cluster_map.csv` | 10.8 GB | Same membership, CSV form. |
| `survivors_0_*.fasta` | ~3.8 GB × 32 | Intermediate: the 384M pre-clustering survivor sequences (per scan shard). Not needed for training; kept for debugging. |
| `pipeline.log` | 4.9 MB | Full run log (per-stage/per-chunk timings, RAM/disk). |
| `_progress_novelty/`, `_progress_leakage/` | — | Per-chunk drop-ID checkpoints (resume state). |
| `_DONE`, `_pipeline_ok` | — | Success markers. |

## Validated facts (useful for the next run)

- The scan readahead cap held: RAM avail bottomed at **360 GB** (the wedged 64-worker
  run hit ~0). See `pipeline.py` `_scan_range` (`batch_readahead=4, fragment_readahead=2`).
- Clustering 163M sequences peaked at only ~45 GB used (~723 GB avail), so **768 GB
  was overkill**. A future rerun can drop to **m7i.24xlarge (96 vCPU / 384 GB)** —
  same speed (novelty is core-bound, not RAM-bound), ~20% cheaper.
- The `_progress_*` → S3 checkpoint round-trip works: a failed run resumes via
  `./create_dataset.sh run` (re-scans ~21 min, rebuilds the target DB, skips done chunks).
- Leakage vs the eval set is ~nil (41.5K of 384M), consistent with the
  metagenomic-vs-recent-PDB reasoning.

## What remains

1. **Decide the final size before materializing.** 66.76M reps came from clustering at
   40% id; the eval-dataset design targeted a 10-100M subset, so this fits, but
   materializing all 66.76M structures is multi-TB. Options to shrink: tighten
   `CLUSTER_ID` (coarser clusters), or sample the manifest. **This is a deliberate
   decision, not a default.**

2. **Materialize structures (`materialize.py`).** Decode each selected rep's
   `structure_blob` (Atlas Lance) → mmCIF → parquet in the `contacts-v1 generate`
   schema (`entry_id` + `cif_content`). Size budget:
   - ~150-180 KB of CIF text per structure × 66.76M ≈ **~10 TB uncompressed**
     (~2-4 TB compressed parquet), plus ~0.7 TB to read the blobs from the Atlas.
   - Run **in-region (us-west-2)** to avoid egress; shardable
     (`--num-shards N --shard-id i`). Write to
     `gs://marin-<region>/...` only after a single bulk copy (respect the >10 GB
     cross-region sign-off rule); the Atlas is in AWS, marin is GCS.

3. **Publish design (gated).** Any HF publish is cross-cloud, multi-TB (needs
   sign-off), and the Atlas is **CC BY-SA 4.0** (ShareAlike) — verify exact license
   terms before republishing a derived dataset.

4. **Housekeeping.** The finish upload pushes the intermediate `survivors_0_*.fasta`
   (~122 GB) every run; a future edit to `run_aws.py`'s finish() excludes could skip
   those and save ~20 min + a few $ per run.

## Failures & fixes (for reproducibility)

- **Run 1** (`i-060bb48dd327f9dce`): novelty stage `mmseqs easy-search` over ~384M
  queries fell into split-prefilter mode, generated 2.3 TB+ scratch, thrashed. Fixed
  by chunking the query set (`--query-chunk-seqs`) with per-chunk scratch deletion.
- **Run 2** (`i-0394e34f2b6368f68`): scan OOM-wedged — 64 workers × uncapped Lance S3
  readahead × the heavy `per_residue_plddt` column marched RAM 772 → 0 GB. Fixed by
  capping scanner readahead + dropping `SCAN_WORKERS` 64 → 32 (commit `b44dc0d`).
- **Run 3** (`i-073c8a0b3fa1a592f`): succeeded.
