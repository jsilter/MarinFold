# Plan: a PINDER-style split for the AFDB dimer release

Status: phases 0 and 2-pilot complete, 2026-09-05. Phases 1 and 3 onward not run.

Context: [`README.md`](README.md) is the survey that issue #145 asked for. This
plan covers the curation work that followed from it, and the decision recorded in
conversation to follow PINDER's split protocol rather than the sequence-only split
already built.

## What this dataset is for, and why that decides its size

**Primary deliverable: a held-out multimer validation set**, a few thousand
non-redundant complexes usable to validate models trained on contacts-v1
monomers. Not a training corpus. This is the decision that sizes every download
below, so it is stated first rather than left implicit.

The consequence is that the clustering population is the **val-eligible 881,010
complexes**, not all 2,010,800. The other 1,129,790 would be clustered only to
define communities nobody draws from.

One caveat, stated because it is the thing that would change if the goal changed:
PINDER's deleaking step exists to keep test systems away from *train* systems. If
we never train on AFDB dimers there is no train side to deleak against, and the
leakage that does matter (against the contacts-v1 monomers our models have seen)
is handled upstream by the sequence split. Reintroducing a training goal
reintroduces the train-side download and the deleaking pass together.

## What already exists

Five steps have run to completion and their outputs are on disk under
`data/curation/`:

| Artifact | Contents |
|---|---|
| `confident_{homo,hetero}dimers.csv.gz` | 1,930,552 + 80,248 complexes passing `ipSAE >= 0.6 and pDockQ2 >= 0.23`, exact counts from a full stream of both metadata CSVs |
| `val_monomers.fasta` | 41,954 contacts-v1 validation monomer sequences, zero unresolved |
| `subunit_shards/*.fasta` | 1,979,542 dimer subunit sequences, zero unresolved, 767 MB |
| `dimer_split_assignment_id30_cov50.csv.gz` | per-complex train/val on subunit homology to the validation monomers, 881,010 val-eligible |
| `fetch_pilot_summary.json`, `fetch_timings.csv` | the phase 2 rate pilot, 8,012 structures |

The confident set is already a filtered slice: the release is 29,025,020 models
(21,430,663 homodimers + 7,594,357 heterodimers) and the gate keeps 6.9% of it.

`dimer_split_assignment_id30_cov50.csv.gz` is **not** superseded by this plan. It
answers a question PINDER's protocol does not: whether a dimer's chains are
homologous to the monomers our contacts-v1 models were validated on. PINDER
controls leakage *within* a dimer set; this controls leakage *between* the dimer
set and our existing monomer validation split. The two compose, and phase 5
combines them.

## Why the sequence-only split is not enough

Two complexes can share no chain above 30% identity and still present the same
binding geometry, and the interface is what a multimer model learns. The 45-model
smoke test in `README.md` already produced a case: a *Merluccius polli* SNX4
homodimer and a *Glycine max* heterodimer clustered at multimer TM 0.62 with no
sequence relationship between them.

The sequence rule also treats heterodimers unfairly, and that is an artifact
rather than a real difference in risk. It routes a complex to val if *either*
subunit hits, so a heterodimer gets two chances and a homodimer one. At 30%
identity that produces 69.9% of heterodimers in val against 42.7% of homodimers.
Assigning whole interface clusters to splits removes the asymmetry, because the
decision is made once per cluster rather than once per chain.

## What PINDER does

Read on 2026-08-24 from `github.com/pinder-org/pinder` at `--depth 1`, files
`src/pinder-data/pinder/data/{config.py,foldseek_utils.py,get_clusters.py,graph_utils.py}`.
The preprint's full text is not open in Europe PMC (record PPR884594, `inEPMC=N`,
`OA=N`), so the code is the source. Where the code and the abstract disagree, the
code wins.

**PINDER does not use `easy-multimercluster`.** An earlier draft of this plan was
written around it and was wrong. What PINDER actually does:

1. Runs Foldseek **`easy-search`, all against all, on individual chains**, not on
   complexes (`foldseek_utils.py:206`).
2. Builds a graph whose nodes are chains and whose edges are alignments passing
   the filters in the parameter table (`graph_utils.py:25`).
3. Runs **asynchronous label-propagation community detection** on that graph to
   get monomer communities (`get_clusters.py:298`, `cluster_from_graph`).
4. Defines a dimer's interface cluster as the **sorted pair of its two chains'
   community IDs**: `cluster_{min(cR, cL)}_{max(cR, cL)}`
   (`get_clusters.py:111`, `sysid_to_cluster_string`).

So "interface cluster" means "this pair of monomer structural communities", and
the interface enters through a filter (both sides must have at least 7 interface
residues) rather than through a multimer alignment. This is cheaper than
`easy-multimercluster` and it scatters into independent jobs.

Their published outcome: 2,319,564 dimeric systems, split into 1,560,682 training
dimers from 42,220 clusters, **1,958 validation representatives** and 1,955 test
systems with interface leakage removed.

Three properties to copy:

1. Cluster on the **interface**, not on chain sequence.
2. Assign **whole clusters** to splits, never individual complexes.
3. Draw validation from **cluster representatives** (`top_n = 1`), so the set is
   small and non-redundant by construction. PINDER's validation set is 0.08% of
   its corpus; our sequence rule produced 43.8%, which is a quarantine pool, not
   a validation set.

## Parameters

| Stage | Parameter | PINDER value | Source |
|---|---|---|---|
| Foldseek search | `-s` sensitivity | **11.0** (Foldseek's own default is 9.5) | `FoldseekConfig` |
| | `-e` evalue | **0.05** (default 0.001, so more distant hits) | `FoldseekConfig` |
| | `--max-seqs` | 1000 | `FoldseekConfig` |
| | `--alignment-type` | **2** (3Di+AA Gotoh-Smith-Waterman, local and fast) | `FoldseekConfig` |
| | score type | **lddt** (not `alntmscore`) | `FoldseekConfig` |
| Graph edges | score band | **> 0.5 and < 1.1** | `GraphConfig.score_threshold`, `upper_threshold` |
| | min alignment length | 10 residues | `GraphConfig` |
| | min interface length | **7 residues on both sides** | `GraphConfig.min_interface_length` |
| | coverage threshold | 0.5 | `GraphConfig` |
| Clustering | algorithm | AsynLPA community detection, seed 40 | `ClusterConfig` |
| | edge threshold for clustering | **0.70** | `foldseek_cluster_edge_threshold` |
| | canonical method | `foldseek_community` | `ClusterConfig` |
| Deleaking | edge threshold for neighbour search | **0.55** (looser than clustering) | `foldseek_edge_threshold` |
| | depth limit | 2 hops | `ClusterConfig.depth_limit` |
| | max node degree assumed leaky | 1000 | `ClusterConfig.max_node_degree` |
| Splitting | representatives per cluster | **`top_n = 1`** | `ClusterConfig.top_n` |
| Batching | chains per Foldseek sub-database | **50,000** | `ScatterConfig.foldseek_db_size` |
| Secondary (MMseqs) | `-s`, `-e`, `--min-seq-id` | 11.0, 0.05, 0.2 | `MMSeqsConfig` |
| | pident score band | 30 to 110 | `GraphConfig.mmseqs_*_threshold` |

`top_n = 1` is how the validation set gets to 1,958 systems: one representative
per cluster, and the rest of each cluster is simply not drawn. Our phase 5 target
of a few thousand complexes follows the same mechanism, not a separate quota.

Two threshold pairs are worth noticing because they are deliberately asymmetric.
Clustering cuts edges below 0.70, but leakage detection searches neighbours down
to 0.55, so a complex can be excluded from the test set for a similarity that was
too weak to have put it in the same cluster. And the search runs at `-e 0.05`
against Foldseek's default `0.001`, which deliberately admits more distant hits
than a normal search would.

## What does not transfer

PINDER's test-set criteria are properties of PDB depositions and have no analogue
in a predicted set: `resolution_thr 3.5`, `method "X-RAY DIFFRACTION"`,
`interface_atom_gaps_4A 0`, `prodigy_label "BIO"`, `oligomeric_count 2`,
`alphafold_cutoff_date "2021-10-01"`. Our equivalent quality filter is the
confidence gate already applied upstream (`ipSAE >= 0.6 and pDockQ2 >= 0.23`),
and our equivalent of the AF2-cutoff holdout does not exist, since every model in
the release was produced by AlphaFold-Multimer v2.3.0 weights.

`min_chain_length 40` and `min_atom_types 3` do transfer and are cheap to apply.

## Phases

### Phase 0. Read the PINDER methods and fix parameters. DONE 2026-08-24

Deliverable is the parameter table above, with divergences marked.

### Phase 1. Sequence-collapse the val-eligible set before downloading anything

New, and it comes first because it is free and it sizes everything after it.

The val-eligible 881,010 complexes contain **858,360 distinct subunit
accessions**, whose sequences are already on disk. Cluster them with
`mmseqs easy-cluster` and count the clusters. Complexes whose chains fall in the
same sequence clusters will almost always land in the same interface community,
so one complex per distinct sequence-cluster pair is enough to define the
structural clustering.

The risk is real and stateable rather than unknown: near-identical chains can
still present different interface geometry (domain swaps, alternative binding
modes), so this merges some things PINDER would split. For a validation set that
errs toward fewer and more diverse representatives, which is the safe direction.
Record how many complexes the collapse discards so the cost stays visible.

Deliverable: the number of complexes phase 2 actually has to fetch. Until it
exists, treat 881,010 (103 GB) as the ceiling and 2,010,800 (235 GB) as the
figure for a training corpus we are not currently building.

### Phase 2. Fetch the structures

**Rate pilot: DONE 2026-09-05**, 8,012 structures via
`curate_fetch_structures.py --sweep`. Results:

| workers | files/s | speedup | 429s | retries |
|---|---|---|---|---|
| 1 | 1.44 | 1.00x | 0 | 0 |
| 2 | 2.85 | 1.98x | 0 | 0 |
| 4 | 5.96 | 4.14x | 0 | 0 |
| 8 | 11.54 | 8.01x | 0 | 0 |

Throughput is linear in workers with no throttling of any kind. The mechanism:
`alphafold.ebi.ac.uk/files/` is served out of Google Cloud Storage, not EBI's own
web tier (`server: UploadServer`, `x-guploader-uploadid`, `x-goog-storage-class`
headers, a Google load balancer at 34.149.152.8), with
`cache-control: public,max-age=86400`.

**We cap at 8 workers anyway.** No rate limit is published: EMBL-EBI's terms of
use state only that "any attempt to use EMBL-EBI Data Resources and Tools to a
level that prevents, or is likely to prevent, EMBL-EBI providing services to
others, will result in the user being blocked", and their `robots.txt` disallows
`/api` and `/search/sequence` while allowing `/files/`. The only quantitative EBI
figure anywhere is 10 req/s per IP, forum-sourced from Europe PMC staff for a
different service. 11.5 files/s sits at that order and is demonstrated safe. If
we ever need more, the honest route is to email `afdbhelp@ebi.ac.uk` and ask for
the GCS bucket name, which would make this an in-network copy.

**Measured sizes** over 8,012 files: 117.0 kB on the wire, 508.1 kB decoded,
4.34x. Per population:

| population | complexes | wire | decoded | at 11.5 files/s |
|---|---|---|---|---|
| val-eligible | 881,010 | 103 GB | 448 GB | 21 h |
| all confident | 2,010,800 | 235 GB | 1,022 GB | 49 h |

**The release has gaps.** 28 of 2,000 files in one batch returned 404, all
homodimers in the ID band `AF-...74043945` to `...74276307`. Re-requesting all 28
serially returned 404 again, so they are missing rather than load-shed, and a
random 400-model sample across the whole confident set was 400/400 present. The
fetcher records a 404 and continues; **the list of missing models ships as a
published artifact**, because the dataset is defined by what actually downloaded.

Shard the ID list by `index % n_shards` over a sorted list, not in blocks, for the
same reason as the sequence fetch. Record per-file timings to CSV as `AGENTS.md`
requires.

### Phase 2a. Storage and region

Local disk has 58 GB free, so this needs cloud storage regardless of population.
Decision recorded in conversation: **AWS**, which exp91 already used (EC2 workers
reading and writing S3).

Costs at current us-east-1 pricing, first month:

| line item | val-eligible (881,010 / 103 GB) | all confident (2,010,800 / 235 GB) |
|---|---|---|
| EC2, `c7i.large`, 21 h vs 49 h | $2 on demand | $4 |
| ingress to AWS | $0 | $0 |
| S3 PUT at $0.005/1,000 | $4.41 | $10.05 |
| S3 Standard at $0.023/GB-mo | $2.37/mo | $5.41/mo |
| later egress to the HF bucket, $0.09/GB after 100 GB free | $0.27 | $12.15 |
| **total** | **~$9** | **~$32** |

Two things follow. The PUT charge is a third of the total at 2 M objects, so
**pack structures into tar or parquet shards** rather than writing one key per
structure. And the eventual S3-to-HF move is over the `AGENTS.md` 10 GB
cross-region threshold and needs explicit sign-off; the val-eligible population
keeps that leg under the 100 GB free-egress allowance, which the full set does
not.

Per [[ec2-worker-no-boto3]]: on ephemeral EC2 workers, read and write S3 through
`pyarrow.fs.S3FileSystem` with the instance role. Do not `pip install boto3`; it
wedges the box.

### Phase 3. Measure Foldseek scaling on one sub-database pair

The unit is the distinct chain, not the complex. At PINDER's 50,000-chain
sub-database size, 858,360 val-side chains is ~18 sub-databases and ~160 ordered
pairs; the full 1,979,542 would be ~40 and ~1,600. Each pair is an independent
job, so this fans out.

Measure wall time and peak RSS for **one 50k-by-50k pair** at `-s 11.0`,
`--alignment-type 2`, lddt scoring, then multiply. One pair needs only 100,000
structures (12 GB, 2.5 hours of fetching), so this happens **before** the full
download, not after. The only Foldseek numbers we have are from 45 structures and
nothing about this scale follows from them.

Deliverable: a scaling table in this file and a fan-out width.

### Phase 4. Graph, communities, cluster IDs

Build the chain graph with the edge filters in the parameter table, cut edges
below 0.70, run AsynLPA at seed 40, then label each complex `cluster_{min}_{max}`
from its two chains' community IDs. Homodimers get `cluster_c_c`, which falls out
of the construction rather than needing a special case.

**Interface residues are a separate computation.** The 7-residues-both-sides
filter needs per-complex contact residues, which Foldseek does not report. That is
a pass over the downloaded structures (contact residues within a distance cutoff
on each side) and it has to happen before the graph is built. It is not in the
Foldseek cost and was missing from earlier drafts of this plan.

Guard: `smoke_test_foldseek_multimer.py` carries `assert_all_chains_ingested`,
which catches Foldseek silently ingesting a fraction of the input and exiting 0.
That check must run on every sub-database build.

Deliverable: `cluster_id` per `model_entity_id`, plus the cluster-size
distribution.

### Phase 5. Draw the validation set from cluster representatives

For each interface cluster:

1. Confirm every member is flagged val-eligible by
   `dimer_split_assignment_id30_cov50.csv.gz`. A cluster mixing val-eligible and
   train-side members is a leakage path and goes to val whole, as in earlier
   drafts; with a val-only download it should not arise, and if it does, that is
   a finding worth reporting rather than silently resolving.
2. Draw **one representative per cluster** (`top_n = 1`), targeting a few thousand.

Report, so the cost of the criterion stays visible: how many clusters and how many
complexes land in each split, and how many complexes are discarded (in a
val-eligible cluster but not drawn as a representative).

### Phase 6. Publish

Structures and per-complex tables to the HF bucket under
`data/afdb-nvda-dimers/`, following the `AGENTS.md` prefix rule. Small artifacts
(cluster-size distribution, split counts, scaling table, the 404 manifest) stay in
`data/` in git. A dataset README stating the gate, the split rule, the thresholds
and the exact counts ships with it.

## Costs and unknowns

| Item | Status |
|---|---|
| EBI aggregate rate limit | **measured**: linear to 8 workers, 11.5 files/s, no throttling; GCS-backed |
| Structure transfer | **measured**: 117.0 kB wire / 508.1 kB decoded per model over 8,012 files |
| Missing models | **measured**: clustered 404s, 0/400 in a whole-set sample; manifest ships with the data |
| Storage and AWS cost | **estimated**: ~$9 val-eligible, ~$32 full set, first month |
| PINDER's exact parameters | **read**, table above |
| Sequence collapse factor | **unknown**, phase 1 measures it, free |
| Foldseek cost per 50k pair | **unknown**, phase 3 measures it on 100 k structures |
| Interface-residue computation | **unknown**, a pass over downloaded structures, not yet designed |

The sequence-only split needed no structures. This one does, and that is the
material change: phases 1 through 3 exist to make the structure download as small
as the deliverable allows.
