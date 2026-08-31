# Plan: a PINDER-style split for the AFDB dimer release

Status: proposed, 2026-08-24. Nothing in phases 2 onward has been run.

Context: [`README.md`](README.md) is the survey that issue #145 asked for. This
plan covers the curation work that followed from it, and the decision recorded in
conversation to follow PINDER's split protocol rather than the sequence-only split
already built.

## What already exists

Four steps have run to completion and their outputs are on disk under
`data/curation/`:

| Artifact | Contents |
|---|---|
| `confident_{homo,hetero}dimers.csv.gz` | 1,930,552 + 80,248 complexes passing `ipSAE >= 0.6 and pDockQ2 >= 0.23`, exact counts from a full stream of both metadata CSVs |
| `val_monomers.fasta` | 41,954 contacts-v1 validation monomer sequences, zero unresolved |
| `subunit_shards/*.fasta` | 1,979,542 dimer subunit sequences, zero unresolved, 767 MB |
| `dimer_split_assignment_id30_cov50.csv.gz` | per-complex train/val on subunit homology to the validation monomers, 881,010 val-eligible |

That last file is **not** superseded by this plan. It answers a question PINDER's
protocol does not: whether a dimer's chains are homologous to the monomers our
contacts-v1 models were validated on. PINDER controls leakage *within* a dimer
set; this controls leakage *between* the dimer set and our existing monomer
validation split. The two compose, and phase 5 combines them.

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

Recorded in `README.md` from the PINDER preprint (doi 10.1101/2024.07.17.603980):
2,319,564 dimeric systems clustered by interface structure, split at the cluster
level into 1,560,682 training dimers from 42,220 clusters, **1,958 validation
representatives** and 1,955 test systems with interface leakage removed.

Three properties to copy:

1. Cluster on the **interface**, not on chain sequence.
2. Assign **whole clusters** to splits, never individual complexes.
3. Draw validation and test from **cluster representatives**, so they are small
   and non-redundant by construction. PINDER's validation set is 0.08% of its
   corpus; our sequence rule produced 43.8%, which is a quarantine pool, not a
   validation set.

**Open item:** PINDER's exact clustering parameters (Foldseek settings, whether
iAlign refinement is used, how the deleaking graph is built) are not recorded in
`README.md` and I have not read the methods section. Phase 0 is to read it, so the
parameters below are replaced by theirs where they differ.

## Phases

### Phase 0. Read the PINDER methods and fix parameters

Cheap and blocking. Everything downstream depends on the clustering thresholds and
the split-assignment procedure. Deliverable: the parameter table appended to this
file, with anything we deliberately diverge on marked and justified.

### Phase 1. Choose a region, and accept that this needs cloud storage

The structures are **~230 GB gzipped on the wire and ~950 GB decoded**, measured
at 466 kB of mmCIF and 115 kB gzipped per model over ten samples. Local disk has
58 GB free, so this cannot run on the workstation.

The source is EMBL-EBI in Hinxton, UK, so the transatlantic leg is unavoidable if
compute lands in a US region. `europe-west4` is nearest to EBI and would make the
fetch the cheapest leg; `AGENTS.md` warns against `europe-west4` because exp53
spilled workers there and got straggler tails, but that was a job whose *data* was
in us-central1, which is the opposite of this case. Whichever region is chosen,
everything downstream stays in it: a 230 GB cross-region copy needs explicit
human sign-off under `AGENTS.md`, and doing it twice would be a real cost event.

Deliverable: a region, a bucket prefix
(`gs://marin-<region>/protein-structure/MarinFold/exp145-afdb-dimers/`), and a
pinned worker zone.

### Phase 2. Fetch 2,010,800 structures

Per-file GETs against `https://alphafold.ebi.ac.uk/files/AF-<id>-model_v1.cif`
with `Accept-Encoding: gzip`. The FTP tar shards are no help: the confident slice
is 9% of the release spread across 9,892 + 8,209 shards, so pulling it via tars
means pulling all 48.8 TB.

Measured inputs: 0.2 to 0.9 s per file serially, and 4.1x compression over ten
models. **Unmeasured and pilot-blocking:** EBI's aggregate rate limit. UniProt
throttled us per source IP during the sequence fetch, where four concurrent
workers on one host produced *lower* aggregate throughput than one; if the EBI
file endpoint behaves the same way, concurrency has to come from distinct source
addresses rather than from threads.

Pilot first: 20,000 structures, one worker, then four workers on one host, then
four workers on four hosts. Report throughput for each before sizing the fan-out.

Shard the ID list by `index % n_shards` over a sorted list, not in blocks, for the
same reason as the sequence fetch. Write directly to object storage through
fsspec; never stage 950 GB on a worker's local disk. Record per-file timings to a
CSV as `AGENTS.md` requires.

### Phase 3. Measure Foldseek-Multimer scaling before committing to a full run

The only Foldseek-Multimer numbers we have are from 45 structures: 90 chains,
3.6 s to cluster, 3.3 s to search. Nothing about 2 M complexes follows from that.

Run `easy-multimercluster` at 10 k, 50 k and 200 k complexes, recording wall time,
peak RSS and scratch-disk high-water mark at each point, then fit. The failure to
plan for is memory on the all-versus-all prefilter over ~4 M chains, not CPU time.
If the fit says the full run does not fit on one machine, the fallback is to
cluster in taxon-blocked batches and then merge representatives, which is weaker
and should only be reached for with the measurement in hand.

Deliverable: a scaling table in this file and a go/no-go on a single-machine run.

### Phase 4. Cluster the full set

`easy-multimercluster` over all 2,010,800 confident complexes, homodimers and
heterodimers **in one job**. Splitting them into two jobs is the specific mistake
the smoke test found, since the cross-set cluster it produced would be invisible.

Starting thresholds (subject to phase 0): `--multimer-tm-threshold 0.5`,
`--chain-tm-threshold 0.0`, `--interface-lddt-threshold 0.0`.

Guard: `smoke_test_foldseek_multimer.py` already carries
`assert_all_chains_ingested`, which catches Foldseek silently ingesting a fraction
of the input and exiting 0. That check must run here too, at 4,021,600 chains.

Deliverable: `cluster_id` per `model_entity_id`, plus the cluster-size
distribution.

### Phase 5. Assign splits at the cluster level

For each interface cluster:

1. If **any** member is flagged by `dimer_split_assignment_id30_cov50.csv.gz` as
   homologous to a contacts-v1 validation monomer, the whole cluster goes to val.
   This is what removes the homo/heterodimer asymmetry: the decision is made once
   per cluster, not once per chain.
2. Otherwise assign the cluster to train, holding back a fraction for test.
3. Draw the actual validation and test sets from **cluster representatives**,
   targeting a few thousand each rather than the full eligible pool.

Report, so the cost of the criterion stays visible: how many clusters and how many
complexes land in each split, how many clusters were forced to val by rule 1, and
how many complexes are discarded (in a val-eligible cluster but not drawn as a
representative).

### Phase 6. Publish

Structures and per-complex tables to the HF bucket under
`data/afdb-nvda-dimers/`, following the `AGENTS.md` prefix rule. Small artifacts
(cluster-size distribution, split counts, scaling table) stay in `data/` in git.
A dataset README stating the gate, the split rule, the thresholds and the exact
counts ships with it.

## Costs and unknowns, stated up front

| Item | Status |
|---|---|
| Structure transfer | ~230 GB, measured per-model, not yet measured in bulk |
| Storage | ~230 GB compressed, ~950 GB if kept decoded |
| EBI aggregate rate limit | **unknown**, pilot-blocking |
| Foldseek-Multimer cost at 2 M | **unknown**, phase 3 measures it |
| PINDER's exact parameters | **unread**, phase 0 |
| Local disk | 58 GB free, insufficient; cloud storage required |

The sequence-only split needed no structures. This one does, and that is the
material change: phases 1 and 2 are the price of clustering on interfaces rather
than on chains.

---

# Phase 0 result: PINDER's actual parameters

Read on 2026-08-24 from `github.com/pinder-org/pinder` at `--depth 1`, files
`src/pinder-data/pinder/data/{config.py,foldseek_utils.py,get_clusters.py,graph_utils.py}`.
The preprint's full text is not open in Europe PMC (record PPR884594,
`inEPMC=N`, `OA=N`), so the code is the source. Where the code and the abstract
disagree, the code wins.

## The finding that changes the plan

**PINDER does not use `easy-multimercluster`.** Phase 4 above was written around
it and is wrong.

What PINDER actually does:

1. Runs Foldseek **`easy-search`, all against all, on individual chains**, not on
   complexes (`foldseek_utils.py:206`).
2. Builds a graph whose nodes are chains and whose edges are alignments passing
   the filters below (`graph_utils.py:25`).
3. Runs **asynchronous label-propagation community detection** on that graph to
   get monomer communities (`get_clusters.py:298`, `cluster_from_graph`).
4. Defines a dimer's interface cluster as the **sorted pair of its two chains'
   community IDs**: `cluster_{min(cR, cL)}_{max(cR, cL)}`
   (`get_clusters.py:111`, `sysid_to_cluster_string`).

So "interface cluster" means "this pair of monomer structural communities", and
the interface enters through a filter (both sides must have at least 7 interface
residues) rather than through a multimer alignment. This is cheaper than
`easy-multimercluster` and it scatters, which resolves the phase 3 scaling worry
in a way the old plan did not anticipate.

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
per cluster, and the rest of each cluster is simply not drawn.

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

## Revised phases 3 and 4

**Phase 3 (revised): all-against-all Foldseek `easy-search` over distinct chains.**

The unit is the distinct chain, not the complex. Our 2,010,800 complexes contain
**1,979,543 distinct subunit accessions** (a homodimer's two chains are the same
protein), which is the same set we already fetched sequences for. At PINDER's
50,000-chain sub-database size that is ~40 sub-databases and ~1,600 ordered
sub-database pairs, each an independent job. This is embarrassingly parallel and
belongs in a fan-out, which is a better shape than the single large
`easy-multimercluster` run phase 4 originally assumed.

Still to measure: wall time and peak RSS for one 50k-by-50k pair at `-s 11.0`,
`--alignment-type 2`. Multiply by ~1,600 for the total, then decide the fan-out
width. Measure this before fetching all 2 M structures, since one pair needs only
100,000 structures and the answer might change the whole approach.

**Phase 4 (revised): graph, communities, cluster IDs.**

Build the chain graph with the edge filters in the table, cut edges below 0.70,
run AsynLPA at seed 40, then label each complex `cluster_{min}_{max}` from its two
chains' community IDs. Homodimers get `cluster_c_c`, which falls out of the
construction rather than needing a special case.

The interface-length filter needs interface residues per complex, which Foldseek
does not give us. That is a separate computation over the downloaded structures
(contact residues within a distance cutoff on each side), and it has to happen
before the graph is built.
