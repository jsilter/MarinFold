---
marinfold_experiment:
  issue: 145
  title: 'Explore multimer data'
  kind: data
  branch: main
---

# Explore multimer data

**Issue:** [#145](https://github.com/Open-Athena/MarinFold/issues/145) · **Kind:** `data` · **Branch:** `main`

## Question

Our training corpora (AFDB, ESM Atlas) are monomer-only. Being able to operate on
multimers would be valuable. **How much multimer content actually exists in the
databases we already use, and in the other databases we could reach?**

This is a survey, not an experiment. Nothing is trained, nothing is generated.
The numbers below were read out of the databases' own APIs, FTP listings, and
object stores on **2026-07-31**; where a number is an estimate rather than an
exact count, it says so.

## Headline

**The AlphaFold DB acquired a multimer corpus this year and it is enormous.** A
four-way EMBL-EBI / Google DeepMind / NVIDIA / Seoul National University release
(announced March 2026, bulk data landing on the FTP site in June and July) holds
roughly **21 M predicted homodimers and 7.6 M predicted heterodimers (about 29 M
dimers and 48.8 TB)**. That is more predicted complexes than everything else in
this report combined, by roughly two orders of magnitude. Every model is scored
(ipTM, ipSAE, pDockQ, pDockQ2, LIS, clash counts) and cross-referenced to UniProt
on both chains, so it can be filtered to a confident subset without re-running
anything. Under the release's own quality gate that subset is **1.96 M
homodimers and 70 k heterodimers**; it is served from the AFDB API and web UI,
and the rest is FTP bulk download.

Everything else is much as expected:

| Source | Multimer content | Non-redundant scale | Access |
|---|---|---|---|
| AFDB main release (v6) | **none**; 241,070,489 single-chain models | n/a | FTP + GCS, public |
| AFDB complexes (`collaborations/nvda/`) | **~29 M predicted dimers** (21.5 M homo, 7.6 M hetero) | 1.96 M homo, 70 k hetero pass the release's quality gate | API + web UI for the confident slice; 48.8 TB FTP bulk |
| ESM Atlas | **none**; 1,095,530,880 + 6,600,755 single-chain models | n/a | S3, public |
| PDB (experimental ground truth) | 200,372 multi-chain protein assemblies | **~40,000** distinct interfaces (three methods agree) | RCSB, public |
| PINDER / PPIRef (PDB-derived) | 2,319,564 dimers / 322,454 interfaces | 42,220 clusters / 45,553 interfaces | public, CC |
| RCSB computed structure models | 2,063 multi-chain (of 1,062,058) | n/a | RCSB |
| ModelArchive | ~11,900 confirmed multi-chain (of 625,966) | concentrated in ~25 deposits | per-deposit ZIP, CC BY-**SA** |
| BFMD (Foldseek) | 297,570 aggregated multimer predictions | 51,757 representatives | Foldseek `databases` |

## AlphaFold DB

### The main release is strictly monomeric

The current release (v6, changelog dated 2025-09-15, built on UniProt 2025_03)
holds **241,070,489 structures**, up from 214,686,924 in v4. Every one is a single
chain: AFDB is keyed on UniProt accession, one accession per file, and the
changelog names the model "AlphaFold Monomer v2.0".

Verified directly rather than assumed. Downloading
`AF-P68871-F1-model_v6.cif` (haemoglobin beta, an obligate component of the α₂β₂
tetramer, so a case where a multimeric model would be biologically warranted)
gives:

```
_entity.id                       1
_entity.pdbx_number_of_molecules 1
_struct_asym.entity_id 1
_struct_asym.id        A
```

One entity, one copy, one chain. Same for P00918, P0AEX9, P04637.

Bulk access is worth noting because it is not what the FTP listing suggests:
`download_metadata.json` describes only **48 curated archives totalling 1,656,090
structures / 0.15 TB** (Swiss-Prot 1,100,244; 24 reference proteomes 328,602;
global-health set 227,244). That 1.66 M is essentially the
[`afdb-1.6M`](https://huggingface.co/datasets/timodonnell/afdb-1.6M) corpus we
already train on. The remaining ~239 M are reachable only through the flat GCS
bucket `gs://public-datasets-deepmind-alphafold-v4` or per-accession API calls.

### The multimer content: `collaborations/nvda/`

The FTP root has a `collaborations/` area holding five partner datasets (`atbc`,
`bfvd`, `ntdx`, `nvda`, `vr3d`. Four are monomer collections. **`nvda/` is not.**

```
https://ftp.ebi.ac.uk/pub/databases/alphafold/collaborations/nvda/
├── homodimer_metadata.csv       6.00 GB
├── heterodimer_metadata.csv     2.45 GB
├── homodimers/                  9,892 tar shards, 35.2 TB
├── heterodimers/                8,209 tar shards, 13.6 TB
└── msas/                        three sets: 1.2 M, 8.8 M, 13.4 M
```

The directories are new: homodimers landed 2026-06-11, heterodimers 2026-07-08,
and the heterodimer metadata was refreshed 2026-07-17.

Row counts are **estimates**, not exact. The metadata CSVs are too large to
download for a line count, so `profile_confidence.py` reads 200 windows of 200 kB
at evenly spaced byte offsets in each and divides the file size by the mean row
length measured in that sample. Row lengths are tight (279.0 B/row for
homodimers, 323.1 for heterodimers), which puts the counting error near 1%:

| | shards | model data | metadata | rows sampled | estimated rows |
|---|---|---|---|---|---|
| homodimers | 9,892 | 35.2 TB | 6.00 GB | 143,142 | **21,489,407** |
| heterodimers | 8,209 | 13.6 TB | 2.45 GB | 123,597 | **7,595,156** |

Every model carries interface confidence scores (`ipTM`, `ipSAE` in both chain
orders, `pDockQ`, `pDockQ2` in both orders, `LIS`, `numberOfInteractions`,
`N_clash_backbone`, `N_clash_heavyAtom`), plus UniProt accessions, taxon IDs and
gene names for both partners.

#### How much of it survives a confidence filter

Measured with `profile_confidence.py`: 200 windows of 200 kB range-read from each
metadata CSV, giving **143,142 homodimer rows (0.67% of the file) and 123,597
heterodimer rows (1.6%)**. Every threshold, with per-window spread, is in
`data/nvda_confidence_gates.csv`.

**One gate applies to both sets, and it is the authors' own.** The heterodimer
table ships a `passes_quality_threshold` column defined by
`quality_ipsae_threshold = 0.6` and `quality_pdockq2_threshold = 0.23`; the
homodimer table ships no verdict column at all. Recomputing that criterion from
the raw score columns reproduces the declared verdict on **123,597 of 123,597
heterodimer rows, zero disagreements in either direction**, so applying the same
recomputation to the homodimers measures them by the release's own standard
rather than by a proxy. That is the `authors' gate` row below.

Cells are `pass rate · estimated complexes`, against 21.5 M homodimers and 7.6 M
heterodimers.

**Homodimers**

| threshold | ipTM | pDockQ | pDockQ2 | ipSAE (max) | ipSAE (both) | LIS |
|---|---|---|---|---|---|---|
| 0.23 | 39.12% · 8.41 M | 32.06% · 6.89 M | 11.36% · 2.44 M | 22.22% · 4.77 M | 22.05% · 4.74 M | 23.17% · 4.98 M |
| 0.3 | 29.54% · 6.35 M | 26.87% · 5.77 M | 10.13% · 2.18 M | 20.63% · 4.43 M | 20.49% · 4.40 M | 19.48% · 4.19 M |
| 0.4 | 22.19% · 4.77 M | 21.23% · 4.56 M | 8.61% · 1.85 M | 18.26% · 3.92 M | 18.09% · 3.89 M | 14.20% · 3.05 M |
| 0.5 | 17.00% · 3.65 M | 16.01% · 3.44 M | 7.27% · 1.56 M | 15.29% · 3.29 M | 15.14% · 3.25 M | 9.04% · 1.94 M |
| 0.6 | 12.81% · 2.75 M | 10.45% · 2.25 M | 6.09% · 1.31 M | 12.29% · 2.64 M | 12.18% · 2.62 M | 4.29% · 921 k |
| 0.7 | 9.21% · 1.98 M | 3.59% · 771 k | 4.84% · 1.04 M | 9.03% · 1.94 M | 8.94% · 1.92 M | 0.65% · 140 k |
| 0.75 | 7.58% · 1.63 M | none | 4.09% · 879 k | 7.35% · 1.58 M | 7.28% · 1.56 M | 0.14% · 29 k |
| 0.8 | 5.90% · 1.27 M | none | 3.30% · 709 k | 5.68% · 1.22 M | 5.61% · 1.21 M | 0.01% · 2 k |
| 0.9 | 2.69% · 577 k | none | 1.02% · 220 k | 1.68% · 361 k | 1.65% · 356 k | none |

**Heterodimers**

| threshold | ipTM | pDockQ | pDockQ2 | ipSAE (max) | ipSAE (both) | LIS |
|---|---|---|---|---|---|---|
| 0.23 | 33.99% · 2.58 M | 17.32% · 1.32 M | 1.55% · 118 k | 5.37% · 408 k | 3.50% · 266 k | 6.23% · 473 k |
| 0.3 | 18.52% · 1.41 M | 11.33% · 861 k | 1.22% · 93 k | 4.51% · 343 k | 2.92% · 222 k | 4.01% · 305 k |
| 0.4 | 9.85% · 748 k | 6.30% · 479 k | 0.88% · 67 k | 3.46% · 263 k | 2.19% · 166 k | 2.04% · 155 k |
| 0.5 | 5.76% · 437 k | 3.21% · 244 k | 0.60% · 46 k | 2.56% · 194 k | 1.54% · 117 k | 0.83% · 63 k |
| 0.6 | 3.35% · 255 k | 1.33% · 101 k | 0.39% · 30 k | 1.75% · 133 k | 0.99% · 75 k | 0.26% · 19 k |
| 0.7 | 1.76% · 134 k | 0.19% · 14 k | 0.23% · 18 k | 1.01% · 77 k | 0.48% · 36 k | 0.05% · 4 k |
| 0.75 | 1.25% · 95 k | none | 0.16% · 12 k | 0.63% · 48 k | 0.26% · 19 k | 246 |
| 0.8 | 0.75% · 57 k | none | 0.09% · 7 k | 0.34% · 26 k | 0.11% · 8 k | none |
| 0.9 | 0.09% · 7 k | none | 0.02% · 1 k | 0.01% · 860 | 307 | none |

`ipSAE (max)` takes the better of the two chain orders, which is what the
release's own `max_ipSAE` column does; `ipSAE (both)` requires both. LIS at its
own published 0.203 cutoff passes 24.67% of homodimers (5.30 M) and 7.37% of
heterodimers (560 k). The two combined gates:

| Gate | homodimers | heterodimers |
|---|---|---|
| **authors' gate** (ipSAE ≥ 0.6 and pDockQ2 ≥ 0.23) | 9.13% · **1.96 M** | 0.92% · **70 k** |
| declared `passes_quality_threshold` | (column absent) | 0.92% · 70 k |
| ipTM ≥ 0.8 and zero backbone clashes | 5.63% · 1.21 M | 0.52% · 39 k |

Five things in these tables matter for planning:

- **Heterodimers are harder by roughly an order of magnitude at every threshold.**
  ipTM ≥ 0.8 passes 5.90% of homodimers and 0.75% of heterodimers; the authors'
  gate passes 9.13% and 0.92%. This is expected rather than a defect: a homodimer
  interface is constrained by symmetry between two copies of one fold, while an
  arbitrary heterodimer pairing may not interact at all, and with an ipTM
  early-stopping threshold of 0.1 the pipeline deliberately keeps the speculative
  pairings instead of discarding them.
- **The metrics disagree with each other by up to 10x.** At threshold 0.5 the
  heterodimer yield is 437 k by ipTM, 244 k by pDockQ, 194 k by ipSAE, and 46 k by
  pDockQ2. Which score you filter on is a bigger decision than where you put the
  cut, and pDockQ2 is the strictest everywhere.
- **A pDockQ gate above 0.7 selects nothing, for a reason that is not about these
  structures.** pDockQ is a fitted sigmoid whose upper asymptote sits below 0.75;
  the largest value anywhere in the two samples is 0.742 (homodimer) and 0.737
  (heterodimer). Anyone porting a "pDockQ ≥ 0.8" rule from elsewhere will silently
  get an empty set.
- **Requiring ipSAE in both chain orders is free on homodimers and halves the
  heterodimers.** ipSAE is directional. On homodimers max and min are nearly the
  same (15.29% against 15.14% at 0.5) because the complex is symmetric; on
  heterodimers at 0.75 the permissive reading gives 48 k and the strict one 19 k.
- **The scores are heavily zero-inflated, so percentiles mislead.** 59.4% of
  homodimers and 76.7% of heterodimers score ipSAE exactly 0. Heterodimer ipSAE
  has a median of 0 and a p90 of 0.020; its p99 is 0.702. Full distributions in
  `data/nvda_confidence_quantiles.csv`:

| | ipTM | pDockQ | pDockQ2 | ipSAE (max) | LIS |
|---|---|---|---|---|---|
| homodimer median / p90 / max | 0.19 / 0.67 / 0.97 | 0.097 / 0.607 / 0.742 | 0.010 / 0.309 / 0.957 | 0.00 / 0.671 / 0.950 | 0.00 / 0.480 / 0.828 |
| heterodimer median / p90 / max | 0.194 / 0.398 / 0.964 | 0.085 / 0.322 / 0.737 | 0.010 / 0.018 / 0.941 | 0.00 / 0.020 / 0.933 | 0.00 / 0.156 / 0.767 |
| fraction scoring exactly 0 (homo / hetero) | 0% / 0% | 1.2% / 7.3% | 1.2% / 7.3% | 59.4% / 76.7% | 51.6% / 67.4% |

#### Why these totals are order-of-magnitude, not three significant figures

Rows are blocked by organism, so a byte window is a cluster sample rather than
700 independent draws. The per-window record in `data/nvda_confidence_windows.csv`
shows how severe that is: in **176 of the 200 homodimer windows a single organism
accounts for more than 90% of the rows** (1,200 organisms across the whole
sample), against 50 of 200 windows for heterodimers.

So the pooled rates carry real uncertainty, much more for homodimers:

| | pooled | window median | window range |
|---|---|---|---|
| homodimer, authors' gate | 9.13% | 7.10% | 0.1% to 98.0% |
| homodimer, ipTM ≥ 0.8 | 5.90% | 3.84% | 0.0% to 98.0% |
| heterodimer, authors' gate | 0.92% | 0.67% | 0.0% to 4.7% |
| heterodimer, ipTM ≥ 0.8 | 0.75% | 0.63% | 0.0% to 4.1% |

The 98% homodimer window is a real feature of the data rather than a parsing
error. It lands inside a run of one well-predicted enzyme (`hemL1`,
glutamate-1-semialdehyde aminomutase, a genuine obligate homodimer) repeated
across hundreds of bacterial strains at consecutive byte offsets, every copy
scoring ipTM 0.95.

#### Distinct model IDs are not distinct proteins

That `hemL1` run points at something the headline counts hide. Every row does
carry its own UniProt accession, but UniProt assigns one accession per strain, so
a conserved bacterial protein enters the set once per sequenced strain:

| homodimer sample | |
|---|---|
| rows | 143,142 |
| distinct UniProt accessions | 143,142 (100%) |
| distinct gene names | 107,195 (74.9%) |
| distinct (gene, taxon) pairs | 118,157 |
| distinct organisms | 1,200 |
| longest run of one gene name at consecutive offsets | 767 |
| rows sitting in a run of 2 or more identical gene names | 15.8% |

Under the authors' gate, the 13,063 confident homodimers in the sample carry
9,964 distinct gene names, so name-level deduplication alone removes about 24% of
them. Heterodimers are far less repetitive: 123,597 rows give 116,284 distinct
gene pairs (94.1%), and the 1,141 confident ones give 1,079 pairs (94.6%).

Gene name is a loose key (`rplC` in *E. coli* and in *S. aureus* are homologs,
not copies of one protein), so these numbers bound the redundancy from one side
only. The structural count will be **lower** than the gene-distinct count, which
is what the Foldseek-Multimer run below exists to measure.

#### Coverage

A separate pass (10 windows of 2 MB, giving 72,286 homodimer and 61,762
heterodimer rows) shows the two sets are shaped quite differently:

- **Homodimers are one model per accession** (every `uniprotAccession` in the
  sample was unique), so they ran a monomer set through a dimer prediction rather
  than screening pairs. One accession is not one protein, though: see the
  redundancy table above, where 143,142 accessions carry only 107,195 distinct
  gene names.
- **Heterodimers are within-organism interactome screens.** 99.5% of pairs have
  `tax_id_1 == tax_id_2`. In the sample, 61,762 pairs involved 55,046 distinct
  proteins with at most 74 pairs per protein, so the pairing is a sparse
  interaction graph, not all-versus-all.
- **Taxonomic scope is broad and plant-heavy.** The most common heterodimer taxa
  seen were *Glycine max* (25.9% of sampled rows), *Arabidopsis thaliana* (10.3%),
  *Homo sapiens* (8.0%), *Oryza sativa* (7.1%), *Rattus norvegicus* (5.3%), *Zea
  mays*, *Mus musculus*. Because of the organism blocking these shares reflect
  which blocks the sample hit, not true dataset proportions, but they do establish
  that the screen spans dozens of species rather than being human-only.

#### Access: two tiers

The high-confidence subset is a first-class part of the database, served by a
**`/api/complex/{accession}`** endpoint that is separate from the familiar
`/api/prediction/` one (querying a dimer's model entity ID against `/prediction/`
returns `{}`, which is what initially misled me into thinking the release was
FTP-only). Asking for human haemoglobin alpha returns its partners directly:

```
$ curl https://alphafold.ebi.ac.uk/api/complex/P69905
  modelEntityId  AF-0000000211046780        providerId  NVIDIA
  uniprotAccession  [P69905, P68873]        gene  [HBA1, HBB]
  assemblyType  Hetero                      oligomericState  dimer
  globalMetricValue 94.52   ipTM 0.92   ipSAE 0.87   pDockQ 0.57   pDockQ2 0.94   LIS 0.69
```

Records carry stoichiometry, a `complexName`, taxon IDs and all five interface
scores. So the split is: **high-confidence complexes through the web UI and API,
the full low-confidence bulk through FTP only.**

#### What produced them, and under what license

The FTP `README.txt` covers only the directory layout and the score schema, but
the structures themselves are fully self-describing. Pulling the first 8 MB of
`heterodimers/shard_0_batch_0.tar` gives individually zstd-compressed ModelCIF
files (`AF-<16-digit-id>-model_v1.cif.zst`, ~30 kB compressed / ~120 kB raw), and
their `_software` and `_ma_protocol_step` blocks name the pipeline exactly:

| step | tool |
|---|---|
| paired MSA | MMseqs2-GPU release 18, ColabFold-style pairing |
| inference | **AlphaFold-Multimer weights v2.3.0**, run through **OpenFold** with TensorRT + cuEquivariance acceleration |
| settings | no templates, 1 model, max 4 recycles, **ipTM early-stopping threshold 0.1** |
| scoring | ipSAE (ipSAE, pDockQ, pDockQ2, LIS) |
| annotation | PyDSSP secondary structure |

That early-stopping threshold explains the long low-confidence tail: rather than
discarding bad pairings, the pipeline runs them cheaply and keeps them. The
primary citation is listed as *"AlphaFold Database expands to proteome-scale
quaternary structures"*, marked **"To be published"**.

#### Cross-check against the official announcement

EMBL announced this on
[2026-03-16, updated 2026-05-19](https://www.embl.org/news/science-technology/first-complexes-alphafold-database/):
a four-way collaboration between EMBL-EBI, Google DeepMind, NVIDIA and Seoul
National University (Steinegger lab), about 17 million GPU-hours, with candidate
pairs drawn from **20 major studied species plus the WHO priority pathogens list**
(which explains the crop and model-organism heavy taxa I sampled). Their published
counts line up with my sampling closely enough to trust both:

| | announced | measured here |
|---|---|---|
| complexes predicted | 30 M | ~29 M |
| homodimers | 1.7 M high-confidence + 18 M lower | 21.5 M total |
| heterodimers | ~80 k high-confidence + 8.1 M lower | 7.6 M total |
| high-confidence homodimers | 1.7 M | 1.96 M (authors' gate), 1.27 M (ipTM ≥ 0.8) |
| high-confidence heterodimers | ~80 k | 70 k (authors' gate, matching their own column) |

The byte-offset sampling reproduced the official figures to within about 10% on
the totals. On the high-confidence homodimers it brackets the announced 1.7 M
rather than reproducing it: the authors' own gate applied to the homodimer scores
gives 1.96 M and an ipTM ≥ 0.8 cut gives 1.27 M. Since the criterion behind the
announced homodimer figure has not been published, 1.7 M sitting between those
two is the expected result, and it says the announced number was produced by
something close to (but not identical with) the heterodimer gate.

**The license is CC-BY-4.0**, stated in each file's `_pdbx_data_usage` record
alongside the text "AVAILABLE FOR ACADEMIC AND COMMERCIAL PURPOSES, UNDER CC-BY
4.0 LICENCE"), the same terms as the AFDB main release, and more permissive than
the ESM Atlas's CC BY-SA 4.0. Usable for us.

Each model also carries full UniProt cross-references for both chains
(`_ma_target_ref_db_details`: accession, gene name, taxon, aligned range), so the
corpus joins cleanly against anything we already key on UniProt.

### Latent multimer content: what AFDB monomers *should* be

AFDB has no multimers, but a large slice of it models proteins that do not exist
as monomers. Counting UniProt SUBUNIT annotations (queried today via the UniProt
REST API):

| | reviewed (Swiss-Prot) | all entries with an AFDB cross-reference |
|---|---|---|
| total | 575,503 | 99,108,531 |
| has a SUBUNIT comment | 301,706 | 8,357,791 |
| annotated homo-oligomer | **101,703** (17.7%) | **2,925,590** |
| annotated hetero-oligomer | **30,937** (5.4%) | **912,992** |
| annotated monomer | 31,545 | n/a |

(Homo/hetero rows use an explicit union of `homodimer`/`homotrimer`/`homotetramer`/
`homopentamer`/`homohexamer`/`homooctamer`/`homooligomer`/`homomultimer` and the
hetero equivalents, rather than a `homo*` prefix wildcard, which also matches
"homolog".)

So at minimum ~2.9 M AFDB entries model a protein that is annotated as a homo-
oligomer; the monomer prediction is a biologically incomplete model of it. In
Swiss-Prot, where annotation is dense, more than three times as many proteins are
annotated as oligomers than as monomers.

The 99.1 M figure is smaller than AFDB's 241 M because UniProt only cross-
references entries that are still current; AFDB retains accessions that have since
been merged or deleted.

### The other four collaboration datasets

`atbc`, `bfvd`, `ntdx` and `vr3d` are monomer chunk-tar collections in the same
layout as the main release (`models/` plus usually `msas/`), several hundred MB
per chunk. `bfvd` is the Big Fantastic Virus Database and `atbc` is almost
certainly AllTheBacteria; `ntdx` and `vr3d` are not identified in any README on
the FTP site.

## ESM Atlas

**No multimers, and none hiding.** `s3://esm-protein-atlas/` (us-west-2, anonymous
read) holds two fold datasets:

| dataset | rows |
|---|---|
| `v1/folds/folds_1B.lance` | 1,095,530,880 |
| `v1/folds/folds_atlas.lance` | 6,600,755 |

Together 75.7 TB across 110,241 objects. Both schemas carry exactly one `sequence`
and one `structure_blob` per row, with no chain column. ESMFold is a single-chain
model, so this is structural, not incidental.

The interesting wrinkle is that **the on-disk structure format is multimer-capable
and the capability is unused**. Decoding a `structure_blob` (brotli → msgpack-numpy)
gives an AlphaFold3-style record:

```
atom37_positions, atom37_mask, atom37_confidence, confidence,
residue_index, insertion_code, sequence, metadata,
chain_id, entity_id, sym_id, chain_boundaries      <- multimer fields
```

`chain_id`, `entity_id`, `sym_id` and `chain_boundaries` are exactly the fields you
need to represent a complex. To check whether any row uses them I sampled 800 rows
from each dataset at eight evenly spaced offsets (contiguous windows are unbiased
here because both tables are keyed on a content hash). **All 1,600 structures had
exactly one distinct `chain_id`, one `entity_id`, one `sym_id`, and
`len(chain_boundaries) == 1`.** No chain-break characters or poly-glycine linkers in
the sequences either.

The cluster annotations are no help: `v1/clusters/data/representative_proteins.parquet`
(7,723,579 rows) carries Pfam domains, taxonomy, product names and UniRef matches;
no interaction, complex, or stoichiometry field.

### The one multimer-relevant signal: genomic neighbourhood

The Atlas's provenance index is **6,824,676,938 rows** across 16 parquet files, with
schema `(protein_hash, source, accession)`. The accessions turn out to encode
genomic position, which is the classic operon signal for physical interaction.
Sampling 1.2 M provenance rows:

| source | share | accession format | positional? |
|---|---|---|---|
| `spire` | 46.3% | `study\|sample\|k141_<contig>_<orf_index>` | yes: contig + ORF index |
| `mgy` | 35.6% | `MGYP003140358282` | no |
| `uniparc` | 10.2% | `UPI001AF2A499` | no |
| `img_m` | 2.9% | `2809785504` | no |
| `umag_prok` | 2.9% | `Ga0209614_10155761`, `Ga0590176_00001086_160893_161660` | yes: scaffold + coordinates |
| `img_vr` | 1.6% | `Ga0309785_1141497` | yes: scaffold + gene index |
| `uhgg` | 0.5% | `MGYG000238883_01394` | yes: genome + gene number |
| `umag_euk` | 0.14% | `Ga0496077_00007_1026162_1026983` | yes: scaffold + coordinates |

All 555,846 sampled SPIRE accessions parse; **74.1% have ORF index ≥ 2**, meaning
they have at least one upstream neighbour on the same contig, and the largest ORF
index seen was 1894. Roughly half the Atlas is therefore addressable as
gene-neighbourhood pairs, at a scale (hundreds of millions of adjacent ORF pairs)
that nothing else in this report approaches.

To be clear about what this is and is not: the Atlas contains **zero multimer
structures**. What it contains is an enormous set of *candidate interacting pairs*
with a folded monomer already available for each partner. Turning that into
structures would mean running a complex predictor ourselves.

## PDB: the experimental baseline

Measured with `probe_pdb.py` (faceted count queries against the RCSB search API;
kilobytes transferred, no structures downloaded). CSVs in `data/`.

Of 257,179 experimental entries, **359,783 assemblies contain at least one protein
chain**, and of those:

| | assemblies |
|---|---|
| single protein chain | 159,411 |
| **multiple protein chains** | **200,372** (55.7%) |
| ...homomeric (1 distinct entity) | 117,300 |
| ...heteromeric (≥ 2 distinct entities) | 83,072 |
| ...containing nucleic acid too | 12,588 |

Chain-count distribution is dominated by dimers (105,626 two-chain assemblies),
then tetramers (29,705) and trimers (22,423), with a long tail out to 300+ chains.
Point-group symmetry mirrors this: C2 94,330, C1 74,019, D2 14,448, C3 13,318.

By method, cryo-EM is where the multimers are: **92% of EM entries are multi-chain**
(32,751 of 35,608) and 18,942 have six or more chains, versus 58% for X-ray
(117,860 of 202,794) and 14% for NMR.

**Redundancy is the real story.** 477,774 protein entities participate in a
multi-chain assembly, but they collapse hard:

| sequence identity | distinct clusters |
|---|---|
| 100% | 105,440 |
| 95% | 81,979 |
| 90% | 75,050 |
| 70% | 55,313 |
| 50% | 46,147 |
| **30%** | **37,124** |

So the PDB's 200 k multimeric assemblies represent roughly **37 k distinct protein
families**. Growth is steady but not fast: 11,795 multi-chain protein entries were
released in 2025 (6,145 of them heteromeric), and 6,943 so far in 2026.

RCSB also mirrors 1,062,058 computed structure models (AFDB and ModelArchive
deposits). Only **2,063 are multi-chain**, all heteromeric, being the ModelArchive
complex deposits.

## Other databases

### PDB-derived interface datasets, and a useful convergence

Two curated redundancy-reduced interface sets are worth knowing, mostly because
they independently corroborate the scale I measured against RCSB:

- **[PINDER](https://doi.org/10.1101/2024.07.17.603980)** (2024): 2,319,564 dimeric
  PPI systems mined from the PDB, split by structural clustering into a training
  set of **1,560,682 dimers drawn from 42,220 clusters**, a 1,958-representative
  validation set, and 1,955 high-quality test PPIs with interface leakage removed
  (plus a 180-dimer subset clean with respect to AlphaFold-Multimer's training
  data). Notably it already pairs 566,171 of its systems with AFDB monomer
  structures, which is exactly the apo/holo bridge we would otherwise have to
  build.
- **[PPIRef](https://github.com/anton-bushuiev/PPIRef)** (ICLR 2024): PPIRef300K is
  322,454 biophysically valid interfaces; deduplicating with the `iDist` interface
  similarity measure leaves **PPIRef50K at 45,553 non-redundant interfaces**.

Put next to my own RCSB census, three different methods land in the same place:

| Method | Non-redundant multimer units |
|---|---|
| RCSB sequence clustering at 30% identity (this report) | 37,124 |
| PINDER structural interface clustering (training clusters) | 42,220 |
| PPIRef `iDist` interface deduplication | 45,553 |

**Roughly 40,000 distinct protein-protein interfaces exist in the PDB**, whichever
way you count. That is the number to hold in mind: it is what all the experimental
multimer data in the world reduces to, and it is about four orders of magnitude
smaller than our monomer corpora.

### ModelArchive: the established home of predicted complexes, and it is small

ModelArchive is where the AlphaFold-Multimer interactome screens have historically
been deposited. Crawled via its API today: **625,966 models across 3,258
depositions**, of which only **~11,879 are confirmed multi-chain (1.9%)**. A
further 1,490 older single-model deposits carry no chain metadata; extrapolating
those gives ~12,600 complexes, but that number is an estimate.

The complexes are concentrated in about 25 deposits:

| Deposit | Models | Chains | Scope | Paper |
|---|---|---|---|---|
| `ma-bak-evip` | 3,203 | 2 | human bacterial pathogens (RoseTTAFold2-Lite) | [Nat Microbiol 2024](https://doi.org/10.1038/s41564-024-01791-x) |
| `ma-drew-dc2` | 2,573 | 2 | human, **AlphaFold3** (AF3 terms, not CC) | no DOI recorded |
| `ma-bak-cepc` | 1,106 | 2 | *S. cerevisiae* | [Science 2021](https://doi.org/10.1126/science.abm4805) |
| `ma-t3vr3` | 957 | 2 | human cancer interactome | [Protein Sci 2022](https://doi.org/10.1002/pro.4479) |
| `ma-low-csi` | 929 | 2 | human PPI dimers from XL-MS | [PNAS 2023](https://doi.org/10.1073/pnas.2219418120) |
| `ma-dm-prc` | 742 | 2 | human Polycomb PRC1/2 | [Sci Adv 2024](https://doi.org/10.1126/sciadv.adl4529) |
| `ma-dm-hisrep` | 268 | 2–7 | replisome H3-H4 chaperone | [Cell 2024](https://doi.org/10.1016/j.cell.2024.07.006) |

Two things to note. Every large deposit is **pairwise**; even Humphreys et al.
2021, whose paper describes assemblies of up to five subunits, deposited only
2-chain models. And the license is **CC BY-SA 4.0**, with the ShareAlike clause,
which is stricter than AFDB's CC BY and would propagate to anything we
redistribute. Bulk access is per-deposit ZIP only; there is no archive-wide dump.

RCSB mirrors 69,326 ModelArchive models, and my own RCSB census independently
agrees on the complex count: 2,063 multi-chain computed structure models, all
heteromeric, corresponding to `ma-bak-cepc` (1,106) plus `ma-t3vr3` (957).

### BFMD: the aggregated predicted-multimer database

The [Foldseek-Multimer paper](https://doi.org/10.1038/s41592-025-02593-7) (Nature
Methods, 2025) assembled the **Big Fantastic Multimer Database**, which "organizes
297,570 multimer predictions from community efforts into a single database". It is
served as a searchable Foldseek database (`bfmd`, version 20240623, the only
complex-capable database on the Foldseek server besides PDB100). This is the
closest thing to a curated, deduplicated aggregate of the predicted-complex
literature, and at ~300 k it is about 1% the size of the new AFDB dimer release.

## Foldseek-Multimer smoke test

Run on 2026-08-02 against Foldseek `718d4217` to check that the AFDB complex files
can be clustered at all before anyone plans a run over the ~2 M confident models.
Two scripts, both capped and both cheap:

- **`fetch_dimer_sample.py`** pulls 30 heterodimers and 15 homodimers, both gated
  on the authors' criterion validated above, by range-reading the FTP metadata
  CSVs at ten spread-out offsets and downloading each model from
  `https://alphafold.ebi.ac.uk/files/AF-<id>-model_v1.cif`. **45 structures,
  21.1 MB of mmCIF, 12.5 s, zero failures** (the transfer is about a quarter of
  that, since the script asks for gzip and a sample of ten models compressed
  4.1x). `MAX_STRUCTURES = 250` raises rather than truncates, so the script
  cannot be turned into a bulk downloader by accident.
- **`smoke_test_foldseek_multimer.py`** runs `easy-multimercluster` and
  `easy-multimersearch` over that sample.

Foldseek read all 90 chains from the 45 ModelCIF files with no parse warnings, so
AFDB complex files need no conversion. Both commands finished in about 3 seconds.

| | |
|---|---|
| structures / chains ingested | 45 / 90 |
| clusters (`--multimer-tm-threshold 0.5`) | 40 |
| multi-member clusters | 4 (sizes 3, 2, 2, 2) |
| complex pairs in the search report | 650 |
| pairs with multimer TM ≥ 0.5 | 7 |
| pairs where both chains aligned | 199 |

A 1.12x reduction is what a random draw of 45 complexes from a 29 M set should
give, so the number that matters is what the four multi-member clusters contain:

- **bma-lsm-6/bma-lsm-5.1, bma-snr-7/bma-snr-5, bma-lsm-6/bma-snr-6** (*Brugia
  malayi*, tax 6279), the 3-member cluster, multimer TM 0.799 to 0.819. All six
  proteins are Sm/Lsm family subunits, which assemble into a heptameric ring, so
  pairwise decomposition of one ring yields many copies of the same Sm-Sm
  interface.
- **LOC100805631/LOC100789233 with LOC100805631/LOC100796470** (soybean, tax 3847),
  multimer TM 0.876. Shared first partner, paralogous second partners.
- **GDU4/LOG2 with GDU1/LOG2** (*Arabidopsis*, tax 3702), multimer TM 0.670. GDU1
  and GDU4 are paralogous glutamine dumpers bound to the same ubiquitin ligase.
- **SNX4 homodimer (*Merluccius polli*, Benguela hake) with the
  LOC100777531/LOC100803157 heterodimer (soybean)**, multimer TM 0.619, both
  chains aligned. A sorting-nexin BAR domain dimerises into a crescent, and a
  soybean heterodimeric pair reproduces the same interface.

Three lessons from four clusters. **The Sm-ring case defeats sequence-level
filtering**: those three complexes carry six different gene names, so a filter on
gene pairs would have kept all three as distinct training examples, exactly the
lesson exp41 hit on monomers (65 of 99 FoldBench structures overlapped training
structurally while sitting below 30% sequence identity). **The SNX4 case means the
homodimer and heterodimer sets must be deduplicated together**, not separately:
they share interfaces across the two files and across a fish and a legume.

And **the search report is not the clustering**. Seven pairs scored multimer
TM ≥ 0.5 but only six of them formed clusters; the seventh (a soybean
LOC100805119/LOC547737 pair against human STK38L/MOB1A, TM 0.519) aligned a
single chain of each, so it is a fold-level match rather than an interface match
and `easy-multimercluster` correctly declined to merge them. Filter on chain
coverage, not on the TM-score alone.

One trap worth recording, because it fails silently. Passing the structure files as
separate argv entries instead of as a directory makes Foldseek treat all but the
last one or two as something else and **still exit 0**; the first run of this smoke
test reported a clean clustering of 2 chains as if it had covered all 45 structures.
`assert_all_chains_ingested()` now checks the reported database size against
2 × the number of files and raises.

## Conclusion

**Both of our current corpora are strictly monomeric, and this was verified rather
than assumed.** AFDB's 241 M main-release models are one UniProt accession and one
`_struct_asym` record each; the ESM Atlas's 1.10 B + 6.6 M models were sampled at
1,600 structures and every one had a single chain, despite a storage format that
carries `chain_id` / `entity_id` / `sym_id` / `chain_boundaries` fields ready for
complexes.

**The situation changed in 2026, and the change is large.** AFDB now carries a
four-way EMBL-EBI / DeepMind / NVIDIA / Seoul National University release of
21.5 M predicted homodimers and 7.6 M predicted heterodimers (48.8 TB),
produced with AlphaFold-Multimer v2.3.0 weights through OpenFold, scored with
ipTM/ipSAE/pDockQ, carrying UniProt cross-references for both chains, and licensed
CC-BY-4.0 for commercial use. It is roughly 100× larger than every previously
available predicted-complex collection put together. The high-confidence slice is
queryable per accession through `/api/complex/`; the bulk lives on FTP under
`collaborations/nvda/`, which is easy to miss because nothing in the main AFDB
download documentation points at it.

**Filter hard and the heterodimers shrink to something familiar in size.** At the
release's own quality gate (ipSAE ≥ 0.6 and pDockQ2 ≥ 0.23) 0.92% of heterodimers
survive, about **70 k complexes**, which lands in the same range as the ~40 k
distinct interfaces in the entire PDB and BFMD's 297 k aggregated predictions. The
homodimers hold up ten times better: **9.13%, about 1.96 M**, measured by that
same gate after checking it reproduces the authors' declared verdict on all
123,597 heterodimer rows sampled. That is a genuinely new quantity of confident
complex structure, and homodimers are also where the latent demand is: ~2.9 M AFDB
entries model proteins that UniProt annotates as homo-oligomers. Two caveats on
that 1.96 M. Which score you filter on matters more than where you put the cut
(at a fixed threshold of 0.5 the heterodimer yield ranges from 46 k by pDockQ2 to
437 k by ipTM). And it counts models, not distinct proteins: about a quarter of
confident homodimers repeat a gene name already present, because UniProt carries
one accession per strain, and structural clustering will cut deeper than that.

**Experimental multimer data is fixed at about 40,000 distinct interfaces.** RCSB
sequence clustering, PINDER's structural interface clustering and PPIRef's `iDist`
deduplication independently give 37 k, 42 k and 46 k. The PDB adds roughly 6,000
multi-chain entries a year, so this number is not going to move by an order of
magnitude. Any multimer capability has to come mostly from predicted structures,
which is what makes the AFDB release significant.

**The ESM Atlas contributes candidate pairs, not structures.** Roughly half its
provenance records (SPIRE, plus UHGG, IMG/VR and the uMAG sets) encode contig and
ORF position, and 74% of sampled SPIRE entries have an upstream neighbour on the
same contig. That is a very large supply of operon-adjacent candidate interacting
pairs, each with a folded monomer already in hand; but turning it into complex
structures means running a predictor ourselves, not downloading anything.

If we pursue multimers, the practical ordering is: the AFDB homodimer set first
(largest confident yield, permissive license, joins on UniProt), PINDER as the
evaluation anchor (it already ships leakage-controlled splits and AFDB-paired
monomers, so we would not rebuild that), then the AFDB heterodimers under their
quality gate. ModelArchive and BFMD are small enough by comparison to be rounding
errors, and ModelArchive's ShareAlike license makes it the least attractive of
the three. Whatever the order, run Foldseek-Multimer over the homodimers and
heterodimers **together**: the smoke test above found a fish SNX4 homodimer and a
soybean heterodimer sharing an interface at multimer TM 0.62, so deduplicating
the two files separately would leave that redundancy in place.

Two open questions this survey did not settle. The complex release's paper is
still listed as "to be published" (the announcement points at a preprint on
`research.nvidia.com` that I did not retrieve), so the exact pair-selection
procedure inside the 20 species plus WHO-pathogen scope is not documented. And
`ntdx/` and `vr3d/`, two of AFDB's five collaboration datasets, carry no README
and remain unidentified.

## Reproducing the numbers

- `probe_pdb.py`: the PDB census. `python3 probe_pdb.py --out data`. Faceted
  count queries only; writes the CSVs in `data/`.
- `profile_confidence.py`: the confidence-gate tables.
  `python3 profile_confidence.py --windows 200 --window-bytes 200000 --out data`.
  Range requests only, ~80 MB total, about two minutes. Writes
  `nvda_confidence_gates.csv` (every metric at every threshold, with per-window
  spread and estimated totals), `nvda_confidence_quantiles.csv` (the underlying
  distributions), `nvda_confidence_windows.csv` (one row per window, which is the
  evidence for the organism-blocking caveat) and `nvda_confidence_sampling.json`
  (row-count estimates, redundancy counts, and the gate-versus-declared check).
- `fetch_dimer_sample.py`: the capped dimer sampler.
  `python3 fetch_dimer_sample.py --n-hetero 30 --n-homo 15 --out sample`. Writes
  `sample/structures/` (gitignored) plus `sample/manifest.csv` and
  `sample/fetch_summary.json`, which are copied into `data/` as
  `dimer_sample_manifest.csv` and `dimer_sample_fetch_summary.json` so the sample
  is recorded without committing 21 MB of coordinates. `--gate iptm` swaps the
  release's ipSAE/pDockQ2 criterion for a plain ipTM cut.
- `smoke_test_foldseek_multimer.py`: the Foldseek run.
  `python3 smoke_test_foldseek_multimer.py --sample sample --out data`. Needs no
  preinstalled Foldseek; `foldseek_env.py` downloads a static build into
  `~/.cache/marinfold/foldseek` on first use.
- AFDB and ESM Atlas numbers were read interactively; the commands are quoted
  inline above (FTP directory listings, range requests against the metadata CSVs,
  and `lance` / `pyarrow` reads against `s3://esm-protein-atlas` with
  `storage_options={"aws_skip_signature": "true", "region": "us-west-2"}`).
