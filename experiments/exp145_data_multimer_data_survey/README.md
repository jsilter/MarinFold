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

**The AlphaFold DB acquired a multimer corpus in the last two months and it is
enormous.** An NVIDIA collaboration release, published only on the AFDB FTP site
(not the web UI, not the REST API), holds roughly **21 M predicted homodimers and
7.6 M predicted heterodimers (about 29 M dimers and 48.8 TB)**. That is more
predicted complexes than everything else in this report combined, by roughly two
orders of magnitude. It is scored (ipTM, ipSAE, pDockQ, pDockQ2, LIS, clash
counts) and cross-referenced to UniProt on both chains, so it can be filtered to
a confident subset without re-running anything.

Everything else is much as expected:

| Source | Multimer content | Non-redundant scale | Access |
|---|---|---|---|
| AFDB main release (v6) | **none**; 241,070,489 single-chain models | n/a | FTP + GCS, public |
| AFDB `collaborations/nvda/` | **~29 M predicted dimers** (~21 M homo, ~7.6 M hetero) | ~1.5–1.9 M homo, ~70 k hetero pass confidence gates | FTP only, 48.8 TB |
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
download for a line count, so I pulled 64 windows of 300 kB at evenly spaced byte
offsets in each and scaled by mean row length (row lengths are tight: 276 ± 30
bytes for homodimers, 323 ± 29 for heterodimers), which puts the counting error
near 1%:

| | shards | model data | metadata | estimated rows |
|---|---|---|---|---|
| homodimers | 9,892 | 35.2 TB | 6.00 GB | **~21 M** |
| heterodimers | 8,209 | 13.6 TB | 2.45 GB | **~7.6 M** |

Every model carries interface confidence scores (`ipTM`, `ipSAE` in both chain
orders, `pDockQ`, `pDockQ2` in both orders, `LIS`, `numberOfInteractions`,
`N_clash_backbone`, `N_clash_heavyAtom`), plus UniProt accessions, taxon IDs and
gene names for both partners.

One methodological warning that matters for reading the next table: **the CSVs are
blocked by organism**, so a byte-offset window lands inside a single species' run
of rows and confidence correlates strongly with which species you hit. Across the
64 windows (67,801 homodimer and 59,158 heterodimer rows) the pass rates were:

| Confidence gate | homodimers | heterodimers |
|---|---|---|
| ipTM ≥ 0.5 | 19.8% (range 4.3–47.0) → ~4.2 M | 6.0% (range 3.1–8.8) → ~0.45 M |
| ipTM ≥ 0.8 | 7.2% (range 0.5–27.9) → ~1.5 M | 0.9% (range 0.1–3.6) → ~0.07 M |
| pDockQ ≥ 0.5 | 18.4% (range 2.6–53.3) → ~3.9 M | 3.1% (range 0.2–9.7) → ~0.24 M |
| max(ipSAE) ≥ 0.75 | 8.8% (range 0.5–25.8) → ~1.9 M | 0.7% (range 0.0–2.1) → ~0.05 M |
| authors' `passes_quality_threshold` | (column absent) | 1.0% (range 0.2–4.4) → **~75 k** |

The per-window ranges are wide, so treat the totals as order-of-magnitude, not
three-significant-figure.

The heterodimer table ships its own verdict column: `passes_quality_threshold`
with `quality_ipsae_threshold = 0.6` and `quality_pdockq2_threshold = 0.23`. Only
about **1%** of heterodimers pass it, so the confident heterodimer set is on the
order of **70,000 complexes**, the same order as the PDB's non-redundant interface
count, not the 7.6 M headline. Homodimers do much better (no author gate, but
~1.5–1.9 M clear ipTM ≥ 0.8 / ipSAE ≥ 0.75). That gap is expected: a homodimer
interface is far easier to predict than an arbitrary heterodimer pairing, and with
an ipTM early-stopping threshold of 0.1 the heterodimer set deliberately retains a
large mass of speculative pairings rather than filtering them out.

#### Coverage

A separate pass (10 windows of 2 MB, giving 72,286 homodimer and 61,762
heterodimer rows) shows the two sets are shaped quite differently:

- **Homodimers are one model per protein.** Every `uniprotAccession` in the sample
  was unique (72,286 distinct accessions in 72,286 rows); they ran the monomer
  set through a dimer prediction rather than screening pairs.
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

**Access is FTP-only.** The AFDB REST API returns `{}` for these model entity IDs
(`AF-0000000065760001`, `AF-0000000203470222`) and they do not appear in the AFDB
web UI. Nothing but the FTP path knows they exist.

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
quaternary structures"*, marked **"To be published"**, so the corpus is out ahead
of its paper.

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

## Conclusion

**Both of our current corpora are strictly monomeric, and this was verified rather
than assumed.** AFDB's 241 M main-release models are one UniProt accession and one
`_struct_asym` record each; the ESM Atlas's 1.10 B + 6.6 M models were sampled at
1,600 structures and every one had a single chain, despite a storage format that
carries `chain_id` / `entity_id` / `sym_id` / `chain_boundaries` fields ready for
complexes.

**The situation changed in June–July 2026, and the change is large.** AFDB's FTP
site now carries an NVIDIA collaboration release of roughly 21 M predicted
homodimers and 7.6 M predicted heterodimers (48.8 TB), produced with
AlphaFold-Multimer v2.3.0 weights through OpenFold, scored with ipTM/ipSAE/pDockQ,
carrying UniProt cross-references for both chains, and licensed CC-BY-4.0 for
commercial use. It is invisible from the AFDB web UI and API, its paper is
unpublished, and it is roughly 100× larger than every previously available
predicted-complex collection put together.

**Filter hard and it shrinks to something familiar in size.** At the authors' own
heterodimer quality gate only ~1% survive (~70 k complexes), which lands in the
same range as the ~40 k distinct interfaces in the entire PDB and BFMD's 297 k
aggregated predictions. The homodimers hold up much better (on the order of
1.5–1.9 M pass ipTM ≥ 0.8 or ipSAE ≥ 0.75), which is a genuinely new quantity of
confident complex structure, and homodimers are also where the latent demand is:
~2.9 M AFDB entries model proteins that UniProt annotates as homo-oligomers.

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
the three.

Two open questions this survey did not settle: the nvda release has no paper yet,
so how the heterodimer candidate pairs were chosen is unknown (it matters, because
it determines whether the confident 1% is a biased slice); and `ntdx/` and `vr3d/`,
two of AFDB's five collaboration datasets, carry no README and remain
unidentified.

## Reproducing the numbers

- `probe_pdb.py`: the PDB census. `python3 probe_pdb.py --out data`. Faceted
  count queries only; writes the CSVs in `data/`.
- AFDB and ESM Atlas numbers were read interactively; the commands are quoted
  inline above (FTP directory listings, range requests against the metadata CSVs,
  and `lance` / `pyarrow` reads against `s3://esm-protein-atlas` with
  `storage_options={"aws_skip_signature": "true", "region": "us-west-2"}`).
