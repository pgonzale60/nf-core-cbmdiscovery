# nf-core/cbmdiscovery: Output

## Introduction

This document describes the output produced by the `nf-core/cbmdiscovery` pipeline. The pipeline performs reference-free discovery, deconvolution, and structural characterization of Programmed DNA Elimination (PDE) Chromosome Breakage Motifs (CBM) from raw whole-genome sequencing reads (PacBio HiFi or Illumina short reads).

By default or when hierarchical publishing is active (`--per_species_outdir`), outputs are organized per specimen within species-level hierarchies. Cohort-level QC reports and execution metadata are published at the top-level directory.

The directory layout and output files described below will be created in your results directory (`--outdir`) after the pipeline has finished.

```tree
<outdir>/
├── <species>/
│   └── <specimen>/
│       ├── 01_anchor_fishing/
│       │   ├── <prefix>.candidate_anchors.fasta
│       │   └── <prefix>.pass1_summary.json
│       ├── 02_germline_baiting/
│       │   ├── <prefix>.break_windows_300bp.fasta
│       │   ├── <prefix>.junction_catalog.tsv
│       │   └── <prefix>.coverage_and_chromosome_census.json
│       ├── 03_core_motif_grid/
│       │   ├── <prefix>.best_core_motif.meme
│       │   ├── <prefix>.core_grid_summary.tsv
│       │   ├── <prefix>.core_motif_metrics.json
│       │   └── <prefix>.best_background.bfile
│       ├── 04_bipartite_profiling/
│       │   ├── <prefix>.best_auxiliary_motif.meme
│       │   ├── <prefix>.bipartite_spacing_summary.tsv
│       │   ├── <prefix>.cbm_discovery_report.md
│       │   └── <prefix>.core_masked.fasta
│       └── read_capsules/
│           ├── <prefix>.pass1_telomeric_reads.fastq.gz
│           └── <prefix>.pass2_baited_reads.fastq.gz
├── multiqc/
│   ├── multiqc_report.html
│   ├── multiqc_data/
│   └── multiqc_plots/
└── pipeline_info/
    ├── execution_report.html
    ├── execution_timeline.html
    ├── execution_trace.txt
    ├── pipeline_dag.html
    ├── nf_core_cbmdiscovery_software_mqc_versions.yml
    └── params.json
```

---

## 300 bp Flank Radius Convention (`[-300 .. +300]` = 600 bp Total Window)

A central architectural convention in `nf-core/cbmdiscovery` is the definition of sequence windows around programmed chromosome breakage junctions:

```text
       Retained Somatic Flank                     Eliminated Germline Flank
   [-300 bp ....................... 0]  |  [+1 ....................... +300 bp]
<======================================|-------------------------------------->
                                   Position 0
                           (Cleavage & Telomere Addition)
<---------------------------- Total Window = 600 bp --------------------------->
```

- **Position 0**: Represents the exact cleavage site where germline DNA is cleaved and _de novo_ telomeric repeats are synthesized on the retained somatic chromosome end.
- **Retained Somatic Flank (`[-300 .. 0]` bp)**: Up to 300 bp of sequence immediately upstream of Position 0 that is preserved in somatic cells.
- **Eliminated Germline Flank (`[+1 .. +300]` bp)**: Up to 300 bp of sequence immediately downstream of Position 0 reconstructed from unclipped germline read-through reads. This region is eliminated during PDE.
- **Total Window Sequence Length (600 bp)**: The file `*.break_windows_300bp.fasta` contains sequences spanning **600 bp total** (`300 bp retained soma + 300 bp eliminated germline`). The label `300bp` refers strictly to the **flank radius** on either side of Position 0, not the total window width.
- **Multi-Radius Slicing in MEME Grid**: Downstream Tier-1 motif discovery (`MEMECOREGRID`) slices this master 600 bp window at symmetric radii (e.g., 30 bp, 60 bp, 100 bp, 200 bp, 300 bp) centered at Position 0 to systematically evaluate motif enrichment as a function of distance from the breakage junction.

---

## Pipeline Overview

The pipeline processes data through sequential discovery stages:

- [Stage 1: Somatic Anchor Fishing (`01_anchor_fishing/`)](#stage-1-somatic-anchor-fishing) - Identification of candidate somatic breakage junctions with _de novo_ telomeres.
- [Stage 2: Germline Junction Baiting & Census (`02_germline_baiting/`)](#stage-2-germline-junction-baiting--census) - K-mer read baiting, germline reconstruction, microhomology annotation, and chromosome census.
- [Stage 3: Core Motif Grid Search (`03_core_motif_grid/`)](#stage-3-core-motif-grid-search) - Multi-radius slicing and MEME grid search to identify the optimal core CBM motif.
- [Stage 4: Bipartite Profiling & Reporting (`04_bipartite_profiling/`)](#stage-4-bipartite-profiling--reporting) - Core motif masking, auxiliary motif discovery, bipartite spacer profiling, and discovery report generation.
- [Read Capsules (`read_capsules/`)](#read-capsules) - Minimal FASTQ capsules of telomeric and junction-baited reads.
- [MultiQC (`multiqc/`)](#multiqc) - Aggregate quality metrics and run summaries across all samples.
- [Pipeline Information (`pipeline_info/`)](#pipeline-information) - Nextflow execution trace, resource utilization metrics, and software versions.

---

### Stage 1: Somatic Anchor Fishing

<details markdown="1">
<summary>Output files</summary>

- `01_anchor_fishing/`
  - `*.candidate_anchors.fasta`: FASTA file containing candidate somatic anchor sequences oriented at the telomere junction (Position 0).
  - `*.pass1_summary.json`: JSON file summarizing Pass 1 telomeric read fishing metrics.

</details>

`FISHTELOMERICANCHORS` fishes sequencing reads containing terminal or interstitial canonical telomeric repeats (e.g., `(TTAGGC)_n` in nematodes or `(TTAGGG)_n` in vertebrates). For each candidate read, the module trims telomeric repeats to isolate the non-telomeric genomic flank and constructs candidate somatic anchors oriented toward the breakage junction (Position 0).

`*.pass1_summary.json` provides fishing statistics including the total number of telomeric reads detected, candidate anchors generated, and sequencing platform parameters.

---

### Stage 2: Germline Junction Baiting & Census

<details markdown="1">
<summary>Output files</summary>

- `02_germline_baiting/`
  - `*.break_windows_300bp.fasta`: Reconstructed break windows spanning 600 bp total (`[-300 .. +300]` bp flank radius centered at Position 0) for confirmed `PDE_BREAK_SITE` loci.
  - `*.junction_catalog.tsv`: Comprehensive TSV catalog classifying all analyzed candidate junctions.
  - `*.coverage_and_chromosome_census.json`: Machine-readable JSON summary of somatic/germline coverage, break site counts, and estimated germline chromosome numbers.

</details>

`BAITGERMLINEJUNCTIONS` extracts terminal 60 bp k-mers from candidate somatic anchors to bait germline read-through reads directly from raw sequencing data. It reconstructs the unbroken germline consensus across the junction, deconvolves somatic vs. germline bifurcation coverage, annotates microhomology, and computes chromosome census statistics.

#### `*.junction_catalog.tsv`

The junction catalog classifies every evaluated locus into one of seven mutually exclusive classifications:

| Column                | Description                                                                                                                                                                                |
| --------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `anchor_id`           | Identifier of the candidate anchor locus                                                                                                                                                   |
| `classification`      | Locus status (`PDE_BREAK_SITE`, `GERMLINE_NATIVE_TELOMERE`, `SHARED_NATIVE_TELOMERE`, `TANDEM_SATELLITE_ARRAY`, `UNBROKEN_GERMLINE_LOCUS`, `LOW_COHERENCE_JUNCTION`, `CANDIDATE_TELOMERE`) |
| `c_telo`              | Somatic coverage (count of reads supporting the telomeric somatic branch)                                                                                                                  |
| `c_germ`              | Germline coverage (count of unclipped germline read-through reads)                                                                                                                         |
| `coherence`           | Ratio quantifying branching coherence (`c_telo / (c_telo + c_germ)` or branch purity)                                                                                                      |
| `microhomology_len`   | Length in base pairs of microhomology at the somatic-germline transition                                                                                                                   |
| `microhomology_motif` | Nucleotide sequence of the shared microhomology motif                                                                                                                                      |
| `soma_flank_len`      | Length of reconstructed somatic flank sequence                                                                                                                                             |
| `germ_flank_len`      | Length of reconstructed germline flank sequence                                                                                                                                            |
| `collapsed_copies`    | Number of redundant or identical anchor copies collapsed into this locus                                                                                                                   |

#### `*.coverage_and_chromosome_census.json`

Conforms strictly to the `cbm_census.schema.json` contract:

- `median_c_soma`: Median somatic coverage across all candidate break sites.
- `median_c_germ`: Median germline coverage across all candidate break sites.
- `germline_fraction`: Estimated germline-to-total cellular fraction in the sequenced specimen.
- `n_breaks`: Total number of validated `PDE_BREAK_SITE` junctions.
- `n_germ_telo`: Total number of detected germline native telomeres.
- `estimated_germline_chromosomes`: Estimated haploid germline chromosome count (`n_germ_telo / 2`).
- `n_shared_telo`: Count of telomeres shared between somatic and germline states.
- `n_satellite_excluded`: Number of loci excluded due to tandem satellite repeat arrays.

---

### Stage 3: Core Motif Grid Search

<details markdown="1">
<summary>Output files</summary>

- `03_core_motif_grid/`
  - `*.best_core_motif.meme`: MEME format representation of the optimal discovered core CBM motif.
  - `*.core_grid_summary.tsv`: Comprehensive TSV recording evaluation metrics for all tested parameter grid conditions.
  - `*.core_motif_metrics.json`: Detailed quantitative metrics for the winning core motif condition.
  - `*.best_background.bfile`: 0-order Markov background model file corresponding to the winning condition.

</details>

`MEMECOREGRID` slices the 600 bp break windows at multiple radii centered at Position 0 and performs an exhaustive parameter grid search across MEME motif discovery models (`zoops`, `anr`, `oops`), Markov background orders, and motif widths. Candidate motifs are scored and ranked via a composite objective function that rewards high information content, sharp positional clustering relative to Position 0, and high site occupancy while penalizing telomeric echo artifacts.

#### `*.core_grid_summary.tsv`

Documents every evaluated condition in the grid:

| Column                    | Description                                                       |
| ------------------------- | ----------------------------------------------------------------- |
| `phase`                   | Grid execution phase (`phase1_wide` or `phase2_refine`)           |
| `condition_id`            | Unique condition identifier encoding grid parameters              |
| `radius`                  | Window slicing flank radius in base pairs centered at Position 0  |
| `markov_order`            | Background Markov model order (e.g. `0`)                          |
| `mod`                     | MEME site distribution model (`zoops`, `anr`, or `oops`)          |
| `revcomp`                 | Whether reverse-complement strand searching was enabled           |
| `pal`                     | Palindromic motif search flag                                     |
| `w_label`                 | Motif width search label                                          |
| `width`                   | Discovered motif width in base pairs                              |
| `nsites`                  | Number of contributing motif occurrence sites                     |
| `evalue`                  | Statistical E-value reported by MEME                              |
| `consensus`               | IUPAC consensus sequence of the discovered motif                  |
| `total_ic`                | Total information content across all positions (bits)             |
| `mean_ic`                 | Mean information content per column (bits)                        |
| `consensus_entropy`       | Shannon entropy of the consensus sequence                         |
| `occupancy`               | Fraction of break window sequences containing the motif           |
| `median_offset_to_pos0`   | Median offset in base pairs from motif center to Position 0       |
| `positional_iqr`          | Interquartile range of motif offsets relative to Position 0       |
| `telomeric_echo_fraction` | Fraction of motif occurrences matching telomeric repeat hexamers  |
| `rejection_flag`          | Rejection classification code (if filtered out, otherwise `NONE`) |
| `composite_score`         | Final composite ranking score                                     |

---

### Stage 4: Bipartite Profiling & Reporting

<details markdown="1">
<summary>Output files</summary>

- `04_bipartite_profiling/`
  - `*.best_auxiliary_motif.meme`: MEME format file containing the secondary/auxiliary motif partner (if significant).
  - `*.bipartite_spacing_summary.tsv`: TSV recording sequence-level inter-motif edge-to-edge spacing and cleavage offsets.
  - `*.cbm_discovery_report.md`: Markdown summary report documenting discovered CBM architecture and chromosome census.
  - `*.core_masked.fasta`: 600 bp break windows FASTA with core CBM occurrences masked with `N` (also referenced as `*.core_masked_windows.fasta`).

</details>

`BIPARTITEPROFILER` masks all occurrences of the primary core motif within the 600 bp break windows with `N` and conducts a Tier-2 auxiliary motif search to identify paired bipartite partners. It then profiles edge-to-edge physical spacing (`spacer_bp`), relative orientation (`same` vs. `opposite` strand), and cleavage distance offsets relative to Position 0.

#### `*.bipartite_spacing_summary.tsv`

| Column                   | Description                                                                    |
| ------------------------ | ------------------------------------------------------------------------------ |
| `sequence_name`          | Name of the break window FASTA record                                          |
| `core_start`             | 1-based start coordinate of the core motif in the 600 bp window                |
| `core_stop`              | 1-based stop coordinate of the core motif in the 600 bp window                 |
| `core_strand`            | Strand orientation of the core motif (`+` or `-`)                              |
| `core_offset_to_pos0`    | Signed distance in base pairs from core motif to Position 0                    |
| `has_partner`            | Boolean indicating whether an auxiliary motif partner is present               |
| `partner_type`           | Classification of the auxiliary partner (or `NONE`)                            |
| `partner_domain`         | Sequence domain where partner resides (`somatic` or `germline`)                |
| `partner_start`          | 1-based start coordinate of the auxiliary partner                              |
| `partner_stop`           | 1-based stop coordinate of the auxiliary partner                               |
| `partner_strand`         | Strand orientation of the auxiliary partner (`+` or `-`)                       |
| `partner_offset_to_pos0` | Signed distance in base pairs from partner motif to Position 0                 |
| `spacer_bp`              | Edge-to-edge physical distance in base pairs between core and auxiliary motifs |
| `spacer_type`            | Spacer category (`FIXED`, `VARIABLE`, or `NONE`)                               |
| `relative_strand`        | Relative strand configuration (`same` or `opposite`)                           |

#### `*.cbm_discovery_report.md`

A standalone Markdown report summarizing:

- Chromosome and coverage census (estimated germline chromosomes, break site counts, and median depths).
- Discovered core motif consensus, information content, and distance to cleavage junction.
- Auxiliary motif discovery results and bipartite spacing distribution.
- Microhomology profile across break junctions.

---

### Read Capsules

<details markdown="1">
<summary>Output files</summary>

- `read_capsules/`
  - `*.pass1_telomeric_reads.fastq.gz`: Gzipped FASTQ capsule containing raw reads with terminal or interstitial telomeric repeats fished in Stage 1.
  - `*.pass2_baited_reads.fastq.gz`: Gzipped FASTQ capsule containing reads baited across candidate somatic-germline junctions in Stage 2.

</details>

Read capsules represent targeted subsets of sequencing reads enriched for biologically critical junctions. They enable reproducible local re-analysis, IGV inspection, and validation without requiring retention of terabyte-scale raw FASTQ or BAM inputs.

---

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: Interactive standalone HTML report summarizing whole-pipeline QC metrics.
  - `multiqc_data/`: Directory containing parsed data tables and statistics from all pipeline modules.
  - `multiqc_plots/`: Directory containing publication-ready static plots in vector and raster formats.

</details>

[MultiQC](https://multiqc.info) aggregates quality metrics across all executed tools in the pipeline. It visualizes raw read quality (FastQC), software versions, and custom workflow execution statistics into an interactive report.

---

### Pipeline Information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - `execution_report.html`: Nextflow execution report with resource utilization (CPU, memory, walltime) per process.
  - `execution_timeline.html`: Interactive timeline of process start, execution, and finish times.
  - `execution_trace.txt`: Tab-separated execution trace detailing task status, exit codes, and resource metrics.
  - `pipeline_dag.html` / `pipeline_dag.svg`: Directed Acyclic Graph (DAG) visualizing workflow tasks and dependencies.
  - `nf_core_cbmdiscovery_software_mqc_versions.yml`: Collated YAML file reporting software versions for all pipeline steps.
  - `params.json`: Exact parameter map utilized for the pipeline execution.

</details>

Nextflow execution reports provide full traceability and auditability for workflow executions on local workstations, institutional HPC clusters, and cloud environments.
