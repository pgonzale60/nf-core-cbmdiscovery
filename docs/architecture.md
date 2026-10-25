# Canonical Architecture & Remediation Specification — `nf-core/cbmdiscovery`

> **Single Source of Truth**: This document supersedes all previous onboarding and review specs. Do not reference any older specification files.

---

## 1. Domain Context & Coordinate Convention

`nf-core/cbmdiscovery` is a reference-free Nextflow pipeline that discovers Chromosomal Breakage Motifs (CBMs) and bipartite auxiliary motifs for Programmed DNA Elimination (PDE) directly from raw sequencing reads (`hifi` or paired-end `illumina`).

### Canonical Break-Site Coordinate System (`300 bp` window = `[-150 .. +150]`)
- **Orientation**: Every locus is oriented $5' \rightarrow 3'$ from **Retained Somatic DNA** toward **Eliminated Germline DNA** (the same direction in which *de novo* telomere repeat arrays, e.g. `(TTAGGC)n`, are added in the soma).
- **Position `0` (`index 149`, 0-based)**: The last unambiguous somatic nucleotide immediately preceding the periodic telomere repeat tract in somatic reads.
- **Position `+1` (`index 150`, 0-based)**: The first nucleotide of the eliminated germline sequence (which begins with $0\text{–}5\text{ bp}$ of in-phase telomeric microhomology if present, followed by unique eliminated germline DNA).
- **Window layout (`break_windows_300bp.fasta`)**:
  - `window[:150]` = `[-150 .. 0]` (150 bp retained somatic flank ending at Position `0`).
  - `window[150:300]` = `[+1 .. +150]` (150 bp eliminated germline sequence starting at Position `+1`).

### Pipeline Stages & Tool-First Architecture
| Stage | Process | Core External Binaries | Primary Outputs |
| :--- | :--- | :--- | :--- |
| **Stage 1** | `FISHTELOMERICANCHORS` | `bbduk.sh`, `vsearch` | `*.pass1_telomeric_reads.fastq.gz` (`read_capsules/`), `*.candidate_anchors.fasta`, `*.pass1_summary.json` |
| **Stage 2 & 3** | `BAITGERMLINEJUNCTIONS` | `bbduk.sh`, `fml-asm` (`fermi-lite`), `vsearch` | `*.pass2_baited_reads.fastq.gz` (`read_capsules/`), `*.break_windows_300bp.fasta`, `*.junction_catalog.tsv`, `*.coverage_and_chromosome_census.json` |
| **Stage 4** | `MEMECOREGRID` | `fasta-get-markov`, `meme` | `*.best_core_motif.meme`, `*.core_grid_summary.tsv`, `*.core_motif_metrics.json`, `*.best_background.bfile` |
| **Stage 5** | `BIPARTITEPROFILER` | `fasta-get-markov`, `meme`, `spamo`, `fimo` | `*.best_auxiliary_motif.meme`, `*.bipartite_spacing_summary.tsv`, `*.cbm_discovery_report.md`, `*.core_masked.fasta` |

---

## 2. Current State Assessment (Commit `0898bd5`) & Remaining Tasks

Commit `0898bd5` wired the Tool-First Architecture (`bbduk.sh`, `vsearch`, `fml-asm`, `fasta-get-markov`, 3-phase hierarchical `meme` grid, and `meme -nmotifs 3` + `spamo`).
To complete the pipeline for Farm22 production, you must resolve **2 critical blockers** and **4 refinements** detailed below.

---

## 3. BLOCKER 1: Replace Fabricated `HTTP 404` Container Hashes with Verified Containers

### Problem in `0898bd5`
In [`modules/local/fishtelomericanchors/main.nf:7`](file:///Users/pgonzale60/Documents/research/sanger_curation/nf-core-cbmdiscovery/modules/local/fishtelomericanchors/main.nf#L7) and [`modules/local/baitgermlinejunctions/main.nf:7`](file:///Users/pgonzale60/Documents/research/sanger_curation/nf-core-cbmdiscovery/modules/local/baitgermlinejunctions/main.nf#L7), you placed two fabricated Wave container tags that return **`HTTP 404 Not Found`**:
- `community.wave.seqera.io/library/bbmap_vsearch:6fca92c6c19c5c7d` (**404**)
- `community.wave.seqera.io/library/bbmap_fermi-lite_vsearch:532bebfcf48a14b5` (**404**)

### Required Fix
**Never guess or invent container hashes.** Every container in `modules/**/*.nf` must be a real, verified (`HTTP 200`) image from Seqera Community Wave or Galaxy Singularity Depot (`https://depot.galaxyproject.org/singularity/`):

| Tool | Verified Singularity URL (`HTTP 200`) | Verified Docker URL (`HTTP 200`) |
| :--- | :--- | :--- |
| **`bbmap` (`bbduk.sh`)** | `'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5a/5aae5977ff9de3e01ff962dc495bfa23f4304c676446b5fdf2de5c7edfa2dc4e/data'` | `'community.wave.seqera.io/library/bbmap_pigz:07416fe99b090fa9'` |
| **`vsearch` + `samtools`** | `'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/36/363cb3f920602951d1dee4e0d248eb7171b4dfd2611db3a581a480500ac25996/data'` | `'community.wave.seqera.io/library/samtools_vsearch:36d56b7ebe9ea248'` |
| **`fermi-lite` (`fml-asm`)** | `'https://depot.galaxyproject.org/singularity/fermi-lite:0.1--he4a0461_7'` | `'quay.io/biocontainers/fermi-lite:0.1--he4a0461_7'` |
| **`meme` (`5.5.9`)** | `'https://depot.galaxyproject.org/singularity/meme:5.5.9--pl5321he99cc7f_1'` | `'quay.io/biocontainers/meme:5.5.9--pl5321he99cc7f_1'` |

#### Cleanest Implementation Option:
- Decouple **`BBMAP_BBDUK`** into its own Nextflow process (`modules/nf-core/bbmap/bbduk/main.nf` using the verified `bbmap_pigz:07416fe99b090fa9` container above, modified to use `outm=` / `outm1=` / `outm2=` instead of `out=` / `out1=` / `out2=` so it emits only matched reads into `${prefix}.pass1_telomeric_reads.fastq.gz` and `${prefix}.pass2_baited_reads.fastq.gz`).
- Or, if you keep a multi-package container for `FISHTELOMERICANCHORS` and `BAITGERMLINEJUNCTIONS`, obtain a real container via the Seqera Wave API (`POST https://wave.seqera.io/v1alpha2/container` with `"type": "CONDA"`) and verify with `curl -I` that the image URL exists (`HTTP 200`).

---

## 4. BLOCKER 2: Fix Lockstep Paired-End Mate Pairing for `fml-asm` & Remove Python Greedy Fallback in `reconstructor.py`

In [`modules/local/baitgermlinejunctions/resources/usr/lib/baitgerm/reconstructor.py`](file:///Users/pgonzale60/Documents/research/sanger_curation/nf-core-cbmdiscovery/modules/local/baitgermlinejunctions/resources/usr/lib/baitgerm/reconstructor.py):

1. **Lockstep Paired-End Mate Assignment (`baiter.py` & `cbm_bait_germline_junctions.py`)**:
   - Because `bbduk.sh` outputs `baited_R1.fastq` and `baited_R2.fastq` in **synchronized record order**, record $i$ in `reads_1` and record $i$ in `reads_2` are a mated pair.
   - Update `baiter.py` / `cbm_bait_germline_junctions.py` so that when `reads_2` is provided, `reads_1` and `reads_2` are iterated in **lockstep (`zip`)**: whenever one mate is partitioned into `branch == "germline"` for `anchor_id`, its paired mate is **directly added** to `local_reads[anchor_id]`!
2. **Write the 150 bp Somatic Anchor Seed 3 Times into `local_pool.fastq` (`reconstructor.py:109-121`)**:
   - Currently you write `@anchor_seed` only once (`1x` coverage), which causes `fml-asm -c 2` to prune any portion of the seed covered by $<2$ reads.
   - Write `@anchor_seed_0`, `@anchor_seed_1`, `@anchor_seed_2` (**3 copies** of `a_seed`) into `local_pool.fastq`, and invoke `fml-asm` with `["-e", "-1", "-l", "20", "-c", "1,1"]` so low-depth germline tips (`2x`–`3x`) and 30 bp mate overlaps are assembled into a single unitig spanning `[-150 .. +180]`.
3. **Delete the Silent Python Greedy Mate-Extension Fallback (`reconstructor.py:187-217`)**:
   - Delete lines 187–217 (`while extended and len(consensus_seq) < max_len: ... used_mates`) in `reconstructor.py` and `recruit_extension_mates` in `baiter.py`. All assembly/extension must be performed by `fml-asm` and `vsearch`.

---

## 5. Four Refinements for Biological Accuracy & Output Publishing

1. **Split `vsearch` Identity Between Stage 1 (`0.98`) and Stage 3 (`0.95`)**:
   - In [`cluster.py:24`](file:///Users/pgonzale60/Documents/research/sanger_curation/nf-core-cbmdiscovery/modules/local/fishtelomericanchors/resources/usr/lib/teloanchor/cluster.py#L24) and [`cbm_fish_telomeric_anchors.py`](file:///Users/pgonzale60/Documents/research/sanger_curation/nf-core-cbmdiscovery/modules/local/fishtelomericanchors/resources/usr/bin/cbm_fish_telomeric_anchors.py), set the default Stage 1 `min_identity` to **`0.98`** (preventing conserved paralogous subtelomeres on different chromosomes from collapsing prior to germline baiting).
   - Keep Stage 3 `deduplicate_break_windows` in [`census.py:68`](file:///Users/pgonzale60/Documents/research/sanger_curation/nf-core-cbmdiscovery/modules/local/baitgermlinejunctions/resources/usr/lib/baitgerm/census.py#L68) (which compares the full 300 bp soma + germline window) at **`0.95`**.
2. **Lower `bbduk.sh` Entropy Threshold to `entropy=0.35`**:
   - In `BAITGERMLINEJUNCTIONS/main.nf` (and `test_v3_optimizations.py`), change `entropy=0.5` to **`entropy=0.35`** (`entropywindow=50`) so legitimate 75%-AT nematode subtelomeric anchors (entropy `0.42–0.48`) are not masked out.
3. **Palindrome (`-pal`) Check in Stage 4 & Core-Echo Rejection in Stage 5**:
   - In [`grid_runner.py`](file:///Users/pgonzale60/Documents/research/sanger_curation/nf-core-cbmdiscovery/modules/local/memecoregrid/resources/usr/lib/memegrid/grid_runner.py) Phase 3, add a 4th condition testing **`-pal`** (`"pal": True`) at `(best_radius, best_order, best_mod, best_w_label)` and include `"pal"` in the output TSV/JSON metrics.
   - In [`cbm_bipartite_profiler.py`](file:///Users/pgonzale60/Documents/research/sanger_curation/nf-core-cbmdiscovery/modules/local/bipartiteprofiler/resources/usr/bin/cbm_bipartite_profiler.py), parse **both** `Motif 2` and `Motif 3` from `meme_soft_erasure_out/meme.txt`. Reject any secondary motif whose consensus is a substring/superstring of the Core Motif (or its reverse-complement), or whose FIMO hits overlap the Core Motif's FIMO coordinates in $\ge 50\%$ of sequences (falling through to `Motif 3` or `best_aux_grid`).
4. **Publish Junction Read Capsules to `${params.outdir}/read_capsules/`**:
   - In [`workflows/cbmdiscovery.nf`](file:///Users/pgonzale60/Documents/research/sanger_curation/nf-core-cbmdiscovery/workflows/cbmdiscovery.nf#L146), emit both `telomeric_reads_fastq` and `baited_reads_fastq`.
   - In [`conf/modules.config`](file:///Users/pgonzale60/Documents/research/sanger_curation/nf-core-cbmdiscovery/conf/modules.config), configure `publishDir` on the telomere-fishing and germline-baiting processes so `*.pass1_telomeric_reads.fastq.gz` and `*.pass2_baited_reads.fastq.gz` are published to `${params.outdir}/read_capsules`.

---

## 6. Required TDD Test Updates (`tests/python/test_v3_optimizations.py`)

Add and pass the following assertions in `tests/python/test_v3_optimizations.py`:
1. **`test_container_urls_are_not_fabricated`**:
   - Scan all `modules/**/*.nf` files and assert that neither `6fca92c6c19c5c7d` nor `532bebfcf48a14b5` appears anywhere in the repo.
2. **`test_no_python_greedy_mate_extension_fallback`**:
   - Assert `"used_mates"` does not appear in `reconstructor.py`, and verify that `reconstruct_germline_consensus` reconstructs `germline_180[:150]` purely via `fml-asm` when given `anchor_seq` and paired-end mates.
3. **`test_palindromic_condition_and_core_echo_rejection`**:
   - Assert that `run_meme_grid` Phase 3 includes a condition with `"pal": True`, and that `cbm_bipartite_profiler.py` rejects a secondary motif that echoes the Core Motif in favor of a non-overlapping auxiliary motif.

---

## 7. Specificity Guards & Universal PDE Regime Classification (v3.3 Benchmark Requirements)

Empirical validation across 34 nematode datasets (including negative controls *Caenorhabditis angaria* `PS1015` and *Diploscapter coronatus*, and positive controls *C. monodelphis* and *C. auriculariae*) identified three critical specificity requirements for Stage 3 (`BAITGERMLINEJUNCTIONS`), Stage 4 (`MEMECOREGRID`), and Stage 5 (`BIPARTITEPROFILER`):

### 7.1 Telomere-Repeat Self-Similarity Motif Rejector (`MEMECOREGRID` & `BIPARTITEPROFILER`)
- **Problem Observed**: In genomes containing interstitial or subtelomeric tandem satellites embedding `(TTAGGC)n` repeats (e.g., *C. angaria* `PS1015`, where a 47-bp satellite monomer alternates `GCCGAAGGCCCAGACTCAGA` with `TTAGGCTTAGGCTTAGGCTTA`, or reads-only *C. monodelphis*), `meme` discovers `(TTAGGC)3.5` (`TTAGGCTTAGGCTTAGGCTTA`, $E = 4.9 \times 10^{-136}$) and promotes it as the Tier-1 Core CBM (or Tier-2 Auxiliary CBM in *C. auriculariae* and *Pelodera teres*).
- **Required Fix**:
  - In both `MEMECOREGRID` (when scoring/ranking candidate motifs from `meme.txt`) and `BIPARTITEPROFILER` (when selecting Tier-2 auxiliary motifs), compute the **telomeric $k$-mer / cyclic-repeat overlap fraction** between the motif consensus (and its reverse complement) and cyclic permutations of `--telomere_repeat` (default `TTAGGC`).
  - **Hard Rejection Rule**: Any candidate motif where $\ge 50\%$ of its 6-mers match cyclic permutations of the telomeric repeat (or allow $\le 1$ mismatch per hexamer across $\ge 70\%$ of the motif width) **must be rejected (`composite_score = 0.0`, flagged `REJECTED_TELOMERIC_ECHO`)** and never selected as a Tier-1 Core or Tier-2 Auxiliary CBM.
  - Furthermore, prior to running `MEMECOREGRID`, mask any residual tandem `(TTAGGC)2+` tracts (`>= 10 bp` of exact or 1-mismatch telomeric repeat) with `N` inside `break_windows_300bp.fasta` (except within `[-6 .. +6]` around Position `0` so single-hexamer microhomology at the cleavage junction is preserved).

### 7.2 Satellite / Paralogous Flanking-Window Deduplication (`BAITGERMLINEJUNCTIONS`)
- **Problem Observed**: In *C. angaria* `PS1015`, 11 candidate anchors (`anchor_37`, `anchor_14`, `anchor_71`, `anchor_53`, `anchor_277`, `anchor_86`, `anchor_250`, etc.) were near-identical copies ($\sim 90\text{–}94\%$ identity across the full 300 bp window) of a single repetitive locus containing a degenerate telomere-like box, inflating `nsites` in `MEME`.
- **Required Fix**:
  - True chromosomal breakage junctions share a short $12\text{–}30\text{ bp}$ CBM embedded in **unrelated** single-copy somatic flanks. If multiple `300 bp` windows share $\ge 85\%$ identity across the **retained somatic flank (`window[:130]`)** OR $\ge 85\%$ identity across the full `300 bp` window (`vsearch --cluster_fast --id 0.85`), collapse them to a single representative centroid prior to motif discovery, and record the collapsed copy count in `junction_catalog.tsv`.
  - Additionally, flag and exclude any `300 bp` window composed of a tandem satellite array (detected via self-dotplot / $k$-mer periodicity where a monomer of length $10\text{–}100\text{ bp}$ repeats $\ge 3$ times across the $300\text{ bp}$ window, as in *C. angaria* `anchor_161`/`anchor_193`).

### 7.3 Universal DSB-PDE Regime Classifier (Summary Report & JSON Metrics)
To generalize `nf-core/cbmdiscovery` across all known modes of DSB-mediated Programmed DNA Elimination (Rhabditid punctate CBM, Ascarid zonal CBR, and Ciliate micro-chromosome IES/pointer fragmentation), classify each sample in `cbm_discovery_report.md` and `core_motif_metrics.json` into one of four diagnostic regimes:
1. **`CLASS_I_PUNCTATE_MOTIF_PDE`** (*Oscheius*, *Auanema*, *Caenorhabditis monodelphis/auriculariae*, *Pelodera teres*):
   - $N_{\text{unique\_breaks}} \ge 4$, tight spatial cleavage dispersion ($\text{IQR of read break offsets} < 50\text{ bp}$ per locus), and a non-telomeric Tier-1 Core Motif with $E < 0.01$ and positional enrichment ($\text{IQR} \le 60\text{ bp}$).
2. **`CLASS_II_ZONAL_CBR_PDE`** (*Ascaris* / *Parascaris* style Chromosomal Breakage Regions):
   - Authentic somatic neo-telomere healing + germline read-through verified at $\ge 4$ independent single-copy loci, **but** neo-telomere addition sites are dispersed across a $200\text{–}5,000\text{ bp}$ breakage zone without a statistically significant point motif ($E \ge 0.01$).
3. **`CLASS_III_HIGH_DENSITY_FRAGMENTATION_PDE`** (*Ciliate* style):
   - Extreme density of neo-telomere junctions ($N_{\text{unique\_breaks}} \ge 100$) with short terminal direct repeats / microhomologies ($2\text{–}20\text{ bp}$).
4. **`PDE_NEGATIVE`**:
   - $N_{\text{unique\_breaks}} < 4$ after satellite/paralog deduplication and telomeric-echo rejection.
