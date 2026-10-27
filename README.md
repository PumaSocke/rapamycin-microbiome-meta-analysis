# Rapamycin and the mouse gut microbiome: a meta-analysis

Analysis code for:

> Ortega-Matienzo A, Meiners F, Newels P, Zude I, Fuellen G, Barrantes I.
> **Meta-analysis identifies consistent rapamycin-associated shifts in the mouse gut microbiome.**
> *Manuscript under review.*

We re-analysed raw 16S rRNA sequencing data from three public mouse studies (n = 54) covering dietary and
intraperitoneal (IP) rapamycin treatment. Each study was processed with the same pipeline and then pooled using
two meta-analytic approaches: (1) per-study zero-inflated beta regression (GAMLSS, `metamicrobiomeR`) combined by
random-effects meta-analysis, and (2) batch-corrected linear modelling (`MMUPHin`).

## Studies

| Study | Accession | 16S region | Groups | Rapamycin |
|-------|-----------|------------|--------|-----------|
| Bitto et al. (2016) | [PRJEB13256](https://www.ebi.ac.uk/ena/browser/view/PRJEB13256) | V3–V4 | Control / Rapamycin, IP and diet | 8 mg/kg IP (3 mo); 126 ppm diet (90 d) |
| Yang et al. (2022) | [PRJNA717340](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA717340) | V4 | Control / Rapamycin IP | 2.5 mg/kg IP (14 d) |
| Han et al. (2021) | [PRJNA436876](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA436876) | V3–V4 | Control / Rapamycin IP | 0.5 mg/kg IP (84 d) |

Study numbers (2, 5, 7) are the IDs from the systematic search and are used throughout the code.

## Repository layout

```
pipeline/   Cluster steps: raw-read download from ENA and the QIIME2 + DADA2 pipeline
scripts/    R analysis, run in numbered order
data/       QIIME2 exports per study (feature table, taxonomy, tree, metadata, PICRUSt2 pathways)
            and processed/ phyloseq objects that are the input to the meta-analysis
results/    Created by the scripts (not tracked)
```

## Workflow

**1. Upstream processing (HPC cluster, optional).** Raw reads are not stored here. To regenerate the QIIME2
exports in `data/`, download the reads with `pipeline/download_*.sh` and run `pipeline/qiime2_pipeline.sh`
(QIIME2 2020.08 in Singularity, DADA2 denoising, SILVA 138 classifier). PICRUSt2 was run on the same outputs.

**2. Per-study analysis.** Each script builds a phyloseq object from `data/<study>/`, runs the single-study analyses
(diversity, composition, core microbiome, differential abundance, predicted pathways) and writes the cleaned
phyloseq object to `data/processed/`.

```
scripts/01_study2_Bitto2016_PRJEB13256.R
scripts/02_study5_Yang2022_PRJNA717340.R
scripts/03_study7_Han2021_PRJNA436876.R
```

**3. Meta-analysis.** Both scripts start from `data/processed/ps_*_clean.rds`, so they can be run without step 2.

```
scripts/04_meta_MMUPHin.R                  # batch correction, PERMANOVA, MMUPHin lm_meta, dose-response
scripts/05_meta_metamicrobiomeR_GAMLSS.R   # GAMLSS-BEZI per study, random-effects pooling, IP-only sensitivity,
                                           # dose meta-regression (metafor)
```

`scripts/MicrobiomeFunctions.R` holds shared helper functions and is sourced by the per-study scripts.

## Running

Open R in the repository root; paths are resolved with [`here`](https://here.r-lib.org/), so no paths need editing.
Scripts were written to be run interactively section by section in RStudio.

Main R packages: phyloseq, microbiome, vegan, DESeq2, MMUPHin, metamicrobiomeR, gamlss, metafor, meta, tidyverse,
ggplot2, ggpubr, patchwork, pheatmap, ggVennDiagram. Analyses were run with R 4.4.2 and metafor 4.8-0.

## Data availability

No new sequencing data were generated. All raw data are public under the accessions above. The files in `data/`
are derived from those datasets.

## Citation

If you use this code, please cite the paper above.

## Licence

Code is released under the [MIT licence](LICENSE).

## Contact

Asiri Ortega-Matienzo, Institute for Biostatistics and Informatics in Medicine and Ageing Research,
Rostock University Medical Centre. Corresponding author: Israel Barrantes (israel.barrantes@uni-rostock.de).
