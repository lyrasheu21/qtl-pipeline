# cciQTL pipeline

[![CI](https://github.com/lyrasheu21/qtl-pipeline/actions/workflows/ci.yml/badge.svg)](https://github.com/lyrasheu21/qtl-pipeline/actions/workflows/ci.yml)

A reproducible Snakemake pipeline for mapping **cell–cell interaction QTLs (cciQTLs)**: genetic
variants that change how strongly two brain cell types communicate through a ligand–receptor
pair. It also tests whether those variants explain disease risk.

The input is a donor × interaction score matrix from single-nucleus RNA-seq cell–cell
communication inference, for example `Astro:APOE->Micro:TREM2`. The pipeline:

1. **QCs genotypes** (plink: missingness, MAF, HWE, haplotype-based missingness).
2. **Places each interaction twice:** at its ligand gene (*sending*) and at its receptor gene (*receiving*).
3. **Maps cis-QTLs** with QTLtools (permutation and nominal passes), with q-values per cell subclass.
4. **Runs downstream analyses:**
   - GWAS colocalisation (coloc)
   - MAGMA gene-set enrichment
   - π1 replication in single-cell eQTLs
   - mediation through gene expression
   - summary figures

Every step is driven by one config file and runs locally or on LSF/SLURM clusters.

![Pipeline rule graph](docs/rulegraph.svg)

## Quickstart (toy data, ~1 minute)

The repository ships a small synthetic dataset (60 donors, 2 chromosomes, 21 interactions) with
two planted causal variants, so the whole pipeline can be run without real data.

```bash
git clone https://github.com/lyrasheu21/qtl-pipeline.git && cd qtl-pipeline
conda create -n snakemake -c conda-forge -c bioconda "snakemake>=8" && conda activate snakemake

# Linux with QTLtools v1.3 on your PATH:
snakemake --profile profiles/local --configfile config/test.yaml

# No QTLtools (e.g. macOS): use the bundled regression-based stand-in
snakemake --profile profiles/local --configfile config/test.yaml tests/mock.yaml
```

Results land in `results_test/`. The planted variants `rs_causal_1` and `rs_causal_2` should come
out as FDR-significant cciQTLs that colocalise with the toy GWAS, replicate in the toy eQTLs and
show mediation. Unit tests for the shared R code: `Rscript tests/run_tests.R`.

<p align="center">
  <img src="docs/example_locus_plot.png" width="420" alt="Example locus plot: cciQTL, GWAS and colocalisation tracks with gene track">
</p>

## Running on your own data

1. Fill in the input paths and parameters in [`config/config.yaml`](config/config.yaml). Every
   threshold is a config value, e.g. `genotype_qc.hwe`, `phenotype.min_completeness` and
   `qtl.n_permutations`.
2. List the GWAS for colocalisation in [`config/gwas.tsv`](config/gwas.tsv), with a column
   mapping so any summary-statistics layout works.
3. Install QTLtools (static Linux binary) and, for the MAGMA module, MAGMA. Neither is on
   conda. Point `qtl.qtltools_binary` and `magma.binary` at them. Everything else comes from the
   conda envs in [`envs/`](envs/).
4. Run `snakemake --profile profiles/local`, or `profiles/lsf` / `profiles/slurm` after setting
   your account and queue in the profile.

| Input | Format |
|---|---|
| Genotypes | plink bfile or VCF (+ optional keep list / ID update file) |
| Interaction scores | `.parquet` / `.tsv` / `.csv`: one row per interaction, one column per donor |
| Gene annotation | TSV: `gene_name gene_id chrom start end strand` |
| Donor metadata, ancestry PCs | TSV/CSV, one row per donor; covariates and transforms chosen in the config |
| GWAS summary statistics | one row per study in `config/gwas.tsv` |
| Optional | MAGMA `.genes.raw` files, an eQTL resource (per subclass × chromosome), expression residuals (per subclass) |

Downstream modules (coloc, MAGMA, π1, mediation, expression correlation, sensitivity analyses)
can be switched on or off under `modules:` in the config.

## Outputs

| Directory | Contents |
|---|---|
| `genotype/` | QC'd VCF and a per-step variant/sample count summary |
| `phenotype/` | sending/receiving BED files and a report of interactions dropped and why |
| `qtl/`, `qtl_qc/` | QTLtools results annotated with cell types and per-subclass q-values; FDR counts, QQ plots, genomic λ |
| `coloc/` | colocalisation table with Gold/Silver/Bronze tiers and locus plots |
| `magma/`, `pi1/`, `mediation/`, `expression_cor/` | downstream results and plots |
| `figures/` | UpSet plots, colocalisation summary, mediation chord diagrams |

## Repository structure

```
config/        config.yaml (all inputs + parameters), gwas.tsv, cell_types.tsv, test.yaml, JSON schema
workflow/      Snakefile, rules/ (one file per stage, 01-11), scripts/ (one R/bash/Python script per rule)
R/             shared helpers used by the scripts, covered by unit tests
envs/          conda environments (genotype tools, R)
profiles/      execution profiles: local, LSF, SLURM
tests/         testthat unit tests, synthetic toy data + generator, mock QTLtools, CI config override
docs/          rule graph and example figures
.github/       CI: unit tests + full toy-data run on every push
```

## Background: what the refactor changed

The pipeline was rebuilt from a set of exploratory analysis notebooks: interactive R and bash
chunks with hard-coded cluster paths and manual job submission. Besides turning them into a
config-driven workflow, the rewrite fixed several methodological issues:

- **q-values per cell subclass.** A grouping bug meant FDR was computed across all
  interactions at once instead of within each subclass.
- **Covariates.** The covariate file silently dropped the global expression PCs. It is now
  built from an explicit list in the config.
- **Colocalisation inputs.** The notebooks had several diverging copies of the coloc code.
  There is now one implementation, which:
  - uses QTLtools standard errors instead of back-calculating them from p-values
  - takes the QTL sample size from the data instead of hard-coding it
  - runs case-control GWAS as `type = "cc"` with their case fraction
- **Figures from data.** Summary plots that used hand-typed counts are now generated from
  pipeline outputs.
- **Tests.** Unit tests cover the helper functions, and CI runs the full pipeline on synthetic
  data with known answers.

## Notes

- **Shared helper edits:** Snakemake re-runs a rule when its own script changes, but not when a
  shared helper in `R/` changes. After editing `R/`, force the affected rules with
  `--forcerun <rule>`.
- **Mock QTLtools:** `tests/mock/QTLtools` only exists to exercise the workflow on machines
  without QTLtools. It is not a substitute for QTLtools in real analyses.
- **Matching inputs:** SNP IDs must match across the genotypes, GWAS and eQTL files, and all
  inputs must use one genome build.

## License

[MIT](LICENSE) © 2026 Lyra Sheu
