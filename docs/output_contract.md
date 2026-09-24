# Output contract

The pipeline's job ends at a fixed set of files under `results_dir`. The
`qtlpost` R package (and any other downstream tool) only relies on these. A
*set* is one phenotype table mapped as a unit: one per `bed` / `matrix`
entry, two per `cci` entry (`<name>_sending`, `<name>_receiving`).

## Per set

| File | Contents |
|---|---|
| `phenotypes/{set}.annotation.tsv` | One row per phenotype (all phenotypes, before filtering): `phenotype_id set chrom start end strand gene group`, plus `partner_gene partner_group` for sets with an interaction partner (CCI) and any adapter-specific columns (`ligand`, `receptor`, `sending_group`, `receiving_group` for CCI). |
| `phenotypes/{set}.c{fraction}.bed.gz` (+ `.tbi`) | The QTLtools phenotype BED actually mapped (donors = genotyped donors after any sample map). |
| `phenotypes/{set}.c{fraction}.report.tsv` | Phenotypes dropped by the min_nonzero_fraction filter, missing coordinates or chromosome, and the number kept. |
| `covariates/{set}.main.txt` | QTLtools covariate file (`id` + one column per donor). |
| `covariates/{set}.main.donors.txt`, `.report.tsv` | Donors with complete covariates (the QTL sample size) and what was dropped. |
| `qtl/main/{set}.permute.tsv.gz` | **Main results**: one row per phenotype (or per group with `--grp-best`). |
| `qtl/main/{set}.nominal.fdr.tsv.gz` | All nominal associations of the FDR-significant phenotypes (input to colocalisation). |
| `qtl/main/{set}/nominal/chunk*.txt.gz` | Raw QTLtools nominal output, all phenotypes. |

### `qtl/main/{set}.permute.tsv.gz` columns

| Column | Meaning |
|---|---|
| `set`, `phenotype_id` | Phenotype set and ID (with `--grp-best`: the best phenotype in the group) |
| `gene`, `group`, `class`, `superclass` | From the annotation and the group table |
| `grp_id`, `n_phe_in_grp` | Group ID and size for `--grp-best` runs, otherwise NA |
| `var_id`, `var_chr`, `var_from`, `dist_phe_var` | Lead variant |
| `slope`, `slope_se`, `nom_pval`, `r_squared` | Lead variant association (QTLtools `--std-err`) |
| `adj_emp_pval` | Permutation p-value (QTLtools column 19) |
| `adj_beta_pval` | Beta-approximated permutation p-value (column 20); q-values use this |
| `qval`, `fdr_pass` | q-value within the set's `fdr_group_by` strata; `qval <= qtl_qc.fdr` |
| `partner_gene`, `partner_group`, ... | Adapter-specific annotation columns, when present |

### `qtl/main/{set}.nominal.fdr.tsv.gz` columns

`phe_id phe_chr phe_from phe_to phe_strd grp_id n_phe_in_grp n_var_in_cis
dist_phe_var var_id var_chr var_from var_to nom_pval r_squared slope slope_se
best_hit` (QTLtools `cis --nominal --std-err` output with names;
`qtlpost::read_qtltools()` reads the raw files into the same layout).

## Across sets

| File | Contents |
|---|---|
| `genotype/qc.vcf.gz`, `qc_summary.tsv` | Genotypes used for mapping and variants/samples removed per QC step |
| `genotype/pcs.tsv` | Genotype PCs (samples x PCs), when `genotypes.pcs.source` is `file` or `compute` |
| `qtl_qc/fdr_summary.tsv` | Per set and group: phenotypes tested, pi0 (with a warning note), QTLs at FDR 1/5/10% |
| `qtl_qc/nominal_lambda.tsv` | Genomic inflation per set from sampled nominal p-values |
| `qtl_qc/sensitivity_summary.tsv` | FDR counts for every sensitivity variant, including the hidden-factor sweep (`k`) |
| `qtl_qc/*.pdf` | Dot plot, QQ plots, p-value histograms, nominal QQ, sensitivity / K-sweep plots |
| `coloc/coloc_master.tsv` | Colocalised phenotypes (PP.H4 >= threshold) with tier, lead SNP and GWAS / QTL statistics |
| `pi1/pi1.tsv`, `magma/*`, `expression_cor/*`, `figures/*` | Downstream module results |
