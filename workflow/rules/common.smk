"""Shared constants, config lookups and target lists."""

from pathlib import Path

import pandas as pd

R_DIR = str(Path(workflow.basedir).parent / "R")
SCRIPTS = str(Path(workflow.basedir) / "scripts")
RES = config.get("results_dir", "results")
LOGS = "logs"

DIRECTIONS = ["sending", "receiving"]
CHUNKS = list(range(1, int(config["qtl"]["n_chunks"]) + 1))
INPUTS = config["inputs"]

CELL_TYPES = pd.read_csv(INPUTS["cell_types"], sep="\t", comment="#")
SUBCLASSES = CELL_TYPES["subclass"].tolist()


wildcard_constraints:
    direction="sending|receiving",
    chunk=r"\d+",
    completeness=r"[0-9.]+",
    variant=r"[A-Za-z0-9_.]+",
    covset=r"[A-Za-z0-9_.]+",
    gwas=r"[A-Za-z0-9_]+",
    subclass=r"[A-Za-z0-9_]+",


def module(name):
    return bool(config["modules"].get(name, False))


def fmt_num(x):
    """Canonical text form of a threshold, used in file names (0.25 -> '0.25', 0 -> '0')."""
    return f"{float(x):g}"


# ---------------------------------------------------------------------------
# QTL model variants. "main" is the primary analysis; sensitivity analyses
# add one variant per completeness threshold, donor fraction or covariate set.
# Each variant maps to a phenotype BED (completeness), a covariate file and an
# optional list of donors to exclude.
# ---------------------------------------------------------------------------
def build_variants():
    main_c = fmt_num(config["phenotype"]["min_completeness"])
    variants = {
        "main": dict(
            completeness=main_c,
            covset="main",
            donor_fraction=None,
            n_perm=config["qtl"]["n_permutations"],
        )
    }
    if module("sensitivity"):
        s = config["sensitivity"]
        for c in s.get("completeness", []):
            variants[f"completeness{fmt_num(c)}"] = dict(
                completeness=fmt_num(c), covset="main", donor_fraction=None, n_perm=s["n_permutations"]
            )
        for f in s.get("donor_fraction", []):
            variants[f"donors{fmt_num(f)}"] = dict(
                completeness=main_c, covset="main", donor_fraction=float(f), n_perm=s["n_permutations"]
            )
        for name in s.get("covariate_sets", {}):
            variants[f"cov_{name}"] = dict(
                completeness=main_c, covset=name, donor_fraction=None, n_perm=s["n_permutations"]
            )
    return variants


VARIANTS = build_variants()

COVARIATE_SETS = {"main": config["covariates"]["use"]}
if module("sensitivity"):
    COVARIATE_SETS.update(config["sensitivity"].get("covariate_sets", {}))


def variant(wildcards):
    return VARIANTS[wildcards.variant]



# ---------------------------------------------------------------------------
# Input functions
# ---------------------------------------------------------------------------
GENO = INPUTS["genotype"]
if bool(GENO.get("bfile")) == bool(GENO.get("vcf")):
    raise ValueError("Set exactly one of inputs.genotype.bfile or inputs.genotype.vcf")


def genotype_inputs(wildcards):
    if GENO.get("bfile"):
        return {ext: f"{GENO['bfile']}.{ext}" for ext in ("bed", "bim", "fam")}
    return {"vcf": GENO["vcf"]}


QTLTOOLS = config["qtl"].get("qtltools_binary", "QTLtools")


def variant_bed(wildcards):
    c = VARIANTS[wildcards.variant]["completeness"]
    return {
        "bed": f"{RES}/phenotype/{wildcards.direction}.c{c}.bed.gz",
        "tbi": f"{RES}/phenotype/{wildcards.direction}.c{c}.bed.gz.tbi",
    }


def variant_cov(wildcards):
    return f"{RES}/covariates/{VARIANTS[wildcards.variant]['covset']}.txt"


def variant_exclude(wildcards):
    if VARIANTS[wildcards.variant]["donor_fraction"] is None:
        return []
    return f"{RES}/samples/exclude_{wildcards.variant}.txt"


# ---------------------------------------------------------------------------
# GWAS table (coloc) and MAGMA gene results
# ---------------------------------------------------------------------------
def read_gwas_table():
    path = INPUTS.get("gwas_table", "")
    if not module("coloc") or not path:
        return pd.DataFrame()
    df = pd.read_csv(path, sep="\t", comment="#", dtype={"n_col": str, "s": float})
    return df.set_index("name", drop=False)


GWAS = read_gwas_table()
GWAS_NAMES = GWAS.index.tolist()
MAGMA_GWAS = list(config.get("magma", {}).get("gene_results", {}).keys()) if module("magma") else []

if module("coloc") and not GWAS_NAMES:
    raise ValueError(f"modules.coloc is on but no GWAS are listed in {INPUTS.get('gwas_table')}")
if module("magma") and not MAGMA_GWAS:
    raise ValueError("modules.magma is on but magma.gene_results is empty")
if module("mediation") and not module("coloc"):
    raise ValueError("modules.mediation needs modules.coloc")
for mod in ("mediation", "expression_cor"):
    if module(mod) and not INPUTS.get("expression_pattern"):
        raise ValueError(f"modules.{mod} needs inputs.expression_pattern")
if module("pi1") and not config["pi1"].get("eqtl_pattern"):
    raise ValueError("modules.pi1 needs pi1.eqtl_pattern")


def gwas_spec(name):
    """GWAS table row as a plain dict with NaN turned into None."""
    row = GWAS.loc[name].to_dict()
    return {k: (None if pd.isna(v) else v) for k, v in row.items()}


# ---------------------------------------------------------------------------
# Targets
# ---------------------------------------------------------------------------
def all_targets():
    t = [
        f"{RES}/genotype/qc.vcf.gz.tbi",
        f"{RES}/genotype/qc_summary.tsv",
        f"{RES}/samples/sample_report.tsv",
        f"{RES}/phenotype/phenotype_report.c{VARIANTS['main']['completeness']}.tsv",
        f"{RES}/qtl_qc/fdr_counts.tsv",
        f"{RES}/qtl_qc/fdr_dotplot.pdf",
        f"{RES}/qtl_qc/permute_qq.pdf",
        f"{RES}/qtl_qc/shared_snp_hits.tsv",
        f"{RES}/qtl_qc/nominal_lambda.tsv",
        f"{RES}/figures/upset_genes_by_class.pdf",
    ]
    if module("sensitivity"):
        t.append(f"{RES}/qtl_qc/sensitivity_fdr_counts.tsv")
    if module("coloc"):
        t += [
            f"{RES}/coloc/coloc_master.tsv",
            f"{RES}/coloc/coloc_gold_unique.tsv",
            f"{RES}/figures/coloc_summary.pdf",
        ]
    if module("magma"):
        t += expand(
            f"{RES}/magma/heatmap.{{threshold}}.{{level}}.pdf",
            threshold=["fdr", "nominal"],
            level=config["magma"]["levels"],
        )
    if module("pi1"):
        t += [f"{RES}/pi1/pi1.tsv", f"{RES}/pi1/beta_scatter.pdf"]
    if module("mediation"):
        t += [f"{RES}/mediation/mediation_results.tsv", f"{RES}/mediation/interaction_results.tsv"]
        t += [f"{RES}/figures/chord_{target.replace(':', '_')}.pdf" for target in config["figures"]["chord_targets"]]
    if module("expression_cor"):
        t += [f"{RES}/expression_cor/cci_expression_cor.tsv", f"{RES}/expression_cor/cci_expression_cor.pdf"]
    return t
