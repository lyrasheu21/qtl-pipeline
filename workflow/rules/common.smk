"""Shared constants, config lookups, input functions and target lists."""

import hashlib
from pathlib import Path

import pandas as pd

REPO = Path(workflow.basedir).parent
QTLPOST_R = str(REPO / "qtlpost" / "R")
SCRIPTS = str(Path(workflow.basedir) / "scripts")
RES = config.get("results_dir", "results")
LOGS = f"{RES}/logs"
INPUTS = config.get("inputs", {})
CHUNKS = list(range(1, int(config["qtl"]["n_chunks"]) + 1))


def _md5_of_dir(path):
    h = hashlib.md5()
    for f in sorted(Path(path).glob("*.R")):
        h.update(f.read_bytes())
    return h.hexdigest()


# Passed as a param to the QC and downstream R rules: editing qtlpost/R
# changes it, so Snakemake re-runs them. Deliberately not used by the
# genotype / phenotype / covariate rules, whose outputs feed the expensive
# QTLtools jobs: a code change there must not force QTL mapping to rerun.
QTLPOST_MD5 = _md5_of_dir(QTLPOST_R)


wildcard_constraints:
    set=r"[A-Za-z0-9_]+",
    chunk=r"\d+",
    fraction=r"[0-9.]+",
    variant=r"[A-Za-z0-9_.]+",
    covset=r"[A-Za-z0-9_.]+",
    gwas=r"[A-Za-z0-9_]+",


def module(name):
    return bool(config.get("modules", {}).get(name, False))


def fmt_num(x):
    """Canonical text form of a threshold, used in file names (0.25 -> '0.25', 0 -> '0')."""
    return f"{float(x):g}"


# ---------------------------------------------------------------------------
# Phenotype sets. Each phenotypes entry produces one set (bed, matrix) or two
# (cci: <name>_sending, <name>_receiving).
# ---------------------------------------------------------------------------
ADAPTERS = ("bed", "matrix", "cci")


def build_sets():
    sets = {}
    for name, spec in (config.get("phenotypes") or {}).items():
        adapter = spec.get("adapter")
        if adapter not in ADAPTERS:
            raise ValueError(f"phenotypes.{name}: adapter must be one of {ADAPTERS}, got {adapter!r}")
        if adapter == "cci":
            for direction in ("sending", "receiving"):
                sets[f"{name}_{direction}"] = dict(entry=name, adapter=adapter, direction=direction, spec=spec)
        else:
            sets[name] = dict(entry=name, adapter=adapter, direction=None, spec=spec)
    if not sets:
        raise ValueError("No phenotypes configured: add at least one entry under `phenotypes:`")
    return sets


SETS = build_sets()


def set_spec(s):
    return SETS[s]["spec"]


def set_fraction(s):
    return fmt_num(set_spec(s).get("min_nonzero_fraction", 0))


def set_window(s):
    return int(set_spec(s).get("window", config["qtl"]["window"]))


def set_grouped(s):
    return bool(set_spec(s).get("grouped_permutation", False))


def grp_flag(wildcards):
    return "--grp-best" if set_grouped(wildcards.set) else ""


def r_params(**extra):
    """Params shared by every R rule: where qtlpost lives and its checksum."""
    return dict(rdir=QTLPOST_R, qtlpost_md5=QTLPOST_MD5, **extra)


def set_fdr_group(s):
    default = "group" if SETS[s]["adapter"] == "cci" else "none"
    return set_spec(s).get("fdr_group_by", default)


# ---------------------------------------------------------------------------
# Covariate sets: metadata covariates + genotype PCs + hidden factors (k).
# ---------------------------------------------------------------------------
COV = config["covariates"]
HIDDEN = COV.get("hidden_factors", {}) or {}
SENS = config.get("sensitivity", {}) or {}


def build_covsets():
    main = dict(use=list(COV["metadata"].get("use", [])), genotype_pcs=int(COV.get("genotype_pcs", 0)),
                hidden_k=int(HIDDEN.get("k", 0)))
    covsets = {"main": main}
    if module("sensitivity"):
        for name, use in (SENS.get("covariate_sets") or {}).items():
            covsets[f"cov_{name}"] = dict(main, use=list(use))
        for k in SENS.get("sweep_k") or []:
            if int(k) != main["hidden_k"]:  # the main run already covers its own K
                covsets[f"k{int(k)}"] = dict(main, hidden_k=int(k))
    return covsets


COVSETS = build_covsets()


# ---------------------------------------------------------------------------
# QTL model variants: "main" plus one per sensitivity entry.
# fraction=None means "the set's own min_nonzero_fraction".
# ---------------------------------------------------------------------------
def build_variants():
    main_perm = config["qtl"]["n_permutations"]
    variants = {"main": dict(fraction=None, covset="main", donor_fraction=None, n_perm=main_perm)}
    if module("sensitivity"):
        n = SENS.get("n_permutations", main_perm)
        for f in SENS.get("min_nonzero_fraction") or []:
            variants[f"nonzero{fmt_num(f)}"] = dict(fraction=fmt_num(f), covset="main", donor_fraction=None, n_perm=n)
        for f in SENS.get("donor_fraction") or []:
            variants[f"donors{fmt_num(f)}"] = dict(fraction=None, covset="main", donor_fraction=float(f), n_perm=n)
        for name in COVSETS:
            if name != "main":
                variants[name] = dict(fraction=None, covset=name, donor_fraction=None, n_perm=n)
    return variants


VARIANTS = build_variants()


def variant_fraction(variant, s):
    return VARIANTS[variant]["fraction"] or set_fraction(s)


# ---------------------------------------------------------------------------
# Input functions
# ---------------------------------------------------------------------------
GENO = config["genotypes"]
if bool(GENO.get("bfile")) == bool(GENO.get("vcf")):
    raise ValueError("Set exactly one of genotypes.bfile or genotypes.vcf")
PCS_SOURCE = (GENO.get("pcs") or {}).get("source", "none")


def genotype_inputs(wildcards):
    if GENO.get("bfile"):
        return {ext: f"{GENO['bfile']}.{ext}" for ext in ("bed", "bim", "fam")}
    return {"vcf": GENO["vcf"]}


def adapter_inputs(wildcards):
    """Input files of the adapter behind a set (so Snakemake tracks them)."""
    info = SETS[wildcards.set]
    spec = info["spec"]
    files = {"vcf_samples": f"{RES}/genotype/samples.txt"}
    keys = {"bed": ["bed", "annotation"], "matrix": ["matrix", "coordinates"], "cci": ["matrix"]}[info["adapter"]]
    for key in keys:
        if spec.get(key):
            files[key] = spec[key]
    if info["adapter"] == "cci":
        files["gene_annotation"] = INPUTS["gene_annotation"]
    if (spec.get("sample_map") or {}).get("file"):
        files["sample_map"] = spec["sample_map"]["file"]
    return files


def set_bed(wildcards, variant="main"):
    variant = getattr(wildcards, "variant", variant)
    f = variant_fraction(variant, wildcards.set)
    return {"bed": f"{RES}/phenotypes/{wildcards.set}.c{f}.bed.gz",
            "tbi": f"{RES}/phenotypes/{wildcards.set}.c{f}.bed.gz.tbi"}


def covariate_inputs(wildcards):
    files = dict(set_bed(wildcards, "main"))
    if COV["metadata"].get("use") or SENS.get("covariate_sets"):
        files["metadata"] = INPUTS["metadata"]
    if PCS_SOURCE != "none":
        files["pcs"] = f"{RES}/genotype/pcs.tsv"
    if HIDDEN.get("source") == "file":
        files["factors"] = HIDDEN["file"]
        if (HIDDEN.get("sample_map") or {}).get("file"):
            files["factor_map"] = HIDDEN["sample_map"]["file"]
    return files


def variant_cov(wildcards):
    return f"{RES}/covariates/{wildcards.set}.{VARIANTS[wildcards.variant]['covset']}.txt"


def variant_exclude(wildcards):
    if VARIANTS[wildcards.variant]["donor_fraction"] is None:
        return []
    return f"{RES}/samples/{wildcards.set}.exclude_{wildcards.variant}.txt"


# ---------------------------------------------------------------------------
# Downstream inputs and checks
# ---------------------------------------------------------------------------
def read_gwas_table():
    path = INPUTS.get("gwas_table", "")
    if not module("coloc") or not path:
        return pd.DataFrame()
    df = pd.read_csv(path, sep="\t", comment="#", dtype={"n_col": str, "s": float})
    return df.set_index("name", drop=False)


GWAS = read_gwas_table()
GWAS_NAMES = GWAS.index.tolist()
MAGMA_GWAS = list((config.get("magma", {}) or {}).get("gene_results", {}).keys()) if module("magma") else []

if module("coloc") and not GWAS_NAMES:
    raise ValueError(f"modules.coloc is on but no GWAS are listed in {INPUTS.get('gwas_table')}")
if module("magma") and (not MAGMA_GWAS or not INPUTS.get("gene_annotation")):
    raise ValueError("modules.magma needs magma.gene_results and inputs.gene_annotation")
if module("expression_cor") and not INPUTS.get("expression_pattern"):
    raise ValueError("modules.expression_cor needs inputs.expression_pattern")
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
        f"{RES}/qtl_qc/fdr_summary.tsv",
        f"{RES}/qtl_qc/fdr_dotplot.pdf",
        f"{RES}/qtl_qc/permute_qq.pdf",
        f"{RES}/qtl_qc/pvalue_histograms.pdf",
        f"{RES}/qtl_qc/shared_lead_snps.tsv",
        f"{RES}/qtl_qc/nominal_lambda.tsv",
        f"{RES}/figures/upset_genes_by_class.pdf",
    ]
    for s in SETS:
        t += [
            f"{RES}/phenotypes/{s}.annotation.tsv",
            f"{RES}/phenotypes/{s}.c{set_fraction(s)}.report.tsv",
            f"{RES}/covariates/{s}.main.report.tsv",
            f"{RES}/qtl/main/{s}.permute.tsv.gz",
            f"{RES}/qtl/main/{s}.nominal.fdr.tsv.gz",
        ]
    if module("sensitivity"):
        t.append(f"{RES}/qtl_qc/sensitivity_summary.tsv")
    if module("coloc"):
        t += [f"{RES}/coloc/coloc_master.tsv", f"{RES}/coloc/coloc_gold_unique.tsv",
              f"{RES}/figures/coloc_summary.pdf"]
    if module("magma"):
        t += expand(f"{RES}/magma/heatmap.{{threshold}}.{{level}}.pdf",
                    threshold=["fdr", "nominal"], level=config["magma"]["levels"])
    if module("pi1"):
        t += [f"{RES}/pi1/pi1.tsv", f"{RES}/pi1/beta_scatter.pdf"]
    if module("expression_cor"):
        t += [f"{RES}/expression_cor/phenotype_expression_cor.tsv",
              f"{RES}/expression_cor/phenotype_expression_cor.pdf"]
    return t
