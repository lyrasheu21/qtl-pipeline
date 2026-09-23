"""03 Covariates: QTLtools covariate files and donor-exclusion lists."""


rule covariates:
    input:
        donors=f"{RES}/samples/donors.txt",
        metadata=INPUTS["metadata"],
        ancestry_pcs=INPUTS["ancestry_pcs"] if INPUTS.get("ancestry_pcs") else [],
    output:
        f"{RES}/covariates/{{covset}}.txt",
    params:
        rdir=R_DIR,
        metadata_id=INPUTS["metadata_id_column"],
        pcs_id=INPUTS.get("ancestry_pcs_id_column", ""),
        definitions=config["covariates"]["definitions"],
        use=lambda wc: COVARIATE_SETS[wc.covset],
        n_pcs=config["covariates"]["n_ancestry_pcs"],
        pc_prefix=config["covariates"].get("ancestry_pc_prefix", "PC"),
    log:
        f"{LOGS}/covariates.{{covset}}.log",
    conda:
        "../../envs/r.yaml"
    script:
        "../scripts/covariates.R"


rule exclude_samples:
    """Random subset of donors to drop for a donor-fraction sensitivity variant."""
    input:
        f"{RES}/samples/donors.txt",
    output:
        f"{RES}/samples/exclude_{{variant}}.txt",
    params:
        keep_fraction=lambda wc: VARIANTS[wc.variant]["donor_fraction"],
        seed=config.get("sensitivity", {}).get("seed", 42),
    log:
        f"{LOGS}/exclude_samples.{{variant}}.log",
    script:
        "../scripts/exclude_samples.py"
