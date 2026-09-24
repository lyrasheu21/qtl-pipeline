"""03 Covariates: one QTLtools covariate file per set and covariate set
(main, sensitivity covariate sets, hidden-factor sweep)."""


rule covariates:
    input:
        unpack(covariate_inputs),
    output:
        cov=f"{RES}/covariates/{{set}}.{{covset}}.txt",
        donors=f"{RES}/covariates/{{set}}.{{covset}}.donors.txt",
        report=f"{RES}/covariates/{{set}}.{{covset}}.report.tsv",
    params:
        rdir=QTLPOST_R,
        qtlpost_md5=QTLPOST_MD5,
        covset=lambda wc: COVSETS[wc.covset],
        definitions=COV["metadata"].get("definitions", {}),
        metadata_id=INPUTS.get("metadata_id_column", "SubID"),
        hidden=lambda wc: HIDDEN,
        min_donors=COV.get("min_donors", 30),
    log:
        f"{LOGS}/covariates.{{set}}.{{covset}}.log",
    conda:
        "../../envs/r.yaml"
    script:
        "../scripts/covariates.R"


rule exclude_samples:
    """Random subset of donors to drop for a donor-fraction sensitivity variant."""
    input:
        f"{RES}/covariates/{{set}}.main.donors.txt",
    output:
        f"{RES}/samples/{{set}}.exclude_{{variant}}.txt",
    params:
        keep_fraction=lambda wc: VARIANTS[wc.variant]["donor_fraction"],
        seed=SENS.get("seed", 42),
    log:
        f"{LOGS}/exclude_samples.{{set}}.{{variant}}.log",
    script:
        "../scripts/exclude_samples.py"
