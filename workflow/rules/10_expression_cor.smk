"""10 Correlation of CCI scores with ligand / receptor expression, per subclass."""


rule expression_cor_subclass:
    input:
        cci=f"{RES}/phenotype/cci_matrix.rds",
    output:
        temp(f"{RES}/expression_cor/by_subclass/{{subclass}}.tsv"),
    params:
        rdir=R_DIR,
        expression=lambda wc: INPUTS["expression_pattern"].replace("{subclass}", wc.subclass),
        min_donors=config["expression_cor"]["min_donors"],
    resources:
        mem_mb=16000,
    log:
        f"{LOGS}/expression_cor/{{subclass}}.log",
    conda:
        "../../envs/r.yaml"
    script:
        "../scripts/expression_cor_subclass.R"


rule expression_cor_merge:
    input:
        expand(f"{RES}/expression_cor/by_subclass/{{subclass}}.tsv", subclass=SUBCLASSES),
    output:
        table=f"{RES}/expression_cor/cci_expression_cor.tsv",
        plot=f"{RES}/expression_cor/cci_expression_cor.pdf",
    params:
        rdir=R_DIR,
    log:
        f"{LOGS}/expression_cor/merge.log",
    conda:
        "../../envs/r.yaml"
    script:
        "../scripts/expression_cor_merge.R"
