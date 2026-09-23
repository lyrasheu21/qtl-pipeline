"""02 Phenotype: CCI matrix -> donor list -> sending/receiving BED files."""


rule cci_matrix:
    input:
        INPUTS["cci_matrix"],
    output:
        f"{RES}/phenotype/cci_matrix.rds",
    params:
        rdir=R_DIR,
        id_column=INPUTS["cci_id_column"],
    resources:
        mem_mb=32000,
    log:
        f"{LOGS}/cci_matrix.log",
    conda:
        "../../envs/r.yaml"
    script:
        "../scripts/cci_matrix.R"


rule sample_list:
    input:
        vcf_samples=f"{RES}/genotype/samples.txt",
        cci=f"{RES}/phenotype/cci_matrix.rds",
        metadata=INPUTS["metadata"],
    output:
        donors=f"{RES}/samples/donors.txt",
        report=f"{RES}/samples/sample_report.tsv",
    params:
        rdir=R_DIR,
        metadata_id=INPUTS["metadata_id_column"],
        min_donors=config["phenotype"]["min_donors"],
    log:
        f"{LOGS}/sample_list.log",
    conda:
        "../../envs/r.yaml"
    script:
        "../scripts/sample_list.R"


rule make_bed:
    input:
        cci=f"{RES}/phenotype/cci_matrix.rds",
        donors=f"{RES}/samples/donors.txt",
        annotation=INPUTS["gene_annotation"],
        cell_types=INPUTS["cell_types"],
    output:
        sending=temp(f"{RES}/phenotype/sending.c{{completeness}}.bed"),
        receiving=temp(f"{RES}/phenotype/receiving.c{{completeness}}.bed"),
        report=f"{RES}/phenotype/phenotype_report.c{{completeness}}.tsv",
    params:
        rdir=R_DIR,
        min_completeness=lambda wc: float(wc.completeness),
        zero_as_missing=config["phenotype"]["zero_as_missing"],
        chromosomes=[str(c) for c in config["chromosomes"]],
    resources:
        mem_mb=32000,
    log:
        f"{LOGS}/make_bed.c{{completeness}}.log",
    conda:
        "../../envs/r.yaml"
    script:
        "../scripts/make_bed.R"


rule index_bed:
    input:
        f"{RES}/phenotype/{{direction}}.c{{completeness}}.bed",
    output:
        bed=f"{RES}/phenotype/{{direction}}.c{{completeness}}.bed.gz",
        tbi=f"{RES}/phenotype/{{direction}}.c{{completeness}}.bed.gz.tbi",
    log:
        f"{LOGS}/index_bed.{{direction}}.c{{completeness}}.log",
    conda:
        "../../envs/genotype.yaml"
    shell:
        "(bgzip -c {input} > {output.bed} && tabix -f -p bed {output.bed}) > {log} 2>&1"
