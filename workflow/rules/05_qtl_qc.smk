"""05 QTL QC and standard outputs: annotated permutation results with
q-values, diagnostics, genomic inflation, nominal results for FDR hits and
sensitivity / hidden-factor sweep summaries."""

MAIN = f"{RES}/qtl/main"


rule annotate_permute:
    input:
        permute=f"{RES}/qtl/{{variant}}/{{set}}.permute.raw.txt.gz",
        annotation=f"{RES}/phenotypes/{{set}}.annotation.tsv",
        group_table=INPUTS["group_table"] if INPUTS.get("group_table") else [],
    output:
        table=f"{RES}/qtl/{{variant}}/{{set}}.permute.tsv.gz",
        fdr_keys=f"{RES}/qtl/{{variant}}/{{set}}.fdr_keys.txt",
    params:
        rdir=QTLPOST_R,
        qtlpost_md5=QTLPOST_MD5,
        fdr=config["qtl_qc"]["fdr"],
        fdr_group_by=lambda wc: set_fdr_group(wc.set),
    log:
        f"{LOGS}/qtl/{{variant}}.{{set}}.annotate.log",
    conda:
        "../../envs/r.yaml"
    script:
        "../scripts/annotate_permute.R"


rule qtl_summary:
    input:
        tables=expand(f"{MAIN}/{{set}}.permute.tsv.gz", set=list(SETS)),
        group_table=INPUTS["group_table"] if INPUTS.get("group_table") else [],
    output:
        summary=f"{RES}/qtl_qc/fdr_summary.tsv",
        dotplot=f"{RES}/qtl_qc/fdr_dotplot.pdf",
        qq=f"{RES}/qtl_qc/permute_qq.pdf",
        hist=f"{RES}/qtl_qc/pvalue_histograms.pdf",
        shared=f"{RES}/qtl_qc/shared_lead_snps.tsv",
    params:
        rdir=QTLPOST_R,
        qtlpost_md5=QTLPOST_MD5,
    log:
        f"{LOGS}/qtl_summary.log",
    conda:
        "../../envs/r.yaml"
    script:
        "../scripts/qtl_summary.R"


rule nominal_pvalue_sample:
    """Random sample of nominal p-values for lambda and QQ plots (column 12,
    or 14 for grouped runs, which have two extra columns)."""
    input:
        expand(f"{MAIN}/{{{{set}}}}/nominal/chunk{{chunk}}.txt.gz", chunk=CHUNKS),
    output:
        f"{RES}/qtl_qc/nominal_pvalues.{{set}}.txt",
    params:
        fraction=config["qtl_qc"]["nominal_sample_fraction"],
        seed=config["qtl_qc"].get("seed", 42),
        col=lambda wc: 14 if set_grouped(wc.set) else 12,
    log:
        f"{LOGS}/nominal_pvalue_sample.{{set}}.log",
    conda:
        "../../envs/genotype.yaml"
    shell:
        """
        for f in {input}; do gzip -dc "$f"; done \
          | awk -v frac={params.fraction} -v seed={params.seed} -v col={params.col} \
              'BEGIN {{ srand(seed) }} $1 != "phe_id" && $1 != "grp_id" && rand() < frac {{ print $col }}' \
          > {output} 2> {log}
        """


rule nominal_lambda:
    input:
        expand(f"{RES}/qtl_qc/nominal_pvalues.{{set}}.txt", set=list(SETS)),
    output:
        table=f"{RES}/qtl_qc/nominal_lambda.tsv",
        qq=f"{RES}/qtl_qc/nominal_qq.pdf",
    params:
        rdir=QTLPOST_R,
        qtlpost_md5=QTLPOST_MD5,
        sets=list(SETS),
    log:
        f"{LOGS}/nominal_lambda.log",
    conda:
        "../../envs/r.yaml"
    script:
        "../scripts/nominal_lambda.R"


rule extract_fdr_nominal:
    """All nominal associations of the FDR-significant phenotypes, with named
    columns (the input to colocalisation)."""
    input:
        fdr_keys=f"{MAIN}/{{set}}.fdr_keys.txt",
        chunks=expand(f"{MAIN}/{{{{set}}}}/nominal/chunk{{chunk}}.txt.gz", chunk=CHUNKS),
    output:
        raw=temp(f"{MAIN}/{{set}}.nominal.fdr.raw.txt"),
        table=f"{MAIN}/{{set}}.nominal.fdr.tsv.gz",
    params:
        rdir=QTLPOST_R,
        qtlpost_md5=QTLPOST_MD5,
    log:
        f"{LOGS}/extract_fdr_nominal.{{set}}.log",
    conda:
        "../../envs/r.yaml"
    script:
        "../scripts/extract_fdr_nominal.R"


rule sensitivity_summary:
    input:
        expand(f"{RES}/qtl/{{variant}}/{{set}}.permute.tsv.gz", variant=list(VARIANTS), set=list(SETS)),
    output:
        table=f"{RES}/qtl_qc/sensitivity_summary.tsv",
        plot=f"{RES}/qtl_qc/sensitivity_summary.pdf",
    params:
        rdir=QTLPOST_R,
        qtlpost_md5=QTLPOST_MD5,
        variants=[v for v in VARIANTS for s in SETS],
        sets=[s for v in VARIANTS for s in SETS],
        variant_info=VARIANTS,
        covsets=COVSETS,
    log:
        f"{LOGS}/sensitivity_summary.log",
    conda:
        "../../envs/r.yaml"
    script:
        "../scripts/sensitivity_summary.R"
