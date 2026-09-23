"""05 QTL QC: annotate permutation results, per-subclass q-values, summary
plots, genomic inflation, and nominal results for FDR-significant CCIs."""

MAIN = f"{RES}/qtl/main"


rule annotate_permute:
    input:
        permute=f"{RES}/qtl/{{variant}}/{{direction}}.permute.txt.gz",
        cell_types=INPUTS["cell_types"],
    output:
        table=f"{RES}/qtl/{{variant}}/{{direction}}.permute.annotated.tsv.gz",
        fdr_ids=f"{RES}/qtl/{{variant}}/{{direction}}.fdr_ids.txt",
    params:
        rdir=R_DIR,
        fdr=config["qtl_qc"]["fdr"],
    log:
        f"{LOGS}/qtl/{{variant}}.{{direction}}.annotate.log",
    conda:
        "../../envs/r.yaml"
    script:
        "../scripts/annotate_permute.R"


rule qtl_summary:
    input:
        sending=f"{MAIN}/sending.permute.annotated.tsv.gz",
        receiving=f"{MAIN}/receiving.permute.annotated.tsv.gz",
        cell_types=INPUTS["cell_types"],
    output:
        counts=f"{RES}/qtl_qc/fdr_counts.tsv",
        dotplot=f"{RES}/qtl_qc/fdr_dotplot.pdf",
        qq=f"{RES}/qtl_qc/permute_qq.pdf",
        shared=f"{RES}/qtl_qc/shared_snp_hits.tsv",
    params:
        rdir=R_DIR,
    log:
        f"{LOGS}/qtl_summary.log",
    conda:
        "../../envs/r.yaml"
    script:
        "../scripts/qtl_summary.R"


rule nominal_pvalue_sample:
    """Random sample of nominal p-values (column 12) for lambda and QQ plots."""
    input:
        expand(f"{MAIN}/{{{{direction}}}}/nominal/chunk{{chunk}}.txt.gz", chunk=CHUNKS),
    output:
        f"{RES}/qtl_qc/nominal_pvalues.{{direction}}.txt",
    params:
        fraction=config["qtl_qc"]["nominal_sample_fraction"],
        seed=config["qtl_qc"].get("seed", 42),
    log:
        f"{LOGS}/nominal_pvalue_sample.{{direction}}.log",
    conda:
        "../../envs/genotype.yaml"
    shell:
        """
        for f in {input}; do gzip -dc "$f"; done \
          | awk -v frac={params.fraction} -v seed={params.seed} \
              'BEGIN {{ srand(seed) }} $1 != "phe_id" && rand() < frac {{ print $12 }}' \
          > {output} 2> {log}
        """


rule nominal_lambda:
    input:
        sending=f"{RES}/qtl_qc/nominal_pvalues.sending.txt",
        receiving=f"{RES}/qtl_qc/nominal_pvalues.receiving.txt",
    output:
        table=f"{RES}/qtl_qc/nominal_lambda.tsv",
        qq=f"{RES}/qtl_qc/nominal_qq.pdf",
    params:
        rdir=R_DIR,
    log:
        f"{LOGS}/nominal_lambda.log",
    conda:
        "../../envs/r.yaml"
    script:
        "../scripts/nominal_lambda.R"


rule extract_fdr_nominal:
    """All nominal associations of the FDR-significant CCIs (input to coloc)."""
    input:
        ids=f"{MAIN}/{{direction}}.fdr_ids.txt",
        chunks=expand(f"{MAIN}/{{{{direction}}}}/nominal/chunk{{chunk}}.txt.gz", chunk=CHUNKS),
    output:
        f"{MAIN}/{{direction}}.nominal.fdr.txt.gz",
    log:
        f"{LOGS}/extract_fdr_nominal.{{direction}}.log",
    conda:
        "../../envs/genotype.yaml"
    shell:
        """
        if [ -s {input.ids} ]; then
          for f in {input.chunks}; do gzip -dc "$f"; done \
            | awk 'NR == FNR {{ ids[$1]; next }} $1 in ids' {input.ids} - \
            | gzip -c > {output}
        else
          : | gzip -c > {output}
        fi 2> {log}
        """


rule sensitivity_summary:
    input:
        expand(
            f"{RES}/qtl/{{variant}}/{{direction}}.permute.annotated.tsv.gz",
            variant=list(VARIANTS),
            direction=DIRECTIONS,
        ),
    output:
        table=f"{RES}/qtl_qc/sensitivity_fdr_counts.tsv",
        plot=f"{RES}/qtl_qc/sensitivity_fdr_counts.pdf",
    params:
        rdir=R_DIR,
        variants=[v for v in VARIANTS for d in DIRECTIONS],
        variant_info=VARIANTS,
    log:
        f"{LOGS}/sensitivity_summary.log",
    conda:
        "../../envs/r.yaml"
    script:
        "../scripts/sensitivity_summary.R"
