"""06 Downstream analyses: thin wrappers around the qtlpost R package for
batch runs (coloc across all GWAS x sets, MAGMA, pi1, phenotype-expression
correlation) and summary figures. Mediation analysis is not part of the
pipeline; see analyses/mediation/."""

GROUP_TABLE = INPUTS["group_table"] if INPUTS.get("group_table") else []



# ---- colocalisation -----------------------------------------------------------
rule coloc:
    input:
        nominal=f"{MAIN}/{{set}}.nominal.fdr.tsv.gz",
        permute=f"{MAIN}/{{set}}.permute.tsv.gz",
        donors=f"{RES}/covariates/{{set}}.main.donors.txt",
        gwas=lambda wc: GWAS.loc[wc.gwas, "path"],
        annotation=INPUTS["gene_annotation"] if INPUTS.get("gene_annotation") else [],
    output:
        summary=f"{RES}/coloc/{{gwas}}.{{set}}.summary.tsv",
        snps=f"{RES}/coloc/{{gwas}}.{{set}}.snps.tsv",
        plots=directory(f"{RES}/coloc/locus_plots/{{gwas}}.{{set}}"),
    params:
        **r_params(),
        spec=lambda wc: gwas_spec(wc.gwas),
        settings=config["coloc"],
    resources:
        mem_mb=16000,
        runtime=12 * 60,
    log:
        f"{LOGS}/coloc/{{gwas}}.{{set}}.log",
    conda:
        "../../envs/r.yaml"
    script:
        "../scripts/coloc.R"


rule coloc_collect:
    input:
        summaries=expand(f"{RES}/coloc/{{gwas}}.{{set}}.summary.tsv", gwas=GWAS_NAMES, set=list(SETS)),
        permute=expand(f"{MAIN}/{{set}}.permute.tsv.gz", set=list(SETS)),
    output:
        master=f"{RES}/coloc/coloc_master.tsv",
        counts=f"{RES}/coloc/coloc_counts.tsv",
        gold=f"{RES}/coloc/coloc_gold_unique.tsv",
    params:
        **r_params(pp_h4=config["coloc"]["pp_h4"]),
    log:
        f"{LOGS}/coloc/collect.log",
    conda:
        "../../envs/r.yaml"
    script:
        "../scripts/coloc_collect.R"


# ---- MAGMA --------------------------------------------------------------------
rule magma_gene_sets:
    input:
        tables=expand(f"{MAIN}/{{set}}.permute.tsv.gz", set=list(SETS)),
        annotation=INPUTS.get("gene_annotation", ""),
    output:
        f"{RES}/magma/sets/{{threshold}}.{{level}}.txt",
    wildcard_constraints:
        threshold="fdr|nominal",
        level="group|class|superclass",
    params:
        **r_params(gene_id_column=config["magma"]["gene_id_column"], nominal_p=config["magma"]["nominal_p"]),
    log:
        f"{LOGS}/magma/sets.{{threshold}}.{{level}}.log",
    conda:
        "../../envs/r.yaml"
    script:
        "../scripts/magma_gene_sets.R"


rule magma_run:
    input:
        sets=f"{RES}/magma/sets/{{threshold}}.{{level}}.txt",
        genes_raw=lambda wc: config["magma"]["gene_results"][wc.gwas],
    output:
        f"{RES}/magma/{{threshold}}.{{level}}/{{gwas}}.gsa.out",
    params:
        magma=config["magma"]["binary"],
        prefix=lambda wc, output: output[0].removesuffix(".gsa.out"),
    resources:
        runtime=10 * 60,
    log:
        f"{LOGS}/magma/{{threshold}}.{{level}}.{{gwas}}.log",
    shell:
        """
        if [ ! -s {input.sets} ]; then
          echo "Empty gene-set file; skipping MAGMA" > {log}
          : > {output}
          exit 0
        fi
        {params.magma} --gene-results {input.genes_raw} --set-annot {input.sets} col=1,2 \
            --out {params.prefix} > {log} 2>&1
        # MAGMA < 1.08 names the output .sets.out
        if [ ! -f {output} ] && [ -f {params.prefix}.sets.out ]; then
          mv {params.prefix}.sets.out {output}
        fi
        """


rule magma_heatmap:
    input:
        results=expand(f"{RES}/magma/{{{{threshold}}}}.{{{{level}}}}/{{gwas}}.gsa.out", gwas=MAGMA_GWAS),
        group_table=GROUP_TABLE,
    output:
        plot=f"{RES}/magma/heatmap.{{threshold}}.{{level}}.pdf",
        table=f"{RES}/magma/results.{{threshold}}.{{level}}.tsv",
    params:
        **r_params(gwas=MAGMA_GWAS),
    log:
        f"{LOGS}/magma/heatmap.{{threshold}}.{{level}}.log",
    conda:
        "../../envs/r.yaml"
    script:
        "../scripts/magma_heatmap.R"


# ---- pi1 replication ----------------------------------------------------------
rule pi1:
    input:
        tables=expand(f"{MAIN}/{{set}}.permute.tsv.gz", set=list(SETS)),
    output:
        pi1=f"{RES}/pi1/pi1.tsv",
        matched=f"{RES}/pi1/qtl_eqtl_matched.tsv.gz",
        scatter=f"{RES}/pi1/beta_scatter.pdf",
    params:
        **r_params(),
        settings=lambda wc: config["pi1"],  # a function keeps {group}/{chr} literal
    resources:
        mem_mb=16000,
        runtime=6 * 60,
    log:
        f"{LOGS}/pi1.log",
    conda:
        "../../envs/r.yaml"
    script:
        "../scripts/pi1.R"


# ---- phenotype vs expression ----------------------------------------------------
rule expression_cor:
    input:
        unpack(lambda wc: set_bed(wc, "main")),
        annotation=f"{RES}/phenotypes/{{set}}.annotation.tsv",
    output:
        temp(f"{RES}/expression_cor/by_set/{{set}}.tsv"),
    params:
        **r_params(min_donors=config["expression_cor"]["min_donors"]),
        expression_pattern=lambda wc: INPUTS.get("expression_pattern", ""),
    resources:
        mem_mb=16000,
    log:
        f"{LOGS}/expression_cor/{{set}}.log",
    conda:
        "../../envs/r.yaml"
    script:
        "../scripts/expression_cor.R"


rule expression_cor_merge:
    input:
        expand(f"{RES}/expression_cor/by_set/{{set}}.tsv", set=list(SETS)),
    output:
        table=f"{RES}/expression_cor/phenotype_expression_cor.tsv",
        plot=f"{RES}/expression_cor/phenotype_expression_cor.pdf",
    params:
        **r_params(),
    log:
        f"{LOGS}/expression_cor/merge.log",
    conda:
        "../../envs/r.yaml"
    script:
        "../scripts/expression_cor_merge.R"


# ---- figures ------------------------------------------------------------------
rule upset:
    input:
        tables=expand(f"{MAIN}/{{set}}.permute.tsv.gz", set=list(SETS)),
        group_table=GROUP_TABLE,
    output:
        f"{RES}/figures/upset_genes_by_class.pdf",
    params:
        **r_params(),
    log:
        f"{LOGS}/figures/upset.log",
    conda:
        "../../envs/r.yaml"
    script:
        "../scripts/fig_upset.R"


rule coloc_summary_figure:
    input:
        master=f"{RES}/coloc/coloc_master.tsv",
        group_table=GROUP_TABLE,
    output:
        f"{RES}/figures/coloc_summary.pdf",
    params:
        **r_params(gwas=GWAS_NAMES, pp_h4=config["coloc"]["pp_h4"]),
    log:
        f"{LOGS}/figures/coloc_summary.log",
    conda:
        "../../envs/r.yaml"
    script:
        "../scripts/fig_coloc_summary.R"
