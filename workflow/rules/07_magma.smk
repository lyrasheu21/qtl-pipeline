"""07 MAGMA gene-set analysis of cciQTL genes against GWAS gene-level results."""


rule magma_gene_sets:
    input:
        sending=f"{RES}/qtl/main/sending.permute.annotated.tsv.gz",
        receiving=f"{RES}/qtl/main/receiving.permute.annotated.tsv.gz",
        annotation=INPUTS["gene_annotation"],
    output:
        f"{RES}/magma/sets/{{threshold}}.{{level}}.txt",
    wildcard_constraints:
        threshold="fdr|nominal",
        level="subclass|class|superclass",
    params:
        rdir=R_DIR,
        gene_id_column=config["magma"]["gene_id_column"],
        nominal_p=config["magma"]["nominal_p"],
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
        expand(f"{RES}/magma/{{{{threshold}}}}.{{{{level}}}}/{{gwas}}.gsa.out", gwas=MAGMA_GWAS),
    output:
        plot=f"{RES}/magma/heatmap.{{threshold}}.{{level}}.pdf",
        table=f"{RES}/magma/results.{{threshold}}.{{level}}.tsv",
    params:
        rdir=R_DIR,
        gwas=MAGMA_GWAS,
        cell_types=INPUTS["cell_types"],
    log:
        f"{LOGS}/magma/heatmap.{{threshold}}.{{level}}.log",
    conda:
        "../../envs/r.yaml"
    script:
        "../scripts/magma_heatmap.R"
