"""06 Colocalisation of FDR-significant cciQTLs with each GWAS (coloc.abf)."""


rule coloc:
    input:
        nominal=f"{RES}/qtl/main/{{direction}}.nominal.fdr.txt.gz",
        donors=f"{RES}/samples/donors.txt",
        gwas=lambda wc: GWAS.loc[wc.gwas, "path"],
        annotation=INPUTS["gene_annotation"],
    output:
        summary=f"{RES}/coloc/{{gwas}}.{{direction}}.summary.tsv",
        snps=f"{RES}/coloc/{{gwas}}.{{direction}}.snps.tsv",
        plots=directory(f"{RES}/coloc/locus_plots/{{gwas}}.{{direction}}"),
    params:
        rdir=R_DIR,
        spec=lambda wc: gwas_spec(wc.gwas),
        settings=config["coloc"],
    resources:
        mem_mb=16000,
        runtime=12 * 60,
    log:
        f"{LOGS}/coloc/{{gwas}}.{{direction}}.log",
    conda:
        "../../envs/r.yaml"
    script:
        "../scripts/coloc.R"


rule coloc_collect:
    input:
        summaries=expand(f"{RES}/coloc/{{gwas}}.{{direction}}.summary.tsv", gwas=GWAS_NAMES, direction=DIRECTIONS),
    output:
        master=f"{RES}/coloc/coloc_master.tsv",
        counts=f"{RES}/coloc/coloc_counts.tsv",
        gold=f"{RES}/coloc/coloc_gold_unique.tsv",
    params:
        rdir=R_DIR,
        pp_h4=config["coloc"]["pp_h4"],
    log:
        f"{LOGS}/coloc/collect.log",
    conda:
        "../../envs/r.yaml"
    script:
        "../scripts/coloc_collect.R"
