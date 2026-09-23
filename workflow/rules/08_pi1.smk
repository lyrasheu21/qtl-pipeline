"""08 pi1 replication of FDR-significant cciQTLs in an external eQTL resource."""


rule pi1:
    input:
        sending=f"{RES}/qtl/main/sending.permute.annotated.tsv.gz",
        receiving=f"{RES}/qtl/main/receiving.permute.annotated.tsv.gz",
    output:
        pi1=f"{RES}/pi1/pi1.tsv",
        matched=f"{RES}/pi1/cciqtl_eqtl_matched.tsv.gz",
        scatter=f"{RES}/pi1/beta_scatter.pdf",
    params:
        rdir=R_DIR,
        settings=lambda wc: config["pi1"],  # keeps {subclass}/{chr} in eqtl_pattern literal
    resources:
        mem_mb=16000,
        runtime=6 * 60,
    log:
        f"{LOGS}/pi1.log",
    conda:
        "../../envs/r.yaml"
    script:
        "../scripts/pi1.R"
