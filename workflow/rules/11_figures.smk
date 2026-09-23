"""11 Summary figures."""


rule upset:
    input:
        sending=f"{RES}/qtl/main/sending.permute.annotated.tsv.gz",
        receiving=f"{RES}/qtl/main/receiving.permute.annotated.tsv.gz",
        cell_types=INPUTS["cell_types"],
    output:
        f"{RES}/figures/upset_genes_by_class.pdf",
    params:
        rdir=R_DIR,
    log:
        f"{LOGS}/figures/upset.log",
    conda:
        "../../envs/r.yaml"
    script:
        "../scripts/fig_upset.R"


rule coloc_summary_figure:
    input:
        master=f"{RES}/coloc/coloc_master.tsv",
        cell_types=INPUTS["cell_types"],
    output:
        f"{RES}/figures/coloc_summary.pdf",
    params:
        rdir=R_DIR,
        gwas=GWAS_NAMES,
        pp_h4=config["coloc"]["pp_h4"],
    log:
        f"{LOGS}/figures/coloc_summary.log",
    conda:
        "../../envs/r.yaml"
    script:
        "../scripts/fig_coloc_summary.R"


rule chord:
    input:
        mediation=f"{RES}/mediation/mediation_results.tsv",
        cell_types=INPUTS["cell_types"],
    output:
        f"{RES}/figures/chord_{{target}}.pdf",
    wildcard_constraints:
        target=r"[A-Za-z0-9_]+_[A-Za-z0-9.\-]+",
    params:
        rdir=R_DIR,
        target=lambda wc: next(t for t in config["figures"]["chord_targets"] if t.replace(":", "_") == wc.target),
    log:
        f"{LOGS}/figures/chord_{{target}}.log",
    conda:
        "../../envs/r.yaml"
    script:
        "../scripts/fig_chord.R"
