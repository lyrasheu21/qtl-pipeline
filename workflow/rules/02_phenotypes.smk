"""02 Phenotypes: adapters turn each dataset into QTLtools BEDs plus a
phenotype annotation table (phenotype_id, set, chrom, start, end, strand,
gene, group, [partner_gene, partner_group])."""


rule make_phenotype:
    input:
        unpack(adapter_inputs),
    output:
        bed=temp(f"{RES}/phenotypes/{{set}}.c{{fraction}}.bed"),
        annotation=f"{RES}/phenotypes/{{set}}.c{{fraction}}.annotation.tsv",
        report=f"{RES}/phenotypes/{{set}}.c{{fraction}}.report.tsv",
    params:
        rdir=QTLPOST_R,
        adapter=lambda wc: SETS[wc.set]["adapter"],
        direction=lambda wc: SETS[wc.set]["direction"] or "",
        spec=lambda wc: SETS[wc.set]["spec"],
        min_fraction=lambda wc: float(wc.fraction),
        chromosomes=[str(c) for c in config["chromosomes"]],
    resources:
        mem_mb=32000,
    log:
        f"{LOGS}/make_phenotype.{{set}}.c{{fraction}}.log",
    conda:
        "../../envs/r.yaml"
    script:
        "../scripts/make_phenotype.R"


rule index_bed:
    input:
        f"{RES}/phenotypes/{{set}}.c{{fraction}}.bed",
    output:
        bed=f"{RES}/phenotypes/{{set}}.c{{fraction}}.bed.gz",
        tbi=f"{RES}/phenotypes/{{set}}.c{{fraction}}.bed.gz.tbi",
    log:
        f"{LOGS}/index_bed.{{set}}.c{{fraction}}.log",
    conda:
        "../../envs/genotype.yaml"
    shell:
        "(bgzip -c {input} > {output.bed} && tabix -f -p bed {output.bed}) > {log} 2>&1"


rule set_annotation:
    """The annotation of the set's main BED, as the standard output."""
    input:
        lambda wc: f"{RES}/phenotypes/{wc.set}.c{set_fraction(wc.set)}.annotation.tsv",
    output:
        f"{RES}/phenotypes/{{set}}.annotation.tsv",
    log:
        f"{LOGS}/set_annotation.{{set}}.log",
    conda:
        "../../envs/genotype.yaml"
    shell:
        "cp {input} {output} 2> {log}"
