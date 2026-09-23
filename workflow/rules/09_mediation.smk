"""09 Mediation for strong colocalised cciQTLs:
genotype -> anchor gene expression -> partner gene expression."""


rule mediation_select:
    input:
        f"{RES}/coloc/coloc_gold_unique.tsv",
    output:
        hits=f"{RES}/mediation/hits.tsv",
        snps=f"{RES}/mediation/snps.txt",
    params:
        rdir=R_DIR,
        settings=config["mediation"],
    log:
        f"{LOGS}/mediation/select.log",
    conda:
        "../../envs/r.yaml"
    script:
        "../scripts/mediation_select.R"


rule mediation_genotypes:
    input:
        snps=f"{RES}/mediation/snps.txt",
        vcf=f"{RES}/genotype/qc.vcf.gz",
        tbi=f"{RES}/genotype/qc.vcf.gz.tbi",
    output:
        f"{RES}/mediation/genotypes.tsv",
    log:
        f"{LOGS}/mediation/genotypes.log",
    conda:
        "../../envs/genotype.yaml"
    shell:
        """
        (
          printf "CHROM\\tPOS\\tID\\tREF\\tALT"
          bcftools query -l {input.vcf} | awk '{{ printf "\\t%s", $1 }}'
          printf "\\n"
          if [ -s {input.snps} ]; then
            bcftools view -i 'ID=@{input.snps}' {input.vcf} \
              | bcftools query -f '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT[\\t%GT]\\n'
          fi
        ) > {output} 2> {log}
        """


rule mediation_run:
    input:
        hits=f"{RES}/mediation/hits.tsv",
        genotypes=f"{RES}/mediation/genotypes.tsv",
    output:
        mediation=f"{RES}/mediation/mediation_results.tsv",
        interaction=f"{RES}/mediation/interaction_results.tsv",
    params:
        rdir=R_DIR,
        # a function, so Snakemake does not try to fill {subclass} itself
        expression_pattern=lambda wc: INPUTS.get("expression_pattern", ""),
        sims=config["mediation"]["sims"],
        seed=config["mediation"]["seed"],
    resources:
        mem_mb=16000,
        runtime=12 * 60,
    log:
        f"{LOGS}/mediation/run.log",
    conda:
        "../../envs/r.yaml"
    script:
        "../scripts/mediation_run.R"
