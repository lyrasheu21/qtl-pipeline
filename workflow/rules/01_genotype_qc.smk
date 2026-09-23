"""01 Genotype QC: plink filtering, haplotype-missingness test, bgzipped VCF."""


rule genotype_qc:
    input:
        unpack(genotype_inputs),
    output:
        vcf=f"{RES}/genotype/qc.vcf.gz",
        tbi=f"{RES}/genotype/qc.vcf.gz.tbi",
        summary=f"{RES}/genotype/qc_summary.tsv",
    params:
        script=f"{SCRIPTS}/genotype_qc.sh",
        in_type="bfile" if GENO.get("bfile") else "vcf",
        in_path=GENO.get("bfile") or GENO.get("vcf"),
        keep=GENO.get("keep") or "",
        update_ids=GENO.get("update_ids") or "",
        chroms=",".join(str(c) for c in config["chromosomes"]),
        qc=config["genotype_qc"],
        prefix=lambda wc, output: output.vcf.removesuffix(".vcf.gz"),
    threads: 4
    resources:
        mem_mb=16000,
        runtime=240,
    log:
        f"{LOGS}/genotype_qc.log",
    conda:
        "../../envs/genotype.yaml"
    shell:
        """
        bash {params.script} \
            --in-type {params.in_type} --in "{params.in_path}" \
            --keep "{params.keep}" --update-ids "{params.update_ids}" \
            --chr {params.chroms} \
            --geno {params.qc[geno]} --maf {params.qc[maf]} \
            --hwe {params.qc[hwe]} --mishap-p {params.qc[mishap_p]} \
            --threads {threads} --out {params.prefix} > {log} 2>&1
        """


rule vcf_samples:
    input:
        f"{RES}/genotype/qc.vcf.gz",
    output:
        f"{RES}/genotype/samples.txt",
    log:
        f"{LOGS}/vcf_samples.log",
    conda:
        "../../envs/genotype.yaml"
    shell:
        "bcftools query -l {input} > {output} 2> {log}"
