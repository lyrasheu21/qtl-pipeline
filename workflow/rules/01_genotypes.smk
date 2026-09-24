"""01 Genotypes: donor subset, optional plink QC, bgzipped VCF, optional genotype PCs."""


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
        qc=GENO["qc"],
        qc_flag=1 if GENO["qc"].get("enabled", True) else 0,
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
            --chr {params.chroms} --qc {params.qc_flag} \
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


if PCS_SOURCE == "compute":

    rule genotype_pcs:
        """LD-pruned genotype PCs with plink2."""
        input:
            f"{RES}/genotype/qc.vcf.gz",
        output:
            f"{RES}/genotype/pcs.tsv",
        params:
            n=GENO["pcs"].get("n_compute", 10),
            prefix=lambda wc, output: output[0].removesuffix(".tsv"),
        threads: 4
        log:
            f"{LOGS}/genotype_pcs.log",
        conda:
            "../../envs/genotype.yaml"
        shell:
            """
            (
              plink2 --vcf {input} --threads {threads} --allow-extra-chr --set-missing-var-ids @:# \
                  --indep-pairwise 200 50 0.2 --out {params.prefix}.prune
              plink2 --vcf {input} --threads {threads} --allow-extra-chr --set-missing-var-ids @:# \
                  --extract {params.prefix}.prune.prune.in --pca {params.n} --out {params.prefix}
              # eigenvec header is "#IID PC1 ..." or "#FID IID PC1 ..."; keep IID as "sample"
              awk 'BEGIN {{ OFS = "\\t" }}
                   NR == 1 {{ start = ($1 == "#FID") ? 2 : 1; $start = "sample" }}
                   {{ line = $start; for (i = start + 1; i <= NF; i++) line = line OFS $i; print line }}' \
                  {params.prefix}.eigenvec > {output}
              rm -f {params.prefix}.prune.* {params.prefix}.eigenvec {params.prefix}.eigenval
            ) > {log} 2>&1
            """

elif PCS_SOURCE == "file":

    rule genotype_pcs:
        """Genotype PCs from a file (either orientation), rewritten as samples x PCs."""
        input:
            pcs=GENO["pcs"]["file"],
            samples=f"{RES}/genotype/samples.txt",
        output:
            f"{RES}/genotype/pcs.tsv",
        params:
            rdir=QTLPOST_R,
            qtlpost_md5=QTLPOST_MD5,
        log:
            f"{LOGS}/genotype_pcs.log",
        conda:
            "../../envs/r.yaml"
        script:
            "../scripts/genotype_pcs_file.R"
