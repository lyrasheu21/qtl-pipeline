"""04 cis-QTL mapping with QTLtools: permutation pass (every variant) and
nominal pass (main analysis only), each split into chunks.

QTLtools is not packaged for conda; it is called through qtl.qtltools_binary
(https://qtltools.github.io/qtltools/ provides static Linux binaries)."""



rule qtl_permute_chunk:
    input:
        unpack(variant_bed),
        vcf=f"{RES}/genotype/qc.vcf.gz",
        vcf_tbi=f"{RES}/genotype/qc.vcf.gz.tbi",
        cov=variant_cov,
        exclude=variant_exclude,
    output:
        temp(f"{RES}/qtl/{{variant}}/{{direction}}/permute/chunk{{chunk}}.txt"),
    params:
        qtltools=QTLTOOLS,
        n_perm=lambda wc: VARIANTS[wc.variant]["n_perm"],
        n_chunks=config["qtl"]["n_chunks"],
        window=config["qtl"]["window"],
        seed=config["qtl"]["seed"],
        exclude=lambda wc, input: f"--exclude-samples {input.exclude}" if input.exclude else "",
    resources:
        mem_mb=8000,
        runtime=24 * 60,
    log:
        f"{LOGS}/qtl/{{variant}}.{{direction}}.permute.chunk{{chunk}}.log",
    shell:
        """
        {params.qtltools} cis --vcf {input.vcf} --bed {input.bed} --cov {input.cov} {params.exclude} \
            --permute {params.n_perm} --window {params.window} --seed {params.seed} \
            --normal --std-err --chunk {wildcards.chunk} {params.n_chunks} \
            --out {output} > {log} 2>&1
        touch {output}  # chunks without phenotypes produce no file
        """


rule qtl_permute_merge:
    input:
        expand(
            f"{RES}/qtl/{{{{variant}}}}/{{{{direction}}}}/permute/chunk{{chunk}}.txt",
            chunk=CHUNKS,
        ),
    output:
        f"{RES}/qtl/{{variant}}/{{direction}}.permute.txt.gz",
    log:
        f"{LOGS}/qtl/{{variant}}.{{direction}}.permute.merge.log",
    conda:
        "../../envs/genotype.yaml"
    shell:
        "cat {input} | gzip -c > {output} 2> {log}"


rule qtl_nominal_chunk:
    input:
        unpack(variant_bed),
        vcf=f"{RES}/genotype/qc.vcf.gz",
        vcf_tbi=f"{RES}/genotype/qc.vcf.gz.tbi",
        cov=variant_cov,
    output:
        f"{RES}/qtl/{{variant}}/{{direction}}/nominal/chunk{{chunk}}.txt.gz",
    wildcard_constraints:
        variant="main",
    params:
        qtltools=QTLTOOLS,
        threshold=config["qtl"]["nominal_threshold"],
        n_chunks=config["qtl"]["n_chunks"],
        window=config["qtl"]["window"],
        seed=config["qtl"]["seed"],
        tmp=lambda wc, output: output[0].removesuffix(".gz"),
    resources:
        mem_mb=8000,
        runtime=24 * 60,
    log:
        f"{LOGS}/qtl/{{variant}}.{{direction}}.nominal.chunk{{chunk}}.log",
    shell:
        """
        {params.qtltools} cis --vcf {input.vcf} --bed {input.bed} --cov {input.cov} \
            --nominal {params.threshold} --window {params.window} --seed {params.seed} \
            --normal --std-err --chunk {wildcards.chunk} {params.n_chunks} \
            --out {params.tmp} > {log} 2>&1
        touch {params.tmp}
        gzip -c {params.tmp} > {output} && rm {params.tmp}
        """
