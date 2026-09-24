"""04 cis-QTL mapping with QTLtools, per phenotype set: permutation pass for
every variant (main + sensitivity) and nominal pass for the main run, each
split into chunks.

QTLtools is not on conda. Point qtl.qtltools_binary at an installed binary
(e.g. from tools/qtltools/install_qtltools.sh) or set qtl.container and run
with --sdm apptainer."""

QTLTOOLS = config["qtl"].get("qtltools_binary", "QTLtools")
QTLTOOLS_CONTAINER = config["qtl"].get("container") or None



rule qtl_permute_chunk:
    input:
        unpack(set_bed),
        vcf=f"{RES}/genotype/qc.vcf.gz",
        vcf_tbi=f"{RES}/genotype/qc.vcf.gz.tbi",
        cov=variant_cov,
        exclude=variant_exclude,
    output:
        temp(f"{RES}/qtl/{{variant}}/{{set}}/permute/chunk{{chunk}}.txt"),
    params:
        qtltools=QTLTOOLS,
        n_perm=lambda wc: VARIANTS[wc.variant]["n_perm"],
        n_chunks=config["qtl"]["n_chunks"],
        window=lambda wc: set_window(wc.set),
        seed=config["qtl"]["seed"],
        grp=grp_flag,
        exclude=lambda wc, input: f"--exclude-samples {input.exclude}" if input.exclude else "",
    resources:
        mem_mb=8000,
        runtime=24 * 60,
    log:
        f"{LOGS}/qtl/{{variant}}.{{set}}.permute.chunk{{chunk}}.log",
    container:
        QTLTOOLS_CONTAINER
    shell:
        """
        {params.qtltools} cis --vcf {input.vcf} --bed {input.bed} --cov {input.cov} {params.exclude} \
            --permute {params.n_perm} --window {params.window} --seed {params.seed} {params.grp} \
            --normal --std-err --chunk {wildcards.chunk} {params.n_chunks} \
            --out {output} > {log} 2>&1
        touch {output}  # chunks without phenotypes produce no file
        """


rule qtl_permute_merge:
    input:
        expand(f"{RES}/qtl/{{{{variant}}}}/{{{{set}}}}/permute/chunk{{chunk}}.txt", chunk=CHUNKS),
    output:
        f"{RES}/qtl/{{variant}}/{{set}}.permute.raw.txt.gz",
    log:
        f"{LOGS}/qtl/{{variant}}.{{set}}.permute.merge.log",
    conda:
        "../../envs/genotype.yaml"
    shell:
        "cat {input} | gzip -c > {output} 2> {log}"


rule qtl_nominal_chunk:
    input:
        unpack(set_bed),
        vcf=f"{RES}/genotype/qc.vcf.gz",
        vcf_tbi=f"{RES}/genotype/qc.vcf.gz.tbi",
        cov=variant_cov,
    output:
        f"{RES}/qtl/{{variant}}/{{set}}/nominal/chunk{{chunk}}.txt.gz",
    wildcard_constraints:
        variant="main",
    params:
        qtltools=QTLTOOLS,
        threshold=config["qtl"]["nominal_threshold"],
        n_chunks=config["qtl"]["n_chunks"],
        window=lambda wc: set_window(wc.set),
        seed=config["qtl"]["seed"],
        grp=grp_flag,
        tmp=lambda wc, output: output[0].removesuffix(".gz"),
    resources:
        mem_mb=8000,
        runtime=24 * 60,
    log:
        f"{LOGS}/qtl/{{variant}}.{{set}}.nominal.chunk{{chunk}}.log",
    container:
        QTLTOOLS_CONTAINER
    shell:
        """
        {params.qtltools} cis --vcf {input.vcf} --bed {input.bed} --cov {input.cov} \
            --nominal {params.threshold} --window {params.window} --seed {params.seed} {params.grp} \
            --normal --std-err --chunk {wildcards.chunk} {params.n_chunks} \
            --out {params.tmp} > {log} 2>&1
        touch {params.tmp}
        gzip -c {params.tmp} > {output} && rm {params.tmp}
        """
