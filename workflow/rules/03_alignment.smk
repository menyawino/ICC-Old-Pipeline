# A rule to align the trimmed reads to the reference genome using bwa and convert to BAM using sambamba view
rule bwa_mem:
    message:
        "Aligning and converting to BAM for sample {wildcards.sample}_{lane}"
    input:
        fq1=rules.trimming.output.fq1,
        fq2=rules.trimming.output.fq2
    output:
        sam=config["outdir"] + "/analysis/003_alignment/01_bwa/{sample}_{lane}.sam"
    conda:
        "icc_04_alignment"
    threads:
        config["threads_high"]
    params: 
        ref=config["reference_genome"],
        bwa=config["bwa"]
    log:
        bwa=config["outdir"] + "/logs/003_alignment/01_bwa/{sample}_{lane}_bwa.log"
    benchmark:
        config["outdir"] + "/benchmarks/003_alignment/01_bwa/{sample}_{lane}_alignment.txt"
    shell:
        """
        {params.bwa} mem \
        -M \
        -t {threads} \
        {params.ref} \
        {input.fq1} \
        {input.fq2} \
        2> {log.bwa} \
        > {output.sam} \
        2> {log.bwa}
        """

rule convert_bam:
    message:
        "Converting SAM to BAM for sample {wildcards.sample}_{lane}"
    input:
        sam=rules.bwa_mem.output.sam
    output:
        bam=config["outdir"] + "/analysis/003_alignment/01_bwa/{sample}_{lane}.bam"
    conda:
        "icc_04_alignment"
    threads:
        config["threads_high"]
    log:
        samtools_view=config["outdir"] + "/logs/003_alignment/01_bwa/{sample}_{lane}_samtools_view.log",
        samtools_sort=config["outdir"] + "/logs/003_alignment/01_bwa/{sample}_{lane}_samtools_sort.log"
    benchmark:
        config["outdir"] + "/benchmarks/003_alignment/01_bwa/{sample}_{lane}_convert_bam.txt"
    params:
        samtools=config["samtools"]
    shell:
        """
        {params.samtools} view \
        -bS \
        -@ {threads} \
        {input.sam} \
        -o {output.bam}.unsorted \
        2> {log.samtools_view}
        
        {params.samtools} sort \
        {output.bam}.unsorted \
        -@ {threads} \
        -o {output.bam} \
        > {output.bam} \
        2> {log.samtools_sort}
        
        rm {output.bam}.unsorted
        """