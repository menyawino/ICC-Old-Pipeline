rule bam_sort:
    message:
        "Sorting BAM file for sample {wildcards.sample} on lane {wildcards.lane}"
    input:
        bam=rules.convert_bam.output.bam
    output:
        bam=config["outdir"] + "/analysis/003_alignment/01_bwa/{sample}_{lane}.sorted.bam"
    conda:
        "../envs/java.yml"
    threads:
        config["threads_high"]
    log:
        bam_sort=config["outdir"] + "/logs/003_alignment/01_bwa/{sample}_{lane}_bam_sort.log"
    benchmark:
        config["outdir"] + "/benchmarks/003_alignment/01_bwa/{sample}_{lane}_bam_sort.txt"
    params:
        picard=config["picard"]
    shell:
        """
        java -jar \
        {params.picard}/SortSam.jar \
        INPUT={input.bam} \
        OUTPUT={output.bam} \
        SORT_ORDER=coordinate \
        CREATE_INDEX=true \
        VALIDATION_STRINGENCY=SILENT \
        2> {log.bam_sort}
        """

rule add_read_groups:
    message:
        "Adding read groups to BAM file for sample {wildcards.sample} on lane {wildcards.lane}"
    input:
        bam=rules.bam_sort.output.bam
    output:
        bam_rg=config["outdir"] + "/analysis/003_alignment/01_bwa/{sample}_{lane}.rg.bam"
    conda:
        "../envs/java.yml"
    threads:
        config["threads_high"]
    log:
        picard_add_rg=config["outdir"] + "/logs/003_alignment/01_bwa/{sample}_{lane}_picard_add_rg.log"
    benchmark:
        config["outdir"] + "/benchmarks/003_alignment/01_bwa/{sample}_{lane}_add_rg.txt"
    params:
        picard=config["picard"],
        rgid=config["gatk_params"]["AddOrReplaceReadGroups"]["RGID"],
        rglb=config["gatk_params"]["AddOrReplaceReadGroups"]["RGLB"],
        rgpl=config["gatk_params"]["AddOrReplaceReadGroups"]["RGPL"],
        rgpu=config["gatk_params"]["AddOrReplaceReadGroups"]["RGPU"],
        rgsm=config["gatk_params"]["AddOrReplaceReadGroups"]["RGSM"],
        rgcn=config["gatk_params"]["AddOrReplaceReadGroups"]["RGCN"],
        rgds=config["gatk_params"]["AddOrReplaceReadGroups"]["RGDS"],
        validation_stringency=config["gatk_params"]["AddOrReplaceReadGroups"]["validation_stringency"]
    shell:
        """
        java -jar \
        {params.picard}/AddOrReplaceReadGroups.jar \
        I={input.bam} \
        O={output.bam_rg} \
        SORT_ORDER=coordinate \
        RGID={params.rgid} \
        RGLB={params.rglb} \
        RGPL={params.rgpl} \
        RGPU={params.rgpu} \
        RGSM={params.rgsm} \
        RGCN={params.rgcn} \
        RGDS={params.rgds} \
        VALIDATION_STRINGENCY={params.validation_stringency} \
        2> {log.picard_add_rg}
        """

rule mark_duplicates:
    message:
        "Marking duplicates in BAM file for sample {wildcards.sample} on lane {wildcards.lane}"
    input:
        bam= rules.add_read_groups.output.bam_rg
    output:
        bam_md=config["outdir"] + "/analysis/003_alignment/01_bwa/{sample}_{lane}.md.bam"
    conda:
        "../envs/java.yml"
    threads:
        config["threads_high"]
    log:
        picard_markdup=config["outdir"] + "/logs/003_alignment/01_bwa/{sample}_{lane}_picard_markdup.log"
    benchmark:
        config["outdir"] + "/benchmarks/003_alignment/01_bwa/{sample}_{lane}_markdup.txt"
    params:
        picard=config["picard"],
    shell:
        """
        java -jar \
        {params.picard}/MarkDuplicates.jar \
        I={input.bam} \
        O={output.bam_md} \
        CREATE_INDEX=true \
        VALIDATION_STRINGENCY=SILENT \
        M={log.picard_markdup} \
        2> {log.picard_markdup}
        """

rule realigner_target_creator:
    message:
        "Creating realignment targets for sample {wildcards.sample} on lane {wildcards.lane}"
    input:
        bam=rules.mark_duplicates.output.bam_md
    output:
        intervals=config["outdir"] + "/analysis/003_alignment/01_bwa/{sample}_{lane}.md.bam.intervals"
    conda:
        "../envs/java.yml"
    threads:
        config["threads_high"]
    log:
        realigner_target_creator=config["outdir"] + "/logs/003_alignment/01_bwa/{sample}_{lane}_realigner_target_creator.log"
    benchmark:
        config["outdir"] + "/benchmarks/003_alignment/01_bwa/{sample}_{lane}_realigner_target_creator.txt"
    params:
        gatk=config["gatk"],
        ref=config["reference_genome"],
        TargetFile=config["icc_panel"],
        mills=config["mills"],
        tenkindel=config["tenk_indel"]
    shell:
        """
        java -jar \
        {params.gatk} \
        -T RealignerTargetCreator \
        -nt {threads} \
        -I {input.bam} \
        -R {params.ref} \
        -L {params.TargetFile} \
        --out {output.intervals} \
        -known {params.mills} \
        -known {params.tenkindel} \
        2> {log.realigner_target_creator}
        """

rule indel_realigner:
    message:
        "Performing indel realignment for sample {wildcards.sample} on lane {wildcards.lane}"
    input:
        bam=rules.mark_duplicates.output.bam_md,
        intervals=rules.realigner_target_creator.output.intervals
    output:
        bam_realigned=config["outdir"] + "/analysis/003_alignment/01_bwa/{sample}_{lane}.md.Realigned.bam"
    conda:
        "../envs/java.yml"
    threads:
        config["threads_high"]
    log:
        indel_realigner=config["outdir"] + "/logs/003_alignment/01_bwa/{sample}_{lane}_indel_realigner.log"
    benchmark:
        config["outdir"] + "/benchmarks/003_alignment/01_bwa/{sample}_{lane}_indel_realigner.txt"
    params:
        gatk=config["gatk"],
        ref=config["reference_genome"],
        mills=config["mills"],
        tenkindel=config["tenk_indel"]
    shell:
        """
        java -jar \
        {params.gatk} \
        -I {input.bam} \
        -R {params.ref} \
        -T IndelRealigner \
        -targetIntervals {input.intervals} \
        --out {output.bam_realigned} \
        -known {params.mills} \
        -known {params.tenkindel} \
        -model USE_READS \
        2> {log.indel_realigner}
        """

rule base_recalibrator:
    message:
        "Performing base recalibration for sample {wildcards.sample} on lane {wildcards.lane}"
    input:
        bam=rules.indel_realigner.output.bam_realigned
    output:
        recal_table=config["outdir"] + "/analysis/003_alignment/01_bwa/{sample}_{lane}.md.Realigned.bam.recalibTable.grp"
    conda:
        "../envs/java.yml"
    threads:
        config["threads_high"]
    log:
        base_recal=config["outdir"] + "/logs/003_alignment/01_bwa/{sample}_{lane}_base_recal.log"
    benchmark:
        config["outdir"] + "/benchmarks/003_alignment/01_bwa/{sample}_{lane}_base_recal.txt"
    params:
        gatk=config["gatk"],
        ref=config["reference_genome"],
        TargetFile=config["icc_panel"],
        dbSNP=config["dbsnp"],
        mills=config["mills"],
        tenkindel=config["tenk_indel"]
    shell:
        """
        java -jar \
        {params.gatk} \
        -R {params.ref} \
        -I {input.bam} \
        --knownSites {params.dbSNP} \
        --knownSites {params.mills} \
        --knownSites {params.tenkindel} \
        -L {params.TargetFile} \
        -T BaseRecalibrator \
        -nct {threads} \
        -cov ReadGroupCovariate \
        -cov QualityScoreCovariate \
        -cov CycleCovariate \
        -cov ContextCovariate \
        --out {output.recal_table} \
        2> {log.base_recal}
        """

rule print_reads:
    message:
        "Printing reads for sample {wildcards.sample} on lane {wildcards.lane}"
    input:
        bam=rules.indel_realigner.output.bam_realigned,
        recal_table=rules.base_recalibrator.output.recal_table
    output:
        bam_recal=config["outdir"] + "/analysis/003_alignment/01_bwa/{sample}_{lane}.md.Realigned.recalibrated.bam"
    conda:
        "../envs/java.yml"
    threads:
        config["threads_high"]
    log:
        print_reads=config["outdir"] + "/logs/003_alignment/01_bwa/{sample}_{lane}_print_reads.log"
    benchmark:
        config["outdir"] + "/benchmarks/003_alignment/01_bwa/{sample}_{lane}_print_reads.txt"
    params:
        gatk=config["gatk"],
        ref=config["reference_genome"]
    shell:
        """
        java -jar \
        {params.gatk} \
        -R {params.ref} \
        -I {input.bam} \
        -T PrintReads \
        -nct {threads} \
        --out {output.bam_recal} \
        -BQSR {input.recal_table} \
        -S LENIENT \
        2> {log.print_reads}
        """

rule table_recalibration:
    message:
        "Performing table recalibration for sample {wildcards.sample} on lane {wildcards.lane}"
    input:
        bam=rules.indel_realigner.output.bam_realigned,
        recal_table=rules.base_recalibrator.output.recal_table
    output:
        bam_recal=config["outdir"] + "/analysis/003_alignment/01_bwa/{sample}_{lane}.md.Realigned.recalibrated.bam"
    conda:
        "../envs/java.yml"
    threads:
        config["threads_high"]
    log:
        config["outdir"] + "/logs/003_alignment/01_bwa/{sample}_{lane}_table_recalibration.log"
    benchmark:
        config["outdir"] + "/benchmarks/003_alignment/01_bwa/{sample}_{lane}_table_recalibration.txt"
    params:
        gatk=config["gatk"],
        ref=config["reference_genome"]
    shell:
        """
        java -jar \
        {params.gatk} \
        -R {params.ref} \
        -I {input.bam} \
        -T PrintReads \
        -nct {threads} \
        --out {output.bam_recal} \
        -BQSR {input.recal_table} \
        -S LENIENT \
        2> {log}
        """

wildcard_constraints:
    sample = r"[0-9A-Z]+/[0-9A-Z]+_S\d+"

rule merge_samfiles:
    message:
        "Merging BAM files for sample {wildcards.sample}"
    input:
        l1=config["outdir"] + "/analysis/003_alignment/01_bwa/{sample}_L001.bam",
        # l2=config["outdir"] + "/analysis/003_alignment/01_bwa/{sample}_L002.bam",
        # l3=config["outdir"] + "/analysis/003_alignment/01_bwa/{sample}_L003.bam",
        # l4=config["outdir"] + "/analysis/003_alignment/01_bwa/{sample}_L004.bam"
    output:
        merged_bam=config["outdir"] + "/analysis/003_alignment/01_bwa/{sample}.bam"
    conda:
        "../envs/java.yml"
    threads:
        config["threads_high"]
    log:
        config["outdir"] + "/logs/003_alignment/01_bwa/{sample}_merge_samfiles.log"
    benchmark:
        config["outdir"] + "/benchmarks/003_alignment/01_bwa/{sample}_merge_samfiles.txt"
    params:
        picard=config["picard"]
    shell:
        """
        java -jar \
        {params.picard}/MergeSamFiles.jar \
        INPUT={input.l1} \
        OUTPUT={output.merged_bam} \
        CREATE_INDEX=TRUE \
        VALIDATION_STRINGENCY=SILENT \
        2> {log}
        """

rule mark_duplicates_merged:
    message:
        "Marking duplicates in merged BAM file for sample {wildcards.sample}"
    input:
        bam=rules.merge_samfiles.output.merged_bam
    output:
        bam_md=config["outdir"] + "/analysis/003_alignment/01_bwa/{sample}.markDup.bam"
    conda:
        "../envs/java.yml"
    threads:
        config["threads_high"]
    log:
        markdup_log=config["outdir"] + "/logs/003_alignment/01_bwa/{sample}_markdup_merged.log"
    benchmark:
        config["outdir"] + "/benchmarks/003_alignment/01_bwa/{sample}_markdup_merged.txt"
    params:
        picard=config["picard"]
    shell:
        """
        java -jar \
        {params.picard}/MarkDuplicates.jar \
        INPUT={input.bam} \
        OUTPUT={output.bam_md} \
        METRICS_FILE={output.bam_md}.metrics.txt \
        VALIDATION_STRINGENCY=LENIENT \
        2> {log.markdup_log}
        """

rule add_read_groups_merged:
    message:
        "Adding read groups to merged BAM file for sample {wildcards.sample}"
    input:
        bam=rules.mark_duplicates_merged.output.bam_md
    output:
        bam_rg_merged=config["outdir"] + "/analysis/003_alignment/01_bwa/{sample}.rg.bam"
    conda:
        "../envs/java.yml"
    threads:
        config["threads_high"]
    log:
        picard_add_rg_merged=config["outdir"] + "/logs/003_alignment/01_bwa/{sample}_picard_add_rg_merged.log"
    benchmark:
        config["outdir"] + "/benchmarks/003_alignment/01_bwa/{sample}_add_rg_merged.txt"
    params:
        picard=config["picard"],
        rgid=config["gatk_params"]["AddOrReplaceReadGroups"]["RGID"],
        rglb=config["gatk_params"]["AddOrReplaceReadGroups"]["RGLB"],
        rgpl=config["gatk_params"]["AddOrReplaceReadGroups"]["RGPL"],
        rgpu=config["gatk_params"]["AddOrReplaceReadGroups"]["RGPU"],
        rgsm=config["gatk_params"]["AddOrReplaceReadGroups"]["RGSM"],
        rgcn=config["gatk_params"]["AddOrReplaceReadGroups"]["RGCN"],
        rgds=config["gatk_params"]["AddOrReplaceReadGroups"]["RGDS"],
        validation_stringency=config["gatk_params"]["AddOrReplaceReadGroups"]["validation_stringency"]
    shell:
        """
        java -jar \
        {params.picard}/AddOrReplaceReadGroups.jar \
        I={input.bam} \
        O={output.bam_rg_merged} \
        SORT_ORDER=coordinate \
        RGID={params.rgid} \
        RGLB={params.rglb} \
        RGPL={params.rgpl} \
        RGPU={params.rgpu} \
        RGSM={params.rgsm} \
        RGCN={params.rgcn} \
        RGDS={params.rgds} \
        VALIDATION_STRINGENCY={params.validation_stringency} \
        2> {log.picard_add_rg_merged}
        """

rule index_merged_markdup:
    message:
        "Indexing merged BAM file for sample {wildcards.sample}"
    input:
        bam_md=rules.add_read_groups_merged.output.bam_rg_merged
    output:
        bam_bai=config["outdir"] + "/analysis/003_alignment/01_bwa/{sample}.rg.bam.bai"
    conda:
        "../envs/java.yml"
    threads:
        config["threads_high"]
    log:
        index_log=config["outdir"] + "/logs/003_alignment/01_bwa/{sample}_markdup_index.log"
    benchmark:
        config["outdir"] + "/benchmarks/003_alignment/01_bwa/{sample}_markdup_index.txt"
    params:
        samtools=config["samtools"]
    shell:
        """
        {params.samtools} index {input.bam_md} 2> {log.index_log}
        """

rule realigner_target_creator_merged:
    message:
        "Creating realignment targets for merged BAM file for sample {wildcards.sample}"
    input:
        bam_md=rules.add_read_groups_merged.output.bam_rg_merged,
        bam_bai=rules.index_merged_markdup.output.bam_bai
    output:
        intervals=config["outdir"] + "/analysis/003_alignment/01_bwa/{sample}.markDup.bam.intervals"
    conda:
        "../envs/java.yml"
    threads:
        config["threads_high"]
    log:
        realigner_log=config["outdir"] + "/logs/003_alignment/01_bwa/{sample}_markdup_realigner_target_creator.log"
    benchmark:
        config["outdir"] + "/benchmarks/003_alignment/01_bwa/{sample}_markdup_realigner_target_creator.txt"
    params:
        gatk=config["gatk"],
        ref=config["reference_genome"],
        TargetFile=config["icc_panel"],
        mills=config["mills"],
        tenkindel=config["tenk_indel"]
    shell:
        """
        java -jar \
        {params.gatk} \
        -T RealignerTargetCreator \
        -nt {threads} \
        -I {input.bam_md} \
        -R {params.ref} \
        -L {params.TargetFile} \
        --out {output.intervals} \
        -known {params.mills} \
        -known {params.tenkindel} \
        2> {log.realigner_log}
        """

rule indel_realigner_merged:
    message:
        "Performing indel realignment for merged BAM file for sample {wildcards.sample}"
    input:
        bam_md=rules.add_read_groups_merged.output.bam_rg_merged,
        intervals=rules.realigner_target_creator_merged.output.intervals
    output:
        bam_realigned=config["outdir"] + "/analysis/003_alignment/01_bwa/{sample}.markDup.Realigned.recalibrated.bam"
    conda:
        "../envs/java.yml"
    threads:
        config["threads_high"]
    log:
        indel_realigner_log=config["outdir"] + "/logs/003_alignment/01_bwa/{sample}_markdup_indel_realigner.log"
    benchmark:
        config["outdir"] + "/benchmarks/003_alignment/01_bwa/{sample}_markdup_indel_realigner.txt"
    params:
        gatk=config["gatk"],
        samtools=config["samtools"],
        ref=config["reference_genome"],
        mills=config["mills"],
        tenkindel=config["tenk_indel"]
    shell:
        """
        java -jar \
        {params.gatk} \
        -I {input.bam_md} \
        -R {params.ref} \
        -T IndelRealigner \
        -targetIntervals {input.intervals} \
        --out {output.bam_realigned} \
        -known {params.mills} \
        -known {params.tenkindel} \
        -model USE_READS \
        2> {log.indel_realigner_log}

        {params.samtools} index {output.bam_realigned}
        """

rule on_target_map_qual8:
    message:
        "Filtering on target reads with MAPQ >= 8 for sample {wildcards.sample}"
    input:
        bam=rules.indel_realigner_merged.output.bam_realigned
    output:
        out_bam=config["outdir"] + "/analysis/003_alignment/01_bwa/{sample}_onTarget.bam"
    conda:
        "../envs/java.yml"
    threads:
        config["threads_high"]
    log:
        out_log=config["outdir"] + "/logs/003_alignment/01_bwa/{sample}_onTargetMapQ8.log"
    benchmark:
        config["outdir"] + "/benchmarks/003_alignment/01_bwa/{sample}_onTargetMapQ8.txt"
    params:
        samtools=config["samtools"],
        bedtools=config["bedtools"],
        targetFile=config["icc_panel"]
    shell:
        """
        {params.samtools} view -uq 8 {input.bam} | \
        {params.bedtools}/intersectBed -abam stdin -b {params.targetFile} -u > {output.out_bam}

        {params.samtools} index {output.out_bam} 2> {log.out_log}
        """

rule on_prot_coding_target_map_qual8:
    message:
        "Filtering on prot coding target reads with MAPQ >= 8 for sample {wildcards.sample}"
    input:
        bam=rules.indel_realigner_merged.output.bam_realigned
    output:
        out_bam=config["outdir"] + "/analysis/003_alignment/01_bwa/{sample}_onProtCodingTarget.bam"
    conda:
        "../envs/java.yml"
    threads:
        config["threads_high"]
    log:
        out_log=config["outdir"] + "/logs/003_alignment/01_bwa/{sample}_onProtCodingTargetMapQ8.log"
    benchmark:
        config["outdir"] + "/benchmarks/003_alignment/01_bwa/{sample}_onProtCodingTargetMapQ8.txt"
    params:
        samtools=config["samtools"],
        bedtools=config["bedtools"],
        cdsFile=config["cds_panel"]
    shell:
        """
        {params.samtools} view -uq 8 {input.bam} | \
        {params.bedtools}/intersectBed -abam stdin -b {params.cdsFile} -u > {output.out_bam}

        {params.samtools} index {output.out_bam} 2> {log.out_log}
        """

rule on_canonical_target_map_qual8:
    message:
        "Filtering on canonical target reads with MAPQ >= 8 for sample {wildcards.sample}"
    input:
        bam=rules.indel_realigner_merged.output.bam_realigned
    output:
        out_bam=config["outdir"] + "/analysis/003_alignment/01_bwa/{sample}_onCanonicalTarget.bam"
    conda:
        "../envs/java.yml"
    threads:
        config["threads_high"]
    log:
        out_log=config["outdir"] + "/logs/003_alignment/01_bwa/{sample}_onCanonicalTargetMapQ8.log"
    benchmark:
        config["outdir"] + "/benchmarks/003_alignment/01_bwa/{sample}_onCanonicalTargetMapQ8.txt"
    params:
        samtools=config["samtools"],
        bedtools=config["bedtools"],
        canonTranFile=config["canontran_panel"]
    shell:
        """
        {params.samtools} view -uq 8 {input.bam} | \
        {params.bedtools}/intersectBed -abam stdin -b {params.canonTranFile} -u > {output.out_bam}

        {params.samtools} index {output.out_bam} 2> {log.out_log}
        """
