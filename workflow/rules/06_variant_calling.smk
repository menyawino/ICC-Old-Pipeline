rule unified_genotyper_call:
    message:
        "Running Unified Genotyper on sample {wildcards.sample} on target regions"
    input:
        bam=rules.on_target_map_qual8.output.out_bam
    output:
        snp_indel_vcf=config["outdir"] + "/analysis/004_variant_call/01_unified_genotyper/{sample}_onTarget.vcf",
        snp_indel_metrics=config["outdir"] + "/analysis/004_variant_call/01_unified_genotyper/{sample}_onTarget.metrics"
    conda:
        "../envs/java.yml"
    threads:
        config["threads_high"]
    log:
        out_log=config["outdir"] + "/logs/004_variant_call/01_unified_genotyper/{sample}_onTarget.log"
    benchmark:
        config["outdir"] + "/benchmarks/004_variant_call/01_unified_genotyper/{sample}_onTarget.txt"
    params:
        gatk=config["gatk"],
        ref=config["reference_genome"],
        dbSNP=config["dbsnp"],
        targetFile=config["icc_panel"],
        dcovg=config["gatk_params"]["HaplotypeCaller"]["dcovg"]
    shell:
        """
        java -jar \
        {params.gatk} \
        -l INFO \
        -R {params.ref} \
        -L {params.targetFile} \
        -I {input.bam} \
        -T UnifiedGenotyper \
        -baq CALCULATE_AS_NECESSARY \
        -minIndelCnt 4 \
        -glm BOTH \
        --dbsnp {params.dbSNP} \
        -o {output.snp_indel_vcf} \
        -A Coverage \
        -A AlleleBalance \
        -G Standard \
        --min_base_quality_score 10 \
        --metrics_file {output.snp_indel_metrics} \
        --num_threads {threads} \
        -stand_call_conf 30.0 \
        -stand_emit_conf 10 \
        --downsample_to_coverage {params.dcovg} \
        --downsampling_type BY_SAMPLE \
        --computeSLOD \
        > {log.out_log} 2>&1
        """

rule unified_genotyper_select_indel:
    message:
        "Selecting INDEL variants for sample {wildcards.sample}"
    input:
        snp_indel_vcf=rules.unified_genotyper_call.output.snp_indel_vcf
    output:
        indel_vcf=config["outdir"] + "/analysis/004_variant_call/01_unified_genotyper/{sample}_onTarget.indel.vcf"
    conda:
        "../envs/java.yml"
    log:
        out_log=config["outdir"] + "/logs/004_variant_call/01_unified_genotyper/{sample}_select_indel.log"
    benchmark:
        config["outdir"] + "/benchmarks/004_variant_call/01_unified_genotyper/{sample}_select_indel.txt"
    params:
        gatk=config["gatk"],
        ref=config["reference_genome"],
        dcovg=config["gatk_params"]["HaplotypeCaller"]["dcovg"],
        nt=config["threads_high"]
    shell:
        """
        java -jar \
        {params.gatk} \
        -T SelectVariants \
        --downsample_to_coverage {params.dcovg} \
        --downsampling_type BY_SAMPLE \
        -R {params.ref} \
        -nt {params.nt} \
        --variant {input.snp_indel_vcf} \
        --selectTypeToInclude INDEL \
        -o {output.indel_vcf} \
        > {log.out_log} 2>&1
        """

rule unified_genotyper_filter_indel:
    message:
        "Filtering INDEL variants for sample {wildcards.sample}"
    input:
        indel_vcf=rules.unified_genotyper_select_indel.output.indel_vcf
    output:
        indel_filtered_vcf=config["outdir"] + "/analysis/004_variant_call/01_unified_genotyper/{sample}_onTarget.indel.filtered.vcf"
    conda:
        "../envs/java.yml"
    log:
        out_log=config["outdir"] + "/logs/004_variant_call/01_unified_genotyper/{sample}_filter_indel.log"
    benchmark:
        config["outdir"] + "/benchmarks/004_variant_call/01_unified_genotyper/{sample}_filter_indel.txt"
    params:
        gatk=config["gatk"],
        ref=config["reference_genome"],
        dcovg=config["gatk_params"]["HaplotypeCaller"]["dcovg"]
    shell:
        """
        java -jar \
        {params.gatk} \
        -T VariantFiltration \
        --downsample_to_coverage {params.dcovg} \
        --downsampling_type BY_SAMPLE \
        -R {params.ref} \
        --variant {input.indel_vcf} \
        -o {output.indel_filtered_vcf} \
        --filterExpression "QD < 2.0" \
        --filterExpression "ReadPosRankSum < -20.0" \
        --filterExpression "FS > 200.0" \
        --filterName QDFilter \
        --filterName ReadPosFilter \
        --filterName FSFilter \
        > {log.out_log} 2>&1
        """

rule unified_genotyper_select_snp:
    message:
        "Selecting SNP variants for sample {wildcards.sample}"
    input:
        snp_indel_vcf=rules.unified_genotyper_call.output.snp_indel_vcf
    output:
        snp_vcf=config["outdir"] + "/analysis/004_variant_call/01_unified_genotyper/{sample}_onTarget.snp.vcf"
    conda:
        "../envs/java.yml"
    log:
        out_log=config["outdir"] + "/logs/004_variant_call/01_unified_genotyper/{sample}_select_snp.log"
    benchmark:
        config["outdir"] + "/benchmarks/004_variant_call/01_unified_genotyper/{sample}_select_snp.txt"
    params:
        gatk=config["gatk"],
        ref=config["reference_genome"],
        dcovg=config["gatk_params"]["HaplotypeCaller"]["dcovg"],
        nt=config["threads_high"]
    shell:
        """
        java -jar \
        {params.gatk} \
        -T SelectVariants \
        -R {params.ref} \
        --downsample_to_coverage {params.dcovg} \
        --downsampling_type BY_SAMPLE \
        -nt {params.nt} \
        --variant {input.snp_indel_vcf} \
        --selectTypeToInclude SNP \
        -o {output.snp_vcf} \
        > {log.out_log} 2>&1
        """

rule unified_genotyper_filter_snp:
    message:
        "Filtering SNP variants for sample {wildcards.sample}"
    input:
        snp_vcf=rules.unified_genotyper_select_snp.output.snp_vcf,
        indel_vcf=rules.unified_genotyper_select_indel.output.indel_vcf
    output:
        snp_filtered_vcf=config["outdir"] + "/analysis/004_variant_call/01_unified_genotyper/{sample}_onTarget.snp.filtered.vcf"
    conda:
        "../envs/java.yml"
    log:
        out_log=config["outdir"] + "/logs/004_variant_call/01_unified_genotyper/{sample}_filter_snp.log"
    benchmark:
        config["outdir"] + "/benchmarks/004_variant_call/01_unified_genotyper/{sample}_filter_snp.txt"
    params:
        gatk=config["gatk"],
        ref=config["reference_genome"],
        dcovg=config["gatk_params"]["HaplotypeCaller"]["dcovg"]
    shell:
        """
        java -jar \
        {params.gatk} \
        -T VariantFiltration \
        -R {params.ref} \
        --downsample_to_coverage {params.dcovg} \
        --downsampling_type BY_SAMPLE \
        --variant {input.snp_vcf} \
        --mask {input.indel_vcf} \
        --maskName InDel \
        -o {output.snp_filtered_vcf} \
        --filterExpression "QD < 2.0" \
        --filterExpression "MQ < 40.0" \
        --filterExpression "FS > 60.0" \
        --filterExpression "MQRankSum < -12.5" \
        --filterExpression "ReadPosRankSum < -8.0" \
        --filterName QDFilter \
        --filterName MQFilter \
        --filterName FSFilter \
        --filterName MQRankSumFilter \
        --filterName ReadPosFilter \
        > {log.out_log} 2>&1
        """

rule haplotype_caller:
    message:
        "Running HaplotypeCaller for sample {wildcards.sample}"
    input:
        bam=rules.on_target_map_qual8.output.out_bam
    output:
        snp_indel_vcf=config["outdir"] + "/analysis/004_variant_call/02_haplotype_caller/{sample}_onTarget.vcf",
        bam_out=config["outdir"] + "/analysis/004_variant_call/02_haplotype_caller/{sample}_onTarget.bam"
    log:
        config["outdir"] + "/logs/004_variant_call/02_haplotype_caller/{sample}_onTarget.log"
    conda:
        "../envs/java.yml"
    threads:
        config["threads_high"]
    benchmark:
        config["outdir"] + "/benchmarks/004_variant_call/02_haplotype_caller/{sample}_onTarget.txt"
    params:
        gatk=config["gatk"],
        ref=config["reference_genome"],
        targetFile=config["icc_panel"],
        dbSNP=config["dbsnp"],
        dcovg=config["gatk_params"]["HaplotypeCaller"]["dcovg"]
    shell:
        """
        java -jar \
        {params.gatk} \
        -R {params.ref} \
        -L {params.targetFile} \
        -I {input.bam} \
        -T HaplotypeCaller \
        --downsample_to_coverage {params.dcovg} \
        --downsampling_type BY_SAMPLE \
        -A Coverage \
        -A AlleleBalance \
        -G Standard \
        --dbsnp {params.dbSNP} \
        --min_base_quality_score 10 \
        -o {output.snp_indel_vcf} \
        --bamOutput {output.bam_out} \
        -stand_call_conf 30 \
        -stand_emit_conf 10 \
        > {log} 2>&1
        """

rule haplotype_variant_annotator:
    message:
        "Annotating variants for sample {wildcards.sample}"
    input:
        bam=rules.on_target_map_qual8.output.out_bam,
        snp_indel_vcf=rules.haplotype_caller.output.snp_indel_vcf
    output:
        annotated_vcf=config["outdir"] + "/analysis/004_variant_call/02_haplotype_caller/{sample}_onTarget.annotated.vcf"
    log:
        out_log=config["outdir"] + "/logs/004_variant_call/02_haplotype_caller/{sample}_onTarget.annotate.log"
    conda:
        "../envs/java.yml"
    benchmark:
        config["outdir"] + "/benchmarks/004_variant_call/02_haplotype_caller/{sample}_annotate.txt"
    params:
        gatk=config["gatk"],
        ref=config["reference_genome"],
        targetFile=config["icc_panel"],
        dbSNP=config["dbsnp"],
        threads=config["threads_high"],
        dcovg=config["gatk_params"]["HaplotypeCaller"]["dcovg"]
    shell:
        """
        java -jar \
        {params.gatk} \
        -R {params.ref} \
        -T VariantAnnotator \
        --downsample_to_coverage {params.dcovg} \
        --downsampling_type BY_SAMPLE \
        -L {params.targetFile} \
        -I {input.bam} \
        --variant {input.snp_indel_vcf} \
        -o {output.annotated_vcf} \
        -A Coverage \
        -A AlleleBalance \
        -A HaplotypeScore \
        -A InbreedingCoeff \
        -A HomopolymerRun \
        -A HardyWeinberg \
        -A GCContent \
        --dbsnp {params.dbSNP} \
        -nt {params.threads} \
        > {log.out_log} 2>&1
        """

rule haplotype_select_indel:
    message:
        "Selecting INDEL variants for sample {wildcards.sample}"
    input:
        annotated_vcf=rules.haplotype_variant_annotator.output.annotated_vcf
    output:
        indel_vcf=config["outdir"] + "/analysis/004_variant_call/02_haplotype_caller/{sample}_onTarget.indel.vcf"
    log:
        out_log=config["outdir"] + "/logs/004_variant_call/02_haplotype_caller/{sample}_onTarget.select_indel.log"
    conda:
        "../envs/java.yml"
    benchmark:
        config["outdir"] + "/benchmarks/004_variant_call/02_haplotype_caller/{sample}_select_indel.txt"
    params:
        gatk=config["gatk"],
        ref=config["reference_genome"],
        dcovg=config["gatk_params"]["HaplotypeCaller"]["dcovg"],
        threads=config["threads_high"]
    shell:
        """
        java -jar \
        {params.gatk} \
        -T SelectVariants \
        -R {params.ref} \
        --downsample_to_coverage {params.dcovg} \
        --downsampling_type BY_SAMPLE \
        -nt {params.threads} \
        --variant {input.annotated_vcf} \
        --selectTypeToInclude INDEL \
        -o {output.indel_vcf} \
        > {log.out_log} 2>&1
        """

rule haplotype_filter_indel:
    message:
        "Filtering INDEL variants for sample {wildcards.sample}"
    input:
        indel_vcf=rules.haplotype_select_indel.output.indel_vcf
    output:
        filtered_indel_vcf=config["outdir"] + "/analysis/004_variant_call/02_haplotype_caller/{sample}_onTarget.indel.filtered.vcf"
    log:
        out_log=config["outdir"] + "/logs/004_variant_call/02_haplotype_caller/{sample}_onTarget.filter_indel.log"
    conda:
        "../envs/java.yml"
    benchmark:
        config["outdir"] + "/benchmarks/004_variant_call/02_haplotype_caller/{sample}_filter_indel.txt"
    params:
        gatk=config["gatk"],
        ref=config["reference_genome"],
        dcovg=config["gatk_params"]["HaplotypeCaller"]["dcovg"]
    shell:
        """
        java -jar \
        {params.gatk} \
        -T VariantFiltration \
        -R {params.ref} \
        --downsample_to_coverage {params.dcovg} \
        --downsampling_type BY_SAMPLE \
        --variant {input.indel_vcf} \
        -o {output.filtered_indel_vcf} \
        -baq CALCULATE_AS_NECESSARY \
        --filterExpression "QD < 2.0" \
        --filterExpression "ReadPosRankSum < -20.0" \
        --filterExpression "FS > 200.0" \
        --filterName QDFilter \
        --filterName ReadPosFilter \
        --filterName FSFilter \
        > {log.out_log} 2>&1
        """

rule haplotype_select_snp:
    message:
        "Selecting SNP variants for sample {wildcards.sample}"
    input:
        annotated_vcf=rules.haplotype_variant_annotator.output.annotated_vcf
    output:
        snp_vcf=config["outdir"] + "/analysis/004_variant_call/02_haplotype_caller/{sample}_onTarget.snp.vcf"
    log:
        out_log=config["outdir"] + "/logs/004_variant_call/02_haplotype_caller/{sample}_onTarget.select_snp.log"
    conda:
        "../envs/java.yml"
    benchmark:
        config["outdir"] + "/benchmarks/004_variant_call/02_haplotype_caller/{sample}_select_snp.txt"
    params:
        gatk=config["gatk"],
        ref=config["reference_genome"],
        dcovg=config["gatk_params"]["HaplotypeCaller"]["dcovg"],
        threads=config["threads_high"]
    shell:
        """
        java -jar \
        {params.gatk} \
        -T SelectVariants \
        -R {params.ref} \
        --downsample_to_coverage {params.dcovg} \
        --downsampling_type BY_SAMPLE \
        -nt {params.threads} \
        --variant {input.annotated_vcf} \
        --selectTypeToInclude SNP \
        -o {output.snp_vcf} \
        > {log.out_log} 2>&1
        """

rule haplotype_filter_snp:
    message:
        "Filtering SNP variants for sample {wildcards.sample}"
    input:
        snp_vcf=rules.haplotype_select_snp.output.snp_vcf,
        indel_vcf=rules.haplotype_select_indel.output.indel_vcf
    output:
        filtered_snp_vcf=config["outdir"] + "/analysis/004_variant_call/02_haplotype_caller/{sample}_onTarget.snp.filtered.vcf"
    log:
        out_log=config["outdir"] + "/logs/004_variant_call/02_haplotype_caller/{sample}_onTarget.filter_snp.log"
    conda:
        "../envs/java.yml"
    benchmark:
        config["outdir"] + "/benchmarks/004_variant_call/02_haplotype_caller/{sample}_filter_snp.txt"
    params:
        gatk=config["gatk"],
        ref=config["reference_genome"],
        dcovg=config["gatk_params"]["HaplotypeCaller"]["dcovg"]
    shell:
        """
        java -jar \
        {params.gatk} \
        -T VariantFiltration \
        -R {params.ref} \
        --downsample_to_coverage {params.dcovg} \
        --downsampling_type BY_SAMPLE \
        --variant {input.snp_vcf} \
        --mask {input.indel_vcf} \
        --maskName InDel \
        -o {output.filtered_snp_vcf} \
        --filterExpression "QD < 2.0" \
        --filterExpression "MQ < 40.0" \
        --filterExpression "FS > 60.0" \
        --filterExpression "MQRankSum < -12.5" \
        --filterExpression "ReadPosRankSum < -8.0" \
        --filterName QDFilter \
        --filterName MQFilter \
        --filterName FSFilter \
        --filterName MQRankSumFilter \
        --filterName ReadPosFilter \
        > {log.out_log} 2>&1
        """