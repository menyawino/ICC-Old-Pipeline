rule qc_flagstat_target:
    input:
        bam=rules.on_target_map_qual8.output.out_bam
    output:
        flagstat=config["outdir"] + "/analysis/003_alignment/02_bam_qc/{sample}_target_qc_flagstat.txt"
    params:
        samtools=config["samtools"]
    shell:
        """
        {params.samtools} \
        view -u \
        -F 0x100 {input.bam} \
        | {params.samtools} \
        flagstat - > {output.flagstat}
        """

rule qc_flagstat_protein:
    input:
        bam=rules.on_prot_coding_target_map_qual8.output.out_bam
    output:
        flagstat=config["outdir"] + "/analysis/003_alignment/02_bam_qc/{sample}_protein_qc_flagstat.txt"
    params:
        samtools=config["samtools"]
    shell:
        """
        {params.samtools} \
        view -u \
        -F 0x100 {input.bam} \
        | {params.samtools} \
        flagstat - > {output.flagstat}
        """

rule qc_flagstat_canonical:
    input:
        bam=rules.on_canonical_target_map_qual8.output.out_bam
    output:
        flagstat=config["outdir"] + "/analysis/003_alignment/02_bam_qc/{sample}_canonical_qc_flagstat.txt"
    params:
        samtools=config["samtools"]
    shell:
        """
        {params.samtools} \
        view -u \
        -F 0x100 {input.bam} \
        | {params.samtools} \
        flagstat - > {output.flagstat}
        """

rule qc_flagstat_original:
    input:
        bam=rules.indel_realigner_merged.output.bam_realigned
    output:
        flagstat=config["outdir"] + "/analysis/003_alignment/02_bam_qc/{sample}_original_qc_flagstat.txt"
    params:
        samtools=config["samtools"]
    shell:
        """
        {params.samtools} \
        view -u \
        -F 0x100 {input.bam} \
        | {params.samtools} \
        flagstat - > {output.flagstat}
        """

rule qc_coverage_stats_protein:
    input:
        bam=rules.on_prot_coding_target_map_qual8.output.out_bam
    output:
        cov_stats=config["outdir"] + "/analysis/003_alignment/02_bam_qc/{sample}_protein_qc_coverage_stats.txt"
    log:
        config["outdir"] + "/analysis/003_alignment/02_bam_qc/{sample}_protein_qc_coverage_stats.log"
    benchmark:
        config["outdir"] + "/analysis/003_alignment/02_bam_qc/{sample}_protein_qc_coverage_stats.benchmark"
    params:
        samtools=config["samtools"],
        bedtools=config["bedtools"],
        cds_file=config["cds_panel"]
    shell:
        """
        {params.samtools} \
        view -uF 0x400 \
        {input.bam} \
        | {params.bedtools}/coverageBed \
        -abam stdin \
        -b {params.cds_file} \
        > {output.cov_stats} \
        2> {log}
        """

rule qc_coverage_per_base_protein:
    input:
        bam=rules.on_prot_coding_target_map_qual8.output.out_bam
    output:
        per_base=config["outdir"] + "/analysis/003_alignment/02_bam_qc/{sample}_protein_qc_coverage_per_base.txt"
    log:
        config["outdir"] + "/analysis/003_alignment/02_bam_qc/{sample}_protein_qc_coverage_per_base.log"
    benchmark:
        config["outdir"] + "/analysis/003_alignment/02_bam_qc/{sample}_protein_qc_coverage_per_base.benchmark"
    params:
        samtools=config["samtools"],
        bedtools=config["bedtools"],
        cds_file=config["cds_panel"]
    shell:
        """
        {params.samtools} \
        view -uF 0x400 \
        {input.bam} \
        | {params.bedtools}/coverageBed \
        -abam stdin \
        -b {params.cds_file} \
        -d \
        > {output.per_base} \
        2> {log}
        """

rule qc_coverage_histogram_protein:
    input:
        bam=rules.on_prot_coding_target_map_qual8.output.out_bam
    output:
        hist=config["outdir"] + "/analysis/003_alignment/02_bam_qc/{sample}_protein_qc_coverage_histogram.txt"
    log:
        config["outdir"] + "/analysis/003_alignment/02_bam_qc/{sample}_protein_qc_coverage_histogram.log"
    benchmark:
        config["outdir"] + "/analysis/003_alignment/02_bam_qc/{sample}_protein_qc_coverage_histogram.benchmark"
    params:
        samtools=config["samtools"],
        bedtools=config["bedtools"],
        cds_file=config["cds_panel"]
    shell:
        """
        {params.samtools} \
        view -uF 0x400 \
         {input.bam} \
        | {params.bedtools}/coverageBed \
        -abam stdin \
        -b {params.cds_file} \
        -hist \
        > {output.hist} \
        2> {log}
        """

rule qc_sort_coverage_stats:
    input:
        cov_stats=rules.qc_coverage_stats_protein.output.cov_stats
    output:
        cov_stats_sorted=config["outdir"] + "/analysis/003_alignment/02_bam_qc/{sample}_protein_qc_coverage_stats_sorted.txt"
    shell:
        """
        sort -k1,1 -k2,2n {input.cov_stats} > {output.cov_stats_sorted}
        """

rule qc_sort_coverage_per_base:
    input:
        per_base=rules.qc_coverage_per_base_protein.output.per_base
    output:
        per_base_sorted=config["outdir"] + "/analysis/003_alignment/02_bam_qc/{sample}_protein_qc_coverage_per_base_sorted.txt"
    shell:
        """
        sort -k1,1 -k2,2n {input.per_base} > {output.per_base_sorted}
        """

rule qc_callable_loci_protein_min10:
    input:
        bam=rules.on_prot_coding_target_map_qual8.output.out_bam
    output:
        callable=config["outdir"] + "/analysis/003_alignment/02_bam_qc/{sample}_protein_min10_callable.bed"
    conda:
        "../envs/java.yml"
    log:
        config["outdir"] + "/logs/003_alignment/02_bam_qc/{sample}_protein_min10_callable.log"
    params:
        gatk=config["gatk"],
        ref=config["reference_genome"],
        cds_file=config["cds_panel"]
    shell:
        """
        java -jar \
        {params.gatk} \
        -T CallableLoci \
        -I {input.bam} \
        -R {params.ref} \
        --minDepth 10 \
        -L {params.cds_file} \
        -summary {output.callable}.summary \
        -o {output.callable} \
        2> {log}
        """

rule qc_callable_loci_protein_min20:
    input:
        bam=rules.on_prot_coding_target_map_qual8.output.out_bam
    output:
        callable=config["outdir"] + "/analysis/003_alignment/02_bam_qc/{sample}_protein_min20_callable.bed"
    conda:
        "../envs/java.yml"
    log:
        config["outdir"] + "/logs/003_alignment/02_bam_qc/{sample}_protein_min20_callable.log"
    params:
        gatk=config["gatk"],
        ref=config["reference_genome"],
        cds_file=config["cds_panel"]
    shell:
        """
        java -jar \
        {params.gatk} \
        -T CallableLoci \
        -I {input.bam} \
        -R {params.ref} \
        --minDepth 20 \
        -L {params.cds_file} \
        -summary {output.callable}.summary \
        -o {output.callable} \
        2> {log}
        """

rule qc_callable_loci_protein_min30:
    input:
        bam=rules.on_prot_coding_target_map_qual8.output.out_bam
    output:
        callable=config["outdir"] + "/analysis/003_alignment/02_bam_qc/{sample}_protein_min30_callable.bed"
    conda:
        "../envs/java.yml"
    log:
        config["outdir"] + "/logs/003_alignment/02_bam_qc/{sample}_protein_min30_callable.log"
    params:
        gatk=config["gatk"],
        ref=config["reference_genome"],
        cds_file=config["cds_panel"]
    shell:
        """
        java -jar \
        {params.gatk} \
        -T CallableLoci \
        -I {input.bam} \
        -R {params.ref} \
        --minDepth 30 \
        -L {params.cds_file} \
        -summary {output.callable}.summary \
        -o {output.callable} \
        2> {log}
        """

rule qc_depth_of_coverage_protein:
    input:
        bam=rules.on_prot_coding_target_map_qual8.output.out_bam
    output:
        depth=config["outdir"] + "/analysis/003_alignment/02_bam_qc/{sample}_protein_depth_of_coverage.txt"
    conda:
        "../envs/java.yml"
    log:
        config["outdir"] + "/logs/003_alignment/02_bam_qc/{sample}_protein_depth_of_coverage.log"
    params:
        gatk=config["gatk"],
        ref=config["reference_genome"],
        cds_file=config["cds_panel"]
    shell:
        """
        java -jar \
        {params.gatk} \
        -T DepthOfCoverage \
        -I {input.bam} \
        -R {params.ref} \
        -L {params.cds_file} \
        -o {output.depth} \
        2> {log}
        """

rule qc_mean_coverage_exon_protein:
    input:
        depth_summary=config["outdir"] + "/analysis/003_alignment/02_bam_qc/{sample}_protein_depth_of_coverage.txt",
        cov_stats=config["outdir"] + "/analysis/003_alignment/02_bam_qc/{sample}_protein_qc_coverage_stats.txt"
    output:
        mean_cvg_bed=config["outdir"] + "/analysis/003_alignment/02_bam_qc/{sample}_protein_mean_coverage_exon.bed"
    params:
        bedtools=config["bedtools"]
    shell:
        """
        cat {input.depth_summary} \
        | awk '{{print $1, $2, $3, $4}}' \
        | paste - {input.cov_stats} > {output.mean_cvg_bed}
        """

rule qc_collect_alignment_summary_metrics_protein:
    input:
        bam=rules.on_prot_coding_target_map_qual8.output.out_bam
    output:
        align_sum=config["outdir"] + "/analysis/003_alignment/02_bam_qc/{sample}_protein_alignment_summary_metrics.txt"
    conda:
        "../envs/java.yml"
    log:
        config["outdir"] + "/logs/003_alignment/02_bam_qc/{sample}_protein_alignment_metrics.log"
    params:
        picard=config["picard"],
        ref=config["reference_genome"]
    shell:
        """
        java -jar \
        {params.picard}/CollectAlignmentSummaryMetrics.jar \
        I={input.bam} \
        R={params.ref} \
        O={output.align_sum} \
        2> {log}
        """

rule qc_callable_loci_canonical_min10:
    input:
        bam=rules.on_canonical_target_map_qual8.output.out_bam
    output:
        callable=config["outdir"] + "/analysis/003_alignment/02_bam_qc/{sample}_canonical_min10_callable.bed"
    conda:
        "../envs/java.yml"
    log:
        config["outdir"] + "/logs/003_alignment/02_bam_qc/{sample}_canonical_min10_callable.log"
    params:
        gatk=config["gatk"],
        ref=config["reference_genome"],
        cds_file=config["cds_panel"]
    shell:
        """
        java -jar \
        {params.gatk} \
        -T CallableLoci \
        -I {input.bam} \
        -R {params.ref} \
        --minDepth 10 \
        -L {params.cds_file} \
        -summary {output.callable}.summary \
        -o {output.callable} \
        2> {log}
        """

rule qc_depth_of_coverage_canonical:
    input:
        bam=rules.on_canonical_target_map_qual8.output.out_bam
    output:
        depth=config["outdir"] + "/analysis/003_alignment/02_bam_qc/{sample}_canonical_depth_of_coverage.txt"
    conda:
        "../envs/java.yml"
    log:
        config["outdir"] + "/logs/003_alignment/02_bam_qc/{sample}_canonical_depth_of_coverage.log"
    params:
        gatk=config["gatk"],
        ref=config["reference_genome"],
        cds_file=config["cds_panel"]
    shell:
        """
        java -jar \
        {params.gatk} \
        -T DepthOfCoverage \
        -I {input.bam} \
        -R {params.ref} \
        -L {params.cds_file} \
        -o {output.depth} \
        2> {log}
        """

rule qc_mean_coverage_exon_canonical:
    input:
        depth_summary=config["outdir"] + "/analysis/003_alignment/02_bam_qc/{sample}_canonical_depth_of_coverage.txt",
        cov_stats=config["outdir"] + "/analysis/003_alignment/02_bam_qc/{sample}_canonical_qc_coverage_stats.txt"
    output:
        mean_cvg_bed=config["outdir"] + "/analysis/003_alignment/02_bam_qc/{sample}_canonical_mean_coverage_exon.bed"
    params:
        bedtools=config["bedtools"]
    shell:
        """
        cat {input.depth_summary} | awk '{{print $1, $2, $3, $4}}' | paste - {input.cov_stats} > {output.mean_cvg_bed}
        """

rule qc_collect_alignment_summary_metrics_canonical:
    input:
        bam=rules.on_canonical_target_map_qual8.output.out_bam
    output:
        align_sum=config["outdir"] + "/analysis/003_alignment/02_bam_qc/{sample}_canonical_alignment_summary_metrics.txt"
    conda:
        "../envs/java.yml"
    log:
        config["outdir"] + "/logs/003_alignment/02_bam_qc/{sample}_canonical_alignment_metrics.log"
    params:
        picard=config["picard"],
        ref=config["reference_genome"]
    shell:
        """
        java -jar \
        {params.picard}/CollectAlignmentSummaryMetrics.jar \
        I={input.bam} \
        R={params.ref} \
        O={output.align_sum} \
        2> {log}
        """

rule qc_callable_loci_target_min10:
    input:
        bam=rules.on_target_map_qual8.output.out_bam
    output:
        callable=config["outdir"] + "/analysis/003_alignment/02_bam_qc/{sample}_target_min10_callable.bed"
    conda:
        "../envs/java.yml"
    log:
        config["outdir"] + "/logs/003_alignment/02_bam_qc/{sample}_target_min10_callable.log"
    params:
        gatk=config["gatk"],
        ref=config["reference_genome"],
        cds_file=config["cds_panel"]
    shell:
        """
        java -jar \
        {params.gatk} \
        -T CallableLoci \
        -I {input.bam} \
        -R {params.ref} \
        --minDepth 10 \
        -L {params.cds_file} \
        -summary {output.callable}.summary \
        -o {output.callable} \
        2> {log}
        """

rule qc_depth_of_coverage_target:
    input:
        bam=rules.on_target_map_qual8.output.out_bam
    output:
        depth=config["outdir"] + "/analysis/003_alignment/02_bam_qc/{sample}_target_depth_of_coverage.txt"
    conda:
        "../envs/java.yml"
    log:
        config["outdir"] + "/logs/003_alignment/02_bam_qc/{sample}_target_depth_of_coverage.log"
    params:
        gatk=config["gatk"],
        ref=config["reference_genome"],
        cds_file=config["cds_panel"]
    shell:
        """
        java -jar \
        {params.gatk} \
        -T DepthOfCoverage \
        -I {input.bam} \
        -R {params.ref} \
        -L {params.cds_file} \
        -o {output.depth} \
        2> {log}
        """

rule qc_mean_coverage_exon_target:
    message:
        "Calculating mean coverage for target regions in sample {wildcards.sample}"
    input:
        depth_summary=config["outdir"] + "/analysis/003_alignment/02_bam_qc/{sample}_target_depth_of_coverage.txt",
        cov_stats=config["outdir"] + "/analysis/003_alignment/02_bam_qc/{sample}_target_qc_coverage_stats.txt"
    output:
        mean_cvg_bed=config["outdir"] + "/analysis/003_alignment/02_bam_qc/{sample}_target_mean_coverage_exon.bed"
    log:
        config["outdir"] + "/analysis/003_alignment/02_bam_qc/{sample}_target_mean_coverage_exon.log"
    benchmark:
        config["outdir"] + "/analysis/003_alignment/02_bam_qc/{sample}_target_mean_coverage_exon.benchmark"
    shell:
        """
        cat {input.depth_summary} \
        | awk '{{print $1, $2, $3, $4}}' \
        | paste - {input.cov_stats} > {output.mean_cvg_bed} \
        2> {log}
        """

rule qc_collect_alignment_summary_metrics_target:
    input:
        bam=rules.on_target_map_qual8.output.out_bam
    output:
        align_sum=config["outdir"] + "/analysis/003_alignment/02_bam_qc/{sample}_target_alignment_summary_metrics.txt"
    conda:
        "../envs/java.yml"
    log:
        config["outdir"] + "/analysis/003_alignment/02_bam_qc/{sample}_target_alignment_summary_metrics.log"
    params:
        picard=config["picard"],
        ref=config["reference_genome"]
    shell:
        """
        java -jar \
        {params.picard}/CollectAlignmentSummaryMetrics.jar \
        I={input.bam} \
        R={params.ref} \
        O={output.align_sum} \
        2> {log}
        """

rule qc_coverage_stats_canonical:
    input:
        bam=rules.on_canonical_target_map_qual8.output.out_bam
    output:
        cov_stats=config["outdir"] + "/analysis/003_alignment/02_bam_qc/{sample}_canonical_qc_coverage_stats.txt"
    log:
        config["outdir"] + "/analysis/003_alignment/02_bam_qc/{sample}_canonical_qc_coverage_stats.log"
    params:
        samtools=config["samtools"],
        bedtools=config["bedtools"],
        canon_file=config["canontran_panel"]
    shell:
        """
        {params.samtools} \
        view -uF 0x400 \
        {input.bam} \
        | {params.bedtools}/coverageBed \
        -abam stdin \
        -b {params.canon_file} \
        > {output.cov_stats} \
        2> {log}
        """

rule qc_coverage_per_base_canonical:
    input:
        bam=rules.on_canonical_target_map_qual8.output.out_bam
    output:
        per_base=config["outdir"] + "/analysis/003_alignment/02_bam_qc/{sample}_canonical_qc_coverage_per_base.txt"
    log:
        config["outdir"] + "/analysis/003_alignment/02_bam_qc/{sample}_canonical_qc_coverage_per_base.log"
    params:
        samtools=config["samtools"],
        bedtools=config["bedtools"],
        canon_file=config["canontran_panel"]
    shell:
        """
        {params.samtools} \
        view -uF 0x400 \
        {input.bam} \
        | {params.bedtools}/coverageBed \
        -abam stdin \
        -b {params.canon_file} \
        -d \
        > {output.per_base} \
        2> {log}
        """

rule qc_coverage_histogram_canonical:
    input:
        bam=rules.on_canonical_target_map_qual8.output.out_bam
    output:
        hist=config["outdir"] + "/analysis/003_alignment/02_bam_qc/{sample}_canonical_qc_coverage_histogram.txt"
    log:
        config["outdir"] + "/analysis/003_alignment/02_bam_qc/{sample}_canonical_qc_coverage_histogram.log"
    params:
        samtools=config["samtools"],
        bedtools=config["bedtools"],
        canon_file=config["canontran_panel"]
    shell:
        """
        {params.samtools} \
        view -uF 0x400 \
        {input.bam} \
        | {params.bedtools}/coverageBed \
        -abam stdin \
        -b {params.canon_file} \
        -hist \
        > {output.hist} \
        2> {log}
        """

rule qc_coverage_stats_target:
    input:
        bam=rules.on_target_map_qual8.output.out_bam
    output:
        cov_stats=config["outdir"] + "/analysis/003_alignment/02_bam_qc/{sample}_target_qc_coverage_stats.txt"
    log:
        config["outdir"] + "/analysis/003_alignment/02_bam_qc/{sample}_target_qc_coverage_stats.log"
    params:
        samtools=config["samtools"],
        bedtools=config["bedtools"],
        target_file=config["icc_panel"]
    shell:
        """
        {params.samtools} \
        view -uF 0x400 \
        {input.bam} \
        | {params.bedtools}/coverageBed \
        -abam stdin \
        -b {params.target_file} \
        > {output.cov_stats} \
        2> {log}
        """

rule qc_coverage_per_base_target:
    input:
        bam=rules.on_target_map_qual8.output.out_bam
    output:
        per_base=config["outdir"] + "/analysis/003_alignment/02_bam_qc/{sample}_target_qc_coverage_per_base.txt"
    log:
        config["outdir"] + "/analysis/003_alignment/02_bam_qc/{sample}_target_qc_coverage_per_base.log"
    params:
        samtools=config["samtools"],
        bedtools=config["bedtools"],
        target_file=config["icc_panel"]
    shell:
        """
        {params.samtools} \
        view -uF 0x400 \
        {input.bam} \
        | {params.bedtools}/coverageBed \
        -abam stdin \
        -b {params.target_file} \
        -d \
        > {output.per_base} \
        2> {log}
        """

rule qc_coverage_histogram_target:
    input:
        bam=rules.on_target_map_qual8.output.out_bam
    output:
        hist=config["outdir"] + "/analysis/003_alignment/02_bam_qc/{sample}_target_qc_coverage_histogram.txt"
    conda:
        "../envs/java.yml"
    log:
        config["outdir"] + "/logs/003_alignment/02_bam_qc/{sample}_target_qc_coverage_histogram.log"
    params:
        samtools=config["samtools"],
        bedtools=config["bedtools"],
        target_file=config["icc_panel"]
    shell:
        """
        {params.samtools} \
        view -uF 0x400 \
        {input.bam} \
        | {params.bedtools}/coverageBed \
        -abam stdin \
        -b {params.target_file} \
        -hist \
        > {output.hist} \
        2> {log}
        """

rule qc_callable_loci_canonical_min20:
    input:
        bam=rules.on_canonical_target_map_qual8.output.out_bam
    output:
        callable=config["outdir"] + "/analysis/003_alignment/02_bam_qc/{sample}_canonical_min20_callable.bed"
    conda:
        "../envs/java.yml"
    log:
        config["outdir"] + "/logs/003_alignment/02_bam_qc/{sample}_canonical_min20_callable.log"
    params:
        gatk=config["gatk"],
        ref=config["reference_genome"],
        canon_file=config["canontran_panel"]
    shell:
        """
        java -jar \
        {params.gatk} \
        -T CallableLoci \
        -I {input.bam} \
        -R {params.ref} \
        --minDepth 20 \
        -L {params.canon_file} \
        -summary {output.callable}.summary \
        -o {output.callable} \
        2> {log}
        """

rule qc_callable_loci_canonical_min30:
    input:
        bam=rules.on_canonical_target_map_qual8.output.out_bam
    output:
        callable=config["outdir"] + "/analysis/003_alignment/02_bam_qc/{sample}_canonical_min30_callable.bed"
    conda:
        "../envs/java.yml"
    log:
        config["outdir"] + "/logs/003_alignment/02_bam_qc/{sample}_canonical_min30_callable.log"
    params:
        gatk=config["gatk"],
        ref=config["reference_genome"],
        canon_file=config["canontran_panel"]
    shell:
        """
        java -jar \
        {params.gatk} \
        -T CallableLoci \
        -I {input.bam} \
        -R {params.ref} \
        --minDepth 30 \
        -L {params.canon_file} \
        -summary {output.callable}.summary \
        -o {output.callable} \
        2> {log}
        """

rule qc_callable_loci_target_min20:
    input:
        bam=rules.on_target_map_qual8.output.out_bam
    output:
        callable=config["outdir"] + "/analysis/003_alignment/02_bam_qc/{sample}_target_min20_callable.bed"
    conda:
        "../envs/java.yml"
    params:
        gatk=config["gatk"],
        ref=config["reference_genome"],
        target_file=config["icc_panel"]
    log:
        config["outdir"] + "/logs/003_alignment/02_bam_qc/{sample}_target_min20_callable.log"
    benchmark:
        config["outdir"] + "/benchmarks/003_alignment/02_bam_qc/{sample}_target_min20_callable.benchmark"
    shell:
        """
        java -jar \
        {params.gatk} \
        -T CallableLoci \
        -I {input.bam} \
        -R {params.ref} \
        --minDepth 20 \
        -L {params.target_file} \
        -summary {output.callable}.summary \
        -o {output.callable} \
        2> {log}
        """

rule qc_callable_loci_target_min30:
    input:
        bam=rules.on_target_map_qual8.output.out_bam
    output:
        callable=config["outdir"] + "/analysis/003_alignment/02_bam_qc/{sample}_target_min30_callable.bed"
    conda:
        "../envs/java.yml"
    params:
        gatk=config["gatk"],
        ref=config["reference_genome"],
        target_file=config["icc_panel"]
    log:
        config["outdir"] + "/logs/003_alignment/02_bam_qc/{sample}_target_min30_callable.log"
    benchmark:
        config["outdir"] + "/benchmarks/003_alignment/02_bam_qc/{sample}_target_min30_callable.benchmark"
    shell:
        """
        java -jar \
        {params.gatk} \
        -T CallableLoci \
        -I {input.bam} \
        -R {params.ref} \
        --minDepth 30 \
        -L {params.target_file} \
        -summary {output.callable}.summary \
        -o {output.callable} \
        2> {log}
        """

rule collect_all_qc:
    input:
        flagstat=expand(config["outdir"] + "/analysis/003_alignment/02_bam_qc/{sample}_{type}_qc_flagstat.txt", sample="{sample}", type=["target", "protein", "canonical", "original"]),
        hist=rules.qc_coverage_histogram_protein.output.hist,
        sorted_cov_stats=rules.qc_sort_coverage_stats.output.cov_stats_sorted,
        sorted_per_base=rules.qc_sort_coverage_per_base.output.per_base_sorted,
        qc_callable_loci_protein_min10=rules.qc_callable_loci_protein_min10.output.callable,
        qc_callable_loci_protein_min20=rules.qc_callable_loci_protein_min20.output.callable,
        qc_callable_loci_protein_min30=rules.qc_callable_loci_protein_min30.output.callable,
        qc_callable_loci_canonical_min10=rules.qc_callable_loci_canonical_min10.output.callable,
        qc_callable_loci_canonical_min20=rules.qc_callable_loci_canonical_min20.output.callable,
        qc_callable_loci_canonical_min30=rules.qc_callable_loci_canonical_min30.output.callable,
        qc_callable_loci_target_min10=rules.qc_callable_loci_target_min10.output.callable,
        qc_callable_loci_target_min20=rules.qc_callable_loci_target_min20.output.callable,
        qc_callable_loci_target_min30=rules.qc_callable_loci_target_min30.output.callable,
        qc_depth_of_coverage_protein=rules.qc_depth_of_coverage_protein.output.depth,
        qc_depth_of_coverage_canonical=rules.qc_depth_of_coverage_canonical.output.depth,
        qc_depth_of_coverage_target=rules.qc_depth_of_coverage_target.output.depth,
        qc_mean_coverage_exon_protein=rules.qc_mean_coverage_exon_protein.output.mean_cvg_bed,
        qc_mean_coverage_exon_canonical=rules.qc_mean_coverage_exon_canonical.output.mean_cvg_bed,
        qc_mean_coverage_exon_target=rules.qc_mean_coverage_exon_target.output.mean_cvg_bed,
        qc_collect_alignment_summary_metrics_protein=rules.qc_collect_alignment_summary_metrics_protein.output.align_sum,
        qc_collect_alignment_summary_metrics_canonical=rules.qc_collect_alignment_summary_metrics_canonical.output.align_sum,
        qc_collect_alignment_summary_metrics_target=rules.qc_collect_alignment_summary_metrics_target.output.align_sum,
        qc_coverage_stats_protein=rules.qc_coverage_stats_protein.output.cov_stats,
        qc_coverage_stats_canonical=rules.qc_coverage_stats_canonical.output.cov_stats,
        qc_coverage_stats_target=rules.qc_coverage_stats_target.output.cov_stats,
        qc_coverage_per_base_protein=rules.qc_coverage_per_base_protein.output.per_base,
        qc_coverage_per_base_canonical=rules.qc_coverage_per_base_canonical.output.per_base,
        qc_coverage_per_base_target=rules.qc_coverage_per_base_target.output.per_base,
        qc_coverage_histogram_protein=rules.qc_coverage_histogram_protein.output.hist,
        qc_coverage_histogram_canonical=rules.qc_coverage_histogram_canonical.output.hist,
        qc_coverage_histogram_target=rules.qc_coverage_histogram_target.output.hist,
    output:
        all_qc=config["outdir"] + "/analysis/003_alignment/02_bam_qc/{sample}_all_qc.txt"
    shell:
        """
        cat {input.flagstat} \
        {input.hist} \
        {input.sorted_cov_stats} \
        {input.sorted_per_base} \
        > {output.all_qc}
        """