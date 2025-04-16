import glob
import os

# Trim reads using prinseq
rule trimming:
    singularity:
        "docker://flowcraft/prinseq:0.20.4-1"
    message: 
        "Trimming and removing adapters from sample {wildcards.sample}_{lane}"
    conda: 
        "icc_02_trimming"
    input:
        fq1=lambda wildcards: [f for f in glob.glob(config["inputdir"] + f"/{wildcards.sample.split('/')[0]}/{wildcards.sample.split('/')[-1].split('_')[0]}_S*_{wildcards.lane}_R1_001.fastq.gz") or 
                               glob.glob(config["inputdir"] + f"/{wildcards.sample.split('/')[0]}/{wildcards.sample.split('/')[-1].split('_')[0]}_S*_{wildcards.lane}_R1_001.fastq") or
                               glob.glob(config["inputdir"] + f"/{wildcards.sample.split('/')[-1].split('_')[0]}_S*_{wildcards.lane}_R1_001.fastq.gz") or
                               glob.glob(config["inputdir"] + f"/{wildcards.sample.split('/')[-1].split('_')[0]}_S*_{wildcards.lane}_R1_001.fastq")][0],
        fq2=lambda wildcards: [f for f in glob.glob(config["inputdir"] + f"/{wildcards.sample.split('/')[0]}/{wildcards.sample.split('/')[-1].split('_')[0]}_S*_{wildcards.lane}_R2_001.fastq.gz") or 
                               glob.glob(config["inputdir"] + f"/{wildcards.sample.split('/')[0]}/{wildcards.sample.split('/')[-1].split('_')[0]}_S*_{wildcards.lane}_R2_001.fastq") or
                               glob.glob(config["inputdir"] + f"/{wildcards.sample.split('/')[-1].split('_')[0]}_S*_{wildcards.lane}_R2_001.fastq.gz") or
                               glob.glob(config["inputdir"] + f"/{wildcards.sample.split('/')[-1].split('_')[0]}_S*_{wildcards.lane}_R2_001.fastq")][0]
    output:
        fq1=config["outdir"] + "/analysis/002_trimming/{sample}_{lane}_R1.fastq",
        fq1s=config["outdir"] + "/analysis/002_trimming/{sample}_{lane}_R1_singletons.fastq",
        fq2=config["outdir"] + "/analysis/002_trimming/{sample}_{lane}_R2.fastq",
        fq2s=config["outdir"] + "/analysis/002_trimming/{sample}_{lane}_R2_singletons.fastq"
    params:
        path=config["outdir"] + "/analysis/002_trimming/{sample}_{lane}",
        sample_name=lambda wildcards: wildcards.sample.split("/")[-1]
    log:
        config["outdir"] + "/logs/002_trimming/{sample}/{sample}_{lane}.log"
    benchmark:
        config["outdir"] + "/benchmarks/002_trimming/{sample}/{sample}_{lane}.txt"
    shell:
        """
        # Create temp directory for intermediate files
        mkdir -p $(dirname {params.path})/temp

        # Check if input files are gzipped and handle accordingly
        if [[ "{input.fq1}" == *.gz ]]; then
            gunzip -c {input.fq1} > $(dirname {params.path})/temp/{params.sample_name}_{wildcards.lane}_R1.fastq
            FQ1=$(dirname {params.path})/temp/{params.sample_name}_{wildcards.lane}_R1.fastq
        else
            FQ1={input.fq1}
        fi

        if [[ "{input.fq2}" == *.gz ]]; then
            gunzip -c {input.fq2} > $(dirname {params.path})/temp/{params.sample_name}_{wildcards.lane}_R2.fastq
            FQ2=$(dirname {params.path})/temp/{params.sample_name}_{wildcards.lane}_R2.fastq
        else
            FQ2={input.fq2}
        fi

        prinseq-lite.pl \\
        -fastq $FQ1 \\
        -fastq2 $FQ2 \\
        -out_good {params.path} \\
        -out_bad null \\
        -trim_qual_right 20 \\
        -trim_qual_left 20 \\
        -trim_qual_window 5 \\
        -min_len 35 \\
        &> {log}

        # rename output files to suit the rest of the pipeline
        mv {params.path}_1.fastq {output.fq1}
        mv {params.path}_1_singletons.fastq {output.fq1s}
        mv {params.path}_2.fastq {output.fq2}
        mv {params.path}_2_singletons.fastq {output.fq2s}

        # Clean up temp directory
        rm -rf $(dirname {params.path})/temp
        """

