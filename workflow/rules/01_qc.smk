rule raw_fastqc:
    message: 
        "Running Analysis on sample {lane}"
    singularity: 
        "https://depot.galaxyproject.org/singularity/fastqc%3A0.11.2--1"
    input:
        config["inputdir"] + "/{sample_filename}_{lane}_{R}_001.fastq.gz",
        "/mnt/omar/pipelines/icc/input"
    output:
        zip=config["outdir"] + "/analysis/001_QC/01_pretrim/{sample_filename}_{lane}_{R}_pretrim_fastqc.zip"
    threads: 
        config["threads_low"]
    params: 
        path=config["outdir"] + "/analysis/001_QC/01_pretrim/{sample_filename}",
    log:
        config["outdir"] + "/logs/001_QC/01_pretrim/{sample_filename}_{lane}_{R}.log"
    benchmark:
        config["outdir"] + "/benchmarks/001_QC/01_pretrim/{sample_filename}_{lane}_{R}.txt"
    shell:
        """
        # Generate parent directory path
        parent_path=$(dirname {params.path})

        # Create the output directory
        mkdir -p "$parent_path"

        fastqc {input[0]} \
        -t {threads} \
        -o "$parent_path" \
        --noextract \
        &> {log}

        # Move the FastQC outputs to the desired location
        
        mv {params.path}_{wildcards.lane}_{wildcards.R}_001_fastqc.zip {output.zip}
        """