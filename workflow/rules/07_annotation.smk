rule variant_effect_predictor:
    message:
        "Running VEP annotation for {wildcards.sample} {wildcards.caller} {wildcards.variant_type}"
    input:
        vcf=config["outdir"] + "/analysis/004_variant_call/*{caller}/{sample}_onTarget.{variant_type}.filtered.vcf"
    output:
        vep_vcf=config["outdir"] + "/analysis/005_annotation/01_vep/{sample}_{target}.{caller}.{variant_type}.vep.vcf"
    conda:
        "../envs/vep.yml"
    log:
        config["outdir"] + "/logs/005_annotation/01_vep/{sample}_{target}.{caller}.{variant_type}.vep.log"
    benchmark:
        config["outdir"] + "/benchmarks/005_annotation/01_vep/{sample}_{target}.{caller}.{variant_type}.vep.txt"
    threads:
        config["threads_high"]
    params:
        vep_path=config["vep_path"],
        vep_dir=config["vep_dir"],
        vep_fasta=config["vep_fasta"]
    shell:
        """
        perl {params.vep_path}/variant_effect_predictor.pl \
        --dir_cache {params.vep_dir} \
        --assembly GRCh37 \
        --fasta {params.vep_fasta} \
        --cache --offline --no_progress --fork {threads} \
        -i {input.vcf} -o {output.vep_vcf} \
        --allele_number --hgvs --check_existing --canonical --ccds --maf_exac --pubmed --protein --sift b --polyphen b \
        --fields ALLELE_NUM,Consequence,SYMBOL,SYMBOL_SOURCE,Gene,Feature,BIOTYPE,HGVSc,HGVSp,cDNA_position,CDS_position,Protein_position,Amino_acids,Codons,Existing_variation,STRAND,CANONICAL,CCDS,ENSP,SIFT,PolyPhen,ExAC_MAF,PUBMED \
        --vcf --force_overwrite \
        > {log} 2>&1
        """

rule tableize_vep_vcf:
    message:
        "Converting VEP-annotated VCF to tabular format for {wildcards.sample} {wildcards.caller} {wildcards.variant_type}"
    input:
        vep_vcf=rules.variant_effect_predictor.output.vep_vcf
    output:
        vep_table=config["outdir"] + "/analysis/005_annotation/01_vep/{sample}_{target}.{caller}.{variant_type}.tableize.txt"
    conda:
        "../envs/python.yml"
    log:
        config["outdir"] + "/logs/005_annotation/01_vep/{sample}_{target}.{caller}.{variant_type}.tableize.log"
    benchmark:
        config["outdir"] + "/benchmarks/005_annotation/01_vep/{sample}_{target}.{caller}.{variant_type}.tableize.txt"
    params:
        loftee_path=config["loftee_path"]
    shell:
        """
        python {params.loftee_path}/tableize_vcf.py \
        --vcf {input.vep_vcf} \
        --out {output.vep_table} \
        --split_by_transcript \
        --all_csqs \
        --do_not_minrep \
        --include_id \
        --info ABHet,ABHom,AC,AF,AN,DP,FS,HaplotypeScore,MLEAC,MLEAF,MQ,MQ0,QD \
        --vep_info Consequence,SYMBOL,SYMBOL_SOURCE,Gene,Feature,BIOTYPE,HGVSc,HGVSp,cDNA_position,CDS_position,Protein_position,Amino_acids,Codons,Existing_variation,STRAND,CANONICAL,CCDS,ENSP,SIFT,PolyPhen,ExAC_MAF,PUBMED \
        --mysql \
        > {log} 2>&1
        """