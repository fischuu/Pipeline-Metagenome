rule Concatenate_lanes:
    """
    Concatenate the demultiplexed fastq files (BASH).
    """
    input:
        R1=get_fastq_for_concatenating_read1,
        R2=get_fastq_for_concatenating_read2
    output:
        R1=temp("%s/FASTQ/CONCATENATED/{samples}_R1.fastq.gz" % (config["project-folder"])),
        R1Report="%s/FASTQ/CONCATENATED/{samples}_R1.fastq.gz.report" % (config["project-folder"]),
        R2=temp("%s/FASTQ/CONCATENATED/{samples}_R2.fastq.gz" % (config["project-folder"])),
        R2Report="%s/FASTQ/CONCATENATED/{samples}_R2.fastq.gz.report" % (config["project-folder"])
    log:
        "%s/logs/Concatenate_lanes.{samples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/Concatenate_lanes.{samples}.tsv" % (config["project-folder"])
    params:
        outfolder="%s/FASTQ/CONCATENATED" % (config["project-folder"])
    resources:
        time=cluster["Concatenate_lanes"]["time"],
        mem=cluster["Concatenate_lanes"]["mem-per-cpu"]
    threads: cluster["Concatenate_lanes"]["cpus-per-task"]
    shell:"""
        mkdir -p {params.outfolder}
        cat {input.R1} > {output.R1} 2> {log}
        ls {input.R1} > {output.R1Report} 2> {log}
        cat {input.R2} > {output.R2} 2> {log}
        ls {input.R2} > {output.R2Report} 2> {log}
  	"""
  	
rule fastp_trim_reads:
    """
    Trim adapters and low quality reads (FASTP).
    """
    input:
        R1="%s/FASTQ/CONCATENATED/{samples}_R1.fastq.gz" % (config["project-folder"]),
        R2="%s/FASTQ/CONCATENATED/{samples}_R2.fastq.gz" % (config["project-folder"])
    output:
        R1="%s/FASTQ/TRIMMED/{samples}_R1.trimmed.fastq.gz" % (config["project-folder"]),
        R2="%s/FASTQ/TRIMMED/{samples}_R2.trimmed.fastq.gz" % (config["project-folder"])
    params:
        trim_front1=config["params"]["fastp"]["trim_front1"],
        trim_tail1=config["params"]["fastp"]["trim_tail1"],
        trim_front2=config["params"]["fastp"]["trim_front2"],
        trim_tail2=config["params"]["fastp"]["trim_tail2"]
    log:
        "%s/logs/fastp_trim_reads.{samples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/fastp_trim_reads.{samples}.benchmark.tsv" % (config["project-folder"])
    resources:
        time=cluster["fastp_trim_reads"]["time"],
        mem=cluster["fastp_trim_reads"]["mem-per-cpu"]
    threads: cluster["fastp_trim_reads"]["cpus-per-task"]
    singularity: config["singularity"]["fastp"]
    shell:"""
        fastp -i {input.R1} -I {input.R2} -o {output.R1} -O {output.R2} --trim_front1 {params.trim_front1} --trim_tail1 {params.trim_tail1} --trim_front2 {params.trim_front2} --trim_tail2 {params.trim_tail2} --thread {threads} &> {log}
	"""
	
rule concatenate_inputfiles_cat:
    """
    Concatenate input files (CAT).
    """
    input:
        r1=merge_r1_reads,
        r2=merge_r2_reads
    output:
        r1="%s/FASTQ/MERGED/all_merged_R1.fastq.gz" % (config["project-folder"]),
        r2="%s/FASTQ/MERGED/all_merged_R2.fastq.gz" % (config["project-folder"])
    benchmark:
        "%s/benchmark/concatenate_inputfiles_cat.benchmark.tsv" % (config["project-folder"])
    resources:
        time=cluster["concatenate_inputfiles_cat"]["time"],
        mem=cluster["concatenate_inputfiles_cat"]["mem-per-cpu"]
    threads: cluster["concatenate_inputfiles_cat"]["cpus-per-task"]
    shell:"""
       # In case the other approach crashes, bring it back to this way
       #zcat {input.r1} | gzip -c > {output.r1};
       #zcat {input.r2} | gzip -c > {output.r2};

       cat {input.r1} > {output.r1};
       cat {input.r2} > {output.r2};

    """

