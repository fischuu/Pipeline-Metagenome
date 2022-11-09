rule fastqc_quality_control_raw_data_r1:
    """
    Quality control of lane-wise fastq files (FASTQC).
    """
    input:
        get_raw_input_read_bs1
    output:
        "%s/QC/RAW/{reads1}_fastqc.zip" % (config["project-folder"])
    log:
        "%s/logs/fastqc_quality_control_raw_data_r1.{reads1}.log" % (config["project-folder"]),
    benchmark:
        "%s/benchmark/fastqc_quality_control_raw_data_r1.{reads1}.tsv" % (config["project-folder"])
    threads: cluster["fastqc_quality_control_raw_data_r1"]["cpus-per-task"]
    resources:
        time=cluster["fastqc_quality_control_raw_data_r1"]["time"],
        mem=cluster["fastqc_quality_control_raw_data_r1"]["mem-per-cpu"]
    params:
        outfolder="%s/QC/RAW/" % (config["project-folder"])
    singularity: config["singularity"]["gbs"]
    shell:"""
        mkdir -p {params.outfolder};
        fastqc -t {threads} -o {params.outfolder} --extract {input} &> {log};
    """

rule multiqc_quality_control_raw_data_r1:
    """
    Quality control of lane-wise fastq files in lr1(MULTIQC).
    """
    input:
        expand("%s/QC/RAW/{reads1}_fastqc.zip" % (config["project-folder"]), reads1=reads1_trim),
    output:
        directory("%s/QC/RAW/multiqc_R1/" % (config["project-folder"])),
    log:
        "%s/logs/multiqc_quality_control_raw_data_r1.log" % (config["project-folder"]),
    benchmark:
        "%s/benchmark/multiqc_quality_control_raw_data_r1.tsv" % (config["project-folder"])
    params:
       zip= lambda wildcards: "%s/QC/RAW/{wildcards.reads1}_fastqc.zip" % (config["project-folder"]),
       tmpdir=config["tmpdir"]
    threads: cluster["multiqc_quality_control_raw_data_r1"]["cpus-per-task"]
    resources:
        time=cluster["multiqc_quality_control_raw_data_r1"]["time"],
        mem=cluster["multiqc_quality_control_raw_data_r1"]["mem-per-cpu"]
    singularity: config["singularity"]["gbs"]
    shell:"""
        export TMPDIR={params.tmpdir}
        echo $TMPDIR
        multiqc -f --interactive -o {output} {input} &> {log};
    """

rule qualDist_raw_data:
    """
    Quality control of fastq files (FASTQC).
    """
    input:
        R1=get_raw_input_read1,
        R2=get_raw_input_read2
    output:
        R1="%s/QC/RAW/{rawsamples}_R1_qualdist.txt" % (config["project-folder"]),
        R2="%s/QC/RAW/{rawsamples}_R2_qualdist.txt" % (config["project-folder"])
    benchmark:
        "%s/benchmark/qualDist_raw_data.{rawsamples}.tsv" % (config["project-folder"])
    params:
        outfolder="%s/QC/RAW/" % (config["project-folder"]),
        pipeFolder=config["pipeline-folder"]
    shell:"""
        mkdir -p {params.outfolder};
        {params.pipeFolder}/scripts/getQualDist.sh {input.R1} > {output.R1}
        {params.pipeFolder}/scripts/getQualDist.sh {input.R2} > {output.R2}
    """
    
rule fastqc_quality_control_concatenated_data:
    """
    Quality control of trimmed fastq files (FASTQC).
    """
    input:
        R1="%s/FASTQ/CONCATENATED/{samples}_R1.fastq.gz" % (config["project-folder"]),
        R2="%s/FASTQ/CONCATENATED/{samples}_R2.fastq.gz" % (config["project-folder"])
    output:
        R1="%s/QC/CONCATENATED/{samples}_R1_fastqc.zip" % (config["project-folder"]),
        R2="%s/QC/CONCATENATED/{samples}_R2_fastqc.zip" % (config["project-folder"])
    log:
        R1="%s/logs/fastqc_quality_control_concatenated_data_R1.{samples}.log" % (config["project-folder"]),
        R2="%s/logs/fastqc_quality_control_concatenated_data_R2.{samples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/fastqc_quality_control_concatenated_data.{samples}.tsv" % (config["project-folder"])
    threads: cluster["fastqc_quality_control_concatenated_data_r1"]["cpus-per-task"]
    resources:
        time=cluster["fastqc_quality_control_concatenated_data_r1"]["time"],
        mem=cluster["fastqc_quality_control_concatenated_data_r1"]["mem-per-cpu"]
    params:
        outfolder="%s/QC/CONCATENATED/" % (config["project-folder"])
    singularity: config["singularity"]["gbs"]
    shell:"""
        mkdir -p {params.outfolder};
        fastqc -t {threads} -o {params.outfolder} --extract {input.R1} &> {log.R1};
        fastqc -t {threads} -o {params.outfolder} --extract {input.R2} &> {log.R2};
    """
    
rule multiqc_quality_control_concatenated_data:
    """
    Quality control of trimmed fastq files(MULTIQC).
    """
    input:
        R1=expand("%s/QC/CONCATENATED/{samples}_R1_fastqc.zip" % (config["project-folder"]), samples=samples),
        R2=expand("%s/QC/CONCATENATED/{samples}_R2_fastqc.zip" % (config["project-folder"]), samples=samples),
    output:
        R1=directory("%s/QC/CONCATENATED/multiqc_R1/" % (config["project-folder"])),
        R2=directory("%s/QC/CONCATENATED/multiqc_R2/" % (config["project-folder"]))
    log:
        R1="%s/logs/multiqc_quality_control_concatenated_data_R1.log" % (config["project-folder"]),
        R2="%s/logs/multiqc_quality_control_concatenated_data_R2.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/multiqc_quality_control_concatenated_data.tsv" % (config["project-folder"])
    params:
       R1="%s/QC/CONCATENATED/*_R1_fastqc.zip" % (config["project-folder"]),
       R2="%s/QC/CONCATENATED/*_R2_fastqc.zip" % (config["project-folder"]),
       tmpdir=config["tmpdir"]
    singularity: config["singularity"]["gbs"]
    shell:"""
        export TMPDIR={params.tmpdir}
        echo $TMPDIR
        multiqc -f --interactive -o {output.R1} {params.R1} &> {log.R1};
        multiqc -f --interactive -o {output.R2} {params.R2} &> {log.R2};
    """
    
rule qualDist_concateated_data:
    """
    Quality control of fastq files (FASTQC).
    """
    input:
        R1="%s/FASTQ/CONCATENATED/{samples}_R1.fastq.gz" % (config["project-folder"]),
        R2="%s/FASTQ/CONCATENATED/{samples}_R2.fastq.gz" % (config["project-folder"])
    output:
        R1="%s/QC/CONCATENATED/{samples}_R1_qualdist.txt" % (config["project-folder"]),
        R2="%s/QC/CONCATENATED/{samples}_R2_qualdist.txt" % (config["project-folder"])
    benchmark:
        "%s/benchmark/qualDist_concateated_data.{samples}.tsv" % (config["project-folder"])
    params:
        outfolder="%s/QC/CONCATENATED/" % (config["project-folder"]),
        pipeFolder=config["pipeline-folder"]
    shell:"""
        mkdir -p {params.outfolder};
        {params.pipeFolder}/scripts/getQualDist.sh {input.R1} > {output.R1}
        {params.pipeFolder}/scripts/getQualDist.sh {input.R2} > {output.R2}
    """
  
rule fastqc_quality_control_trimmed_data:
    """
    Quality control of fastq files (FASTQC).
    """
    input:
        R1="%s/FASTQ/TRIMMED/{samples}_R1.trimmed.fastq.gz" % (config["project-folder"]),
        R2="%s/FASTQ/TRIMMED/{samples}_R2.trimmed.fastq.gz" % (config["project-folder"])
    output:
        R1="%s/QC/TRIMMED/{samples}.R1.trimmed_fastqc.zip" % (config["project-folder"]),
        R2="%s/QC/TRIMMED/{samples}.R2.trimmed_fastqc.zip" % (config["project-folder"])
    log:
        R1="%s/logs/fastqc_quality_control_trimmed_data_R1.{samples}.log" % (config["project-folder"]),
        R2="%s/logs/fastqc_quality_control_trimmed_data_R2.{samples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/fastqc_quality_control_trimmed_data.{samples}.tsv" % (config["project-folder"])
    threads: cluster["fastqc_quality_control_trimmed_data_r1"]["cpus-per-task"]
    resources:
        time=cluster["fastqc_quality_control_trimmed_data_r1"]["time"],
        mem=cluster["fastqc_quality_control_trimmed_data_r1"]["mem-per-cpu"]
    params:
        outfolder="%s/QC/TRIMMED/" % (config["project-folder"])
    conda:"envs/gbs.yaml"
    singularity: config["singularity"]["gbs"]
    shell:"""
        mkdir -p {params.outfolder};
        fastqc -t {threads} -o {params.outfolder} --extract {input.R1} &> {log.R1};
        fastqc -t {threads} -o {params.outfolder} --extract {input.R2} &> {log.R2};
    """
    
rule multiqc_quality_control_trimmed_data:
    """
    Quality control of trimmed fastq files(MULTIQC).
    """
    input:
        R1=expand("%s/QC/TRIMMED/{samples}.R1.trimmed_fastqc.zip" % (config["project-folder"]), samples=samples),
        R2=expand("%s/QC/TRIMMED/{samples}.R2.trimmed_fastqc.zip" % (config["project-folder"]), samples=samples),
        #R1="%s/QC/TRIMMED/*.R1_fastqc.zip" % (config["project-folder"]),
        #R2="%s/QC/TRIMMED/*.R2_fastqc.zip" % (config["project-folder"])
    output:
        R1=directory("%s/QC/TRIMMED/multiqc_R1/" % (config["project-folder"])),
        R2=directory("%s/QC/TRIMMED/multiqc_R2/" % (config["project-folder"]))
    log:
        R1="%s/logs/multiqc_quality_control_trimmed_data_R1.log" % (config["project-folder"]),
        R2="%s/logs/multiqc_quality_control_trimmed_data_R2.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/multiqc_quality_control_trimmed_data.tsv" % (config["project-folder"])
    params:
       R1="%s/QC/TRIMMED/*.R1_fastqc.zip" % (config["project-folder"]),
       R2="%s/QC/TRIMMED/*.R2_fastqc.zip" % (config["project-folder"]),
       tmpdir=config["tmpdir"]
    singularity: config["singularity"]["gbs"]
    shell:"""
        export TMPDIR={params.tmpdir}
        echo $TMPDIR
        multiqc -f --interactive -o {output.R1} {params.R1} &> {log.R1};
        multiqc -f --interactive -o {output.R2} {params.R2} &> {log.R2};
    """
    
rule qualDist_trimmed_data:
    """
    Quality control of fastq files (FASTQC).
    """
    input:
        R1="%s/FASTQ/TRIMMED/{samples}_R1.trimmed.fastq.gz" % (config["project-folder"]),
        R2="%s/FASTQ/TRIMMED/{samples}_R2.trimmed.fastq.gz" % (config["project-folder"])
    output:
        R1="%s/QC/TRIMMED/{samples}_R1_qualdist.txt" % (config["project-folder"]),
        R2="%s/QC/TRIMMED/{samples}_R2_qualdist.txt" % (config["project-folder"])
    benchmark:
        "%s/benchmark/qualDist_trimmed_data.{samples}.tsv" % (config["project-folder"])
    params:
        outfolder="%s/QC/TRIMMED/" % (config["project-folder"]),
        pipeFolder=config["pipeline-folder"]
    shell:"""
        mkdir -p {params.outfolder};
        {params.pipeFolder}/scripts/getQualDist.sh {input.R1} > {output.R1}
        {params.pipeFolder}/scripts/getQualDist.sh {input.R2} > {output.R2}
    """
