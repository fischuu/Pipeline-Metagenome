rule assemble_contigs_megahit:
    """
    Denovo assemble merged reads (MEGAHIT).
    """
    input:
        r1="%s/FASTQ/MERGED/all_merged_R1.fastq.gz" % (config["project-folder"]),
        r2="%s/FASTQ/MERGED/all_merged_R2.fastq.gz" % (config["project-folder"])
    output:
        temp_dir=temp(directory("%s/MEGAHIT_temp" % (config["project-folder"]))),
        contigs="%s/MEGAHIT/final.contigs.fa" % (config["project-folder"])
    log:
        "%s/logs/assemble_contigs_megahit.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/assemble_contigs_megahit.benchmark.tsv" % (config["project-folder"])
    singularity: config["singularity"]["megahit"]
    resources:
        time=cluster["assemble_contigs_megahit"]["time"],
        mem=cluster["assemble_contigs_megahit"]["mem-per-cpu"]
    threads: cluster["assemble_contigs_megahit"]["cpus-per-task"]
    shell:"""
        megahit -1 {input.r1}  -2 {input.r2} -t {threads}  -o {output.temp_dir} &> {log};
        
        mv {output.temp_dir}/final.contigs.fa {output.contigs};
    """
    
rule assemble_coassembly_contigs_megahit:
    """
    Denovo assemble merged reads (MEGAHIT).
    """
    input:
        r1="%s/FASTQ/MERGED/group_{cagroup}_R1.fastq.gz" % (config["project-folder"]),
        r2="%s/FASTQ/MERGED/group_{cagroup}_R2.fastq.gz" % (config["project-folder"])
    output:
        temp_dir=temp(directory("%s/MEGAHIT_temp_{cagroup}" % (config["project-folder"]))),
        contigs="%s/MEGAHIT/final.contigs.group_{cagroup}.fa" % (config["project-folder"])
    log:
        "%s/logs/assemble_contigs_megahit_group_{cagroup}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/assemble_contigs_megahit_group_{cagroup}.benchmark.tsv" % (config["project-folder"])
    singularity: config["singularity"]["megahit"]
    resources:
        time=cluster["assemble_contigs_megahit"]["time"],
        mem=cluster["assemble_contigs_megahit"]["mem-per-cpu"]
    threads: cluster["assemble_contigs_megahit"]["cpus-per-task"]
    shell:"""
        megahit -1 {input.r1}  -2 {input.r2} -t {threads}  -o {output.temp_dir} &> {log};
        
        mv {output.temp_dir}/final.contigs.fa {output.contigs};
    """    
    
rule filter_length_full_megahit:
    """
    Remove short contigs from the megahit output
    """
    input:
        "%s/MEGAHIT/final.contigs.fa" % (config["project-folder"])
    output:
        filter1k="%s/MEGAHIT/final.contigs.1k.fa" % (config["project-folder"]),
        filter2k="%s/MEGAHIT/final.contigs.2k.fa" % (config["project-folder"])
    log:
        "%s/logs/filter_length_megahit.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/filter_length_megahit.benchmark.tsv" % (config["project-folder"])
    params:
        folder= config["pipeline-folder"]
    resources:
        time=cluster["filter_length_full_megahit"]["time"],
        mem=cluster["filter_length_full_megahit"]["mem-per-cpu"]
    threads: cluster["filter_length_full_megahit"]["cpus-per-task"]
    shell:"""
        {params.folder}/scripts/filterLengthFasta.sh 1000 {input} {output.filter1k}
        {params.folder}/scripts/filterLengthFasta.sh 2000 {input} {output.filter2k}    
    """
    
rule filter_length_coas_megahit:
    """
    Remove short contigs from the coas megahit output
    """
    input:
        "%s/MEGAHIT/final.contigs.group_{cagroup}.fa" % (config["project-folder"])
    output:
        filter1k="%s/MEGAHIT/final.contigs.group_{cagroup}.1k.fa" % (config["project-folder"]),
        filter2k="%s/MEGAHIT/final.contigs.group_{cagroup}.2k.fa" % (config["project-folder"])
    log:
        "%s/logs/filter_length_group_{cagroup}.megahit.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/filter_length_group_{cagroup}_megahit.benchmark.tsv" % (config["project-folder"])
    params:
        folder= config["pipeline-folder"]
    resources:
        time=cluster["filter_length_coas_megahit"]["time"],
        mem=cluster["filter_length_coas_megahit"]["mem-per-cpu"]
    threads: cluster["filter_length_coas_megahit"]["cpus-per-task"]
    shell:"""
        {params.folder}/scripts/filterLengthFasta.sh 1000 {input} {output.filter1k}
        {params.folder}/scripts/filterLengthFasta.sh 2000 {input} {output.filter2k}    
    """