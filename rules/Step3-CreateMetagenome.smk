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
    
rule filter_length_megahit:
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
        time=cluster["filter_length_megahit"]["time"],
        mem=cluster["filter_length_megahit"]["mem-per-cpu"]
    threads: cluster["filter_length_megahit"]["cpus-per-task"]
    shell:"""
        {params.folder}/scripts/filterLengthFasta.sh 1000 {input} {output.filter1k}
        {params.folder}/scripts/filterLengthFasta.sh 2000 {input} {output.filter2k}    
    """
    
rule create_index_bowtie2:
    """
    Create the reference index (BOWTIE 2).
    """
    input:
        "%s/MEGAHIT/final.contigs.fa" % (config["project-folder"]) 
    output:
        "%s/MEGAHIT/final.contigs.fa.1.bt2l" % (config["project-folder"])
    log:
        "%s/logs/create_index_bowtie2.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/create_index_bowtie2.benchmark.tsv" % (config["project-folder"])
    singularity: config["singularity"]["bowtie2"]
    params:
        index="%s/MEGAHIT/final.contigs.fa" % (config["project-folder"])
    resources:
        time=cluster["create_index_bowtie2"]["time"],
        mem=cluster["create_index_bowtie2"]["mem-per-cpu"]
    threads: cluster["create_index_bowtie2"]["cpus-per-task"]
    shell:"""

    # WARNING, THIS CREATES THE WRONG OUTPUT!! IT SHOULD CREATE bt2mega.* , BUT IT CREATES ONLY .* FILES!!!!
    
        bowtie2-build --threads {threads} {input} {params.index} &> {log}
    """
    
rule map_data_metagenome_bowtie2:
    """
    Map the data to metagenome(BOWTIE 2).
    """
    input:
        R1="%s/FASTQ/TRIMMED/{samples}_R1.trimmed.fastq.gz" % (config["project-folder"]),
        R2="%s/FASTQ/TRIMMED/{samples}_R2.trimmed.fastq.gz" % (config["project-folder"]),
        index="%s/MEGAHIT/final.contigs.fa.1.bt2l" % (config["project-folder"])
    output:
        temp("%s/BAM/megahit/{samples}_mega.sam" % (config["project-folder"]))
    log:
        "%s/logs/map_data_metagenome_bowtie2.{samples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/map_data_metagenome_bowtie2.{samples}.benchmark.tsv" % (config["project-folder"])
    params:
        index="%s/MEGAHIT/final.contigs.fa" % (config["project-folder"])
    singularity: config["singularity"]["bowtie2"]
    resources:
        time=cluster["map_data_metagenome_bowtie2"]["time"],
        mem=cluster["map_data_metagenome_bowtie2"]["mem-per-cpu"]
    threads: cluster["map_data_metagenome_bowtie2"]["cpus-per-task"]
    shell:"""
        bowtie2 -p {threads} -x {params.index} -1 {input.R1} -2 {input.R2} > {output}
    """
    
rule samToBam_data_metagenome:
    """
    Map the data to metagenome(BOWTIE 2).
    """
    input:
        "%s/BAM/megahit/{samples}_mega.sam" % (config["project-folder"])    
    output:
        "%s/BAM/megahit/{samples}_mega.bam" % (config["project-folder"])
    log:
        "%s/logs/samToBam_data_metagenome.{samples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/samToBam_data_metagenome.{samples}.benchmark.tsv" % (config["project-folder"])
    params:
        index="%s/MEGAHIT/final.contigs.fa" % (config["project-folder"])
    singularity: config["singularity"]["samtools"]
    resources:
        time=cluster["samToBam_data_metagenome"]["time"],
        mem=cluster["samToBam_data_metagenome"]["mem-per-cpu"]
    threads: cluster["samToBam_data_metagenome"]["cpus-per-task"]
    shell:"""
        samtools view -bS {input} | samtools sort > {output};
        
        samtools index {output}
        
        # There might be a timestamp issue for the downstream analysis, if so, run a touch *.bai on that folder before the problematic rule.
    """


##### Co-assemblies
################################################################################