rule map_data_final_contigs_full_bowtie2:
    """
    Map the data to metagenome(BOWTIE 2).
    """
    input:
        R1="%s/FASTQ/DECONTAMINATED/{samples}_R1.decontaminated.fastq.gz" % (config["project-folder"]),
        R2="%s/FASTQ/DECONTAMINATED/{samples}_R2.decontaminated.fastq.gz" % (config["project-folder"]),
        index="%s/MEGAHIT/bowtie2.indices/final.contigs/final.contigs_full.1.bt2l" % (config["project-folder"])
    output:
        temp("%s/BAM/final.contigs_full/{samples}_mega.sam" % (config["project-folder"]))
    log:
        "%s/logs/map_data_final_contigs_full_bowtie2.{samples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/map_data_final_contigs_full_bowtie2.{samples}.benchmark.tsv" % (config["project-folder"])
    params:
        index="%s/MEGAHIT/bowtie2.indices/final.contigs/final.contigs_full" % (config["project-folder"])
    singularity: config["singularity"]["bowtie2"]
    resources:
        time=cluster["map_data_final_contigs_full_bowtie2"]["time"],
        mem=cluster["map_data_final_contigs_full_bowtie2"]["mem-per-cpu"]
    threads: cluster["map_data_final_contigs_full_bowtie2"]["cpus-per-task"]
    shell:"""
        bowtie2 -p {threads} -x {params.index} -1 {input.R1} -2 {input.R2} > {output} 2> {log};
    """
    
rule samToBam_data_metagenome:
    """
    Map the data to metagenome(BOWTIE 2).
    """
    input:
        "%s/BAM/final.contigs_full/{samples}_mega.sam" % (config["project-folder"])
    output:
        "%s/BAM/final.contigs_full/{samples}_mega.bam" % (config["project-folder"])
    log:
        "%s/logs/samToBam_data_metagenome.{samples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/samToBam_data_metagenome.{samples}.benchmark.tsv" % (config["project-folder"])
    singularity: config["singularity"]["samtools"]
    resources:
        time=cluster["samToBam_data_metagenome"]["time"],
        mem=cluster["samToBam_data_metagenome"]["mem-per-cpu"]
    threads: cluster["samToBam_data_metagenome"]["cpus-per-task"]
    shell:"""
        samtools view -bS {input} | samtools sort > {output};
        
        samtools index {output}
    """

rule map_data_final_contigs_coas_bowtie2:
    """
    Map the data to metagenome(BOWTIE 2).
    """
    input:
        R1="%s/FASTQ/DECONTAMINATED/{samples}_R1.decontaminated.fastq.gz" % (config["project-folder"]),
        R2="%s/FASTQ/DECONTAMINATED/{samples}_R2.decontaminated.fastq.gz" % (config["project-folder"]),
        index="%s/MEGAHIT/bowtie2.indices/final.contigs_coas{cagroup}/final.contigs_coas{cagroup}.1.bt2l" % (config["project-folder"])
    output:
        temp("%s/BAM/final.contigs_coas{cagroup}/{samples}_mega.sam" % (config["project-folder"]))
    log:
        "%s/logs/map_data_final_contigs_coas{cagroup}_bowtie2.{samples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/map_data_final_contigs_coas{cagroup}_bowtie2.{samples}.benchmark.tsv" % (config["project-folder"])
    params:
        index="%s/MEGAHIT/bowtie2.indices/final.contigs_coas{cagroup}/final.contigs_coas{cagroup}" % (config["project-folder"])
    singularity: config["singularity"]["bowtie2"]
    resources:
        time=cluster["map_data_final_contigs_coas_bowtie2"]["time"],
        mem=cluster["map_data_final_contigs_coas_bowtie2"]["mem-per-cpu"]
    threads: cluster["map_data_final_contigs_coas_bowtie2"]["cpus-per-task"]
    shell:"""
        bowtie2 -p {threads} -x {params.index} -1 {input.R1} -2 {input.R2} > {output} 2> {log};
    """
    
rule samToBam_data_metagenome_coas:
    """
    Map the data to metagenome(BOWTIE 2).
    """
    input:
        "%s/BAM/final.contigs_coas{cagroup}/{samples}_mega.sam" % (config["project-folder"])
    output:
        "%s/BAM/final.contigs_coas{cagroup}/{samples}_mega.bam" % (config["project-folder"])
    log:
        "%s/logs/samToBam_data_metagenome_coas{cagroup}.{samples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/samToBam_data_metagenome_coas{cagroup}.{samples}.benchmark.tsv" % (config["project-folder"])
    singularity: config["singularity"]["samtools"]
    resources:
        time=cluster["samToBam_data_metagenome"]["time"],
        mem=cluster["samToBam_data_metagenome"]["mem-per-cpu"]
    threads: cluster["samToBam_data_metagenome"]["cpus-per-task"]
    shell:"""
        samtools view -bS {input} | samtools sort > {output};
        
        samtools index {output}
    """
