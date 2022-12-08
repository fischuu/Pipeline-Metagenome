# Create the bowtie2 indices for the metagenomes

rule create_index_bowtie2_full:
    """
    Create the reference index (BOWTIE 2).
    """
    input:
        "%s/MEGAHIT/final.contigs.fa" % (config["project-folder"])
    output:
        "%s/MEGAHIT/bowtie2.indices/final.contigs/final.contigs_full.1.bt2l" % (config["project-folder"])
    log:
        "%s/logs/create_index_bowtie2_full.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/create_index_bowtie2_full.benchmark.tsv" % (config["project-folder"])
    params:
        out="%s/MEGAHIT/bowtie2.indices/final.contigs/final.contigs_full" % (config["project-folder"])
    resources:
        time=cluster["create_index_bowtie2"]["time"],
        mem=cluster["create_index_bowtie2"]["mem-per-cpu"]
    threads: cluster["create_index_bowtie2"]["cpus-per-task"]
    singularity: config["singularity"]["bowtie2"]
    shell:"""
        bowtie2-build --threads {threads} {input} {params.out}
    """
    
rule create_index_bowtie2_coas:
    """
    Create the reference index (BOWTIE 2).
    """
    input:
        "%s/MEGAHIT/final.contigs.group_{cagroup}.fa" % (config["project-folder"])
    output:
        "%s/MEGAHIT/bowtie2.indices/final.contigs_coas{cagroup}/final.contigs_coas{cagroup}.1.bt2l" % (config["project-folder"])
    log:
        "%s/logs/create_index_bowtie2_coas{cagroup}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/create_index_bowtie2_coas{cagroup}.benchmark.tsv" % (config["project-folder"])
    params:
        out="%s/MEGAHIT/bowtie2.indices/final.contigs_coas{cagroup}/final.contigs_coas{cagroup}" % (config["project-folder"])
    resources:
        time=cluster["create_index_bowtie2"]["time"],
        mem=cluster["create_index_bowtie2"]["mem-per-cpu"]
    threads: cluster["create_index_bowtie2"]["cpus-per-task"]
    singularity: config["singularity"]["bowtie2"]
    shell:"""
        bowtie2-build --threads {threads} {input} {params.out}
    """