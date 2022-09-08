rule star_create_host_index:
    """
    Create genome index for host genome (STAR).
    """
    input:
        fasta="%s" % (config["host"])
    output:
        "%s/chrName.txt" % (config["host-index"])
    log:
        "%s/logs/star_create_host_index.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/star_create_host_index.benchmark.tsv" % (config["project-folder"])
    params:
        hostindex=config["host-index"]
    singularity: config["singularity"]["star"]
    resources:
        time=cluster["star_create_host_index"]["time"],
        mem=cluster["star_create_host_index"]["mem-per-cpu"]
    threads: cluster["star_create_host_index"]["cpus-per-task"]
    shell:"""
        mkdir -p {params.hostindex}
        
        echo "Number of threads used:" {threads}
        
        STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {input.index} --genomeFastaFiles {input.fasta} &> {log}
    """