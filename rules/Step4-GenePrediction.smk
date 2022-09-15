rule prodigal_gene_prediction:
    """
    Predict genes (PRODIGAL).
    """
    input:
        fa="%s/MEGAHIT/final.contigs.fa" % (config["project-folder"])
    output:
        fa="%s/PRODIGAL/final.contigs.prodigal.fa" % (config["project-folder"]),
        gtf="%s/PRODIGAL/final.contigs.prodigal.gtf" % (config["project-folder"]),
    log:
        "%s/logs/prodigal_gene_prediction.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/prodigal_gene_prediction.tsv" % (config["project-folder"])
    singularity: config["singularity"]["prodigal"]
    resources:
        time=cluster["prodigal_gene_prediction"]["time"],
        mem=cluster["prodigal_gene_prediction"]["mem-per-cpu"]
    threads: cluster["prodigal_gene_prediction"]["cpus-per-task"]
    shell:""" 
         prodigal -i {input.fa} \
                  -o {output.gtf} \
                  -a {output.fa} \
                  -p meta -f gff    
    """
    
checkpoint cut_prodigal_bash:
    """
    Cut prodigal output into smaller contigs (BASH).
    """
    input:
        "%s/PRODIGAL/final.contigs.prodigal.fa" % (config["project-folder"])
    output:
        directory("%s/PRODIGAL/Chunks/" % (config["project-folder"]))
    params:
        out="%s/PRODIGAL/Chunks/final.contigs.prodigal.chunk" % (config["project-folder"]),
        folder=config["pipeline-folder"],
        split=1000000
    log:
        "%s/logs/cut_prodigal_bash.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/cut_prodigal_bash.tsv" % (config["project-folder"])
    resources:
        time=cluster["cut_prodigal_bash"]["time"],
        mem=cluster["cut_prodigal_bash"]["mem-per-cpu"]
    threads: cluster["cut_prodigal_bash"]["cpus-per-task"]
    shell:"""
       mkdir -p {output}
       {params.folder}/scripts/cutProdigal.sh {params.split} {params.out} {input}
    """

rule prepare_eggnog_database:
    """
    Prepare the database for eggnog (EGGNOG).
    """
    output:
        directory("%s/PRODIGAL/EGGNOG-DATA/" % (config["project-folder"]))
    log:
        "%s/logs/prepare_eggnog_database.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/prepare_eggnog_database.tsv" % (config["project-folder"])
    resources:
        time=cluster["prepare_eggnog_database"]["time"],
        mem=cluster["prepare_eggnog_database"]["mem-per-cpu"]
    threads: cluster["prepare_eggnog_database"]["cpus-per-task"]
    singularity: config["singularity"]["eggnog"]
    shell:""" 
         mkdir -p {output}
         download_eggnog_data.py -y --data_dir {output}  &> {log};
    """

rule eggnog_find_homology_parallel:
    """
    Find homolog genes in data (EGGNOG).
    """
    input:
        folder="%s/PRODIGAL/EGGNOG-DATA/" % (config["project-folder"]),
        files="%s/PRODIGAL/Chunks/final.contigs.prodigal.chunk.{i}" % (config["project-folder"])
    output:
        "%s/PRODIGAL/Chunks/{i}.emapper.seed_orthologs" % (config["project-folder"])
    log:
        "%s/logs/eggnog_find_homology_parallel.{i}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/eggnog_find_homology_parallel.{i}.tsv" % (config["project-folder"])
    params:
        tmp=config["local-scratch"],
        fa="%s/PRODIGAL/EGGNOG-DATA/eggnog*" % (config["project-folder"]),
        out="%s/PRODIGAL/Chunks/" % (config["project-folder"])
    resources:
        time=cluster["eggnog_find_homology_parallel"]["time"],
        mem=cluster["eggnog_find_homology_parallel"]["mem-per-cpu"]
    threads: cluster["eggnog_find_homology_parallel"]["cpus-per-task"]
    singularity: config["singularity"]["eggnog"]
    shell:""" 
         cp {params.fa} {params.tmp}  &>> {log};
         emapper.py -m diamond --data_dir {params.tmp} --no_annot --no_file_comments --cpu {threads} -i {input.files} -o {params.out}  &>> {log};
    """
    
def aggregate_eggnog_search(wildcards):
    """
    Aggregate the input object for the eggnog search
    """
    checkpoint_outputEGG = checkpoints.cut_prodigal_bash.get(**wildcards).output[0]
    return expand("%s/PRODIGAL/Chunks/{i}.emapper.seed_orthologs" % (config["project-folder"]),
                  i=glob_wildcards(os.path.join(checkpoint_outputEGG, "final.contigs.prodigal.chunk.{i}")).i)        

rule aggregate_eggnog:
    input:
        aggregate_eggnog_search
    output:
        "%s/PRODIGAL/input_file.emapper.seed_orthologs" % (config["project-folder"])
    log:
        "%s/logs/concatenate.homologs.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/concatenate_homologs.benchmark.tsv" % (config["project-folder"])
    shell:"""
       cat {input} > {output}
    """