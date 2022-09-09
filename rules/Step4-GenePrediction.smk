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
        out="%s/PRODIGAL/final.contigs.prodigal.chunk" % (config["project-folder"]),
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
       ./scripts/cutProdigal.sh {params.split} {params.pre} {input}
    """

#############################################################################
## Above should be alrigfht, below I need to adjust still, check from the GBS
## Pipeline how the checkpoints are used with rules.

rule eggnog_find_homology:
    """
    Find homolog genes in data (EGGNOG).
    """
    input:
        "%s/prodigal_out/{prodigalChunks}" % (config["project-folder"])
    output:
        "%s/prodigal_out/{prodigalChunks}.emapper.seed_orthologs" % (config["project-folder"])
    log:
        "%s/logs/EGGNOG/find_homologs_{prodigalChunks}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/EGGNOG/find_homologs_{prodigalChunks}.tsv" % (config["project-folder"])
    params:
        tmp=config["local-scratch"],
        fa="tools/eggnog-mapper/data/eggnog*"
    shell:""" 
         cp {params.fa} {params.tmp}  &>> {log};
         emapper.py -m diamond --data_dir {params.tmp} --no_annot --no_file_comments --cpu {params.threads} -i {input} -o {input}  &>> {log};
    """