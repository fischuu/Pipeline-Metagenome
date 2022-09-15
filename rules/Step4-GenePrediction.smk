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
        file="%s/PRODIGAL/Chunks/chunk.{i}.emapper.seed_orthologs" % (config["project-folder"])
    log:
        "%s/logs/eggnog_find_homology_parallel.{i}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/eggnog_find_homology_parallel.{i}.tsv" % (config["project-folder"])
    params:
        folder="%s/PRODIGAL/Chunks/" % (config["project-folder"]),
        tmp=config["local-scratch"],
        out="chunk.{i}",
        fa="%s/PRODIGAL/EGGNOG-DATA/eggnog*" % (config["project-folder"])
    resources:
        time=cluster["eggnog_find_homology_parallel"]["time"],
        mem=cluster["eggnog_find_homology_parallel"]["mem-per-cpu"]
    threads: cluster["eggnog_find_homology_parallel"]["cpus-per-task"]
    singularity: config["singularity"]["eggnog"]
    shell:""" 
         cp {params.fa} {params.tmp}  &>> {log};
         emapper.py -m diamond --data_dir {params.tmp} --no_annot --no_file_comments --cpu {threads} -i {input.files} --output_dir {params.folder} -o {params.out}  &>> {log};
    """
    
def aggregate_eggnog_search(wildcards):
    """
    Aggregate the input object for the eggnog search
    """
    checkpoint_outputEGG = checkpoints.cut_prodigal_bash.get(**wildcards).output[0]
    return expand("%s/PRODIGAL/Chunks/chunk.{i}.emapper.seed_orthologs" % (config["project-folder"]),
                  i=glob_wildcards(os.path.join(checkpoint_outputEGG, "final.contigs.prodigal.chunk.{i}")).i)        

rule aggregate_eggnog:
    input:
        aggregate_eggnog_search
    output:
        "%s/PRODIGAL/eggnog.emapper.seed_orthologs" % (config["project-folder"])
    log:
        "%s/logs/concatenate.homologs.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/concatenate_homologs.benchmark.tsv" % (config["project-folder"])
    shell:"""
       cat {input} > {output} 2> {log}
    """
    
rule eggnog_orthology:
    """
    Find orthology and annotate (EGGNOG).
    """
    input:
        orthologs="%s/PRODIGAL/eggnog.emapper.seed_orthologs" % (config["project-folder"])
    output:
        "%s/EGGNOG/eggnog_output.emapper.annotations" % (config["project-folder"])
    log:
        "%s/logs/eggnog_orthology.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/eggnog_orthology.benchmark.tsv" % (config["project-folder"])
    params:
        tmp=config["local-scratch"],
        fa="%s/PRODIGAL/EGGNOG-DATA/eggnog*" % (config["project-folder"]),
        out="%s/EGGNOG/eggnog_output" % (config["project-folder"])
    resources:
        time=cluster["eggnog_orthology"]["time"],
        mem=cluster["eggnog_orthology"]["mem-per-cpu"]
    threads: cluster["eggnog_orthology"]["cpus-per-task"]
    singularity: config["singularity"]["eggnog"]
    shell:"""
    # Actually, not sure if that makes any sense, I think the files should not be concatenated here, but treated still here concatenated and then rather be merged afterwards
    # For now it is alright, as we annotated approx 650 per seconds (=run takes 24h), but if you ever rerun this step again, change it that it will be processed on the chunks
    # and then merge the chunks!
       cp {params.fa} {params.tmp}  &>> {log};
       emapper.py --data_dir {params.tmp} --annotate_hits_table {input.orthologs} --no_file_comments -o {params.out} --cpu {threads} &> {log}
    """
    
rule quantify_predictedGenes_featureCounts:
    """
    Quantify the mapped reads after merging (featureCounts).
    """
    input:
        bam="%s/BAM/megahit/{samples}_mega.bam" % (config["project-folder"]),
        gtf="%s/PRODIGAL/final.contigs.prodigal.gtf" % (config["project-folder"])
    output:
        file="%s/QUANTIFICATION/PRODIGAL_FC/{samples}_fc.txt" % (config["project-folder"])
    log:
        "%s/logs/quantify_predictedGenes_featureCounts.{samples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/quantify_predictedGenes_featureCounts.{samples}.benchmark.tsv" % (config["project-folder"])
    resources:
        time=cluster["quantify_predictedGenes_featureCounts"]["time"],
        mem=cluster["quantify_predictedGenes_featureCounts"]["mem-per-cpu"]
    threads: cluster["quantify_predictedGenes_featureCounts"]["cpus-per-task"]
    singularity: config["singularity"]["subread"]
    shell:"""
        featureCounts -p \
                      -T {threads} \
                      -a {input.gtf} \
                      -o {output.file} \
                      -t CDS \
                      -g ID \
                      {input.bam} 2> {log}
    """