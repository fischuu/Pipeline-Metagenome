rule prodigal_full_gene_prediction:
    """
    Predict genes from full assembly (PRODIGAL).
    """
    input:
        fa="%s/MEGAHIT/final.contigs.fa" % (config["project-folder"])
    output:
        fa="%s/PRODIGAL/final.contigs_full/final.contigs_full.prodigal.fa" % (config["project-folder"]),
        gtf="%s/PRODIGAL/final.contigs_full/final.contigs_full.prodigal.gtf" % (config["project-folder"]),
    log:
        "%s/logs/prodigal_full_gene_prediction.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/prodigal_full_gene_prediction.tsv" % (config["project-folder"])
    singularity: config["singularity"]["prodigal"]
    resources:
        time=cluster["prodigal_full_gene_prediction"]["time"],
        mem=cluster["prodigal_full_gene_prediction"]["mem-per-cpu"]
    threads: cluster["prodigal_full_gene_prediction"]["cpus-per-task"]
    shell:""" 
         prodigal -i {input.fa} \
                  -o {output.gtf} \
                  -a {output.fa} \
                  -p meta -f gff    
    """
    
rule prodigal_coas_gene_prediction:
    """
    Predict genes from full assembly (PRODIGAL).
    """
    input:
        fa="%s/MEGAHIT/final.contigs.group_{cagroup}.fa" % (config["project-folder"])
    output:
        fa="%s/PRODIGAL/final.contigs_group_{cagroup}/final.contigs_group_{cagroup}.prodigal.fa" % (config["project-folder"]),
        gtf="%s/PRODIGAL/final.contigs_group_{cagroup}/final.contigs_group_{cagroup}.prodigal.gtf" % (config["project-folder"]),
    log:
        "%s/logs/prodigal_group_{cagroup}_gene_prediction.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/prodigal_group_{cagroup}_gene_prediction.tsv" % (config["project-folder"])
    singularity: config["singularity"]["prodigal"]
    resources:
        time=cluster["prodigal_full_gene_prediction"]["time"],
        mem=cluster["prodigal_full_gene_prediction"]["mem-per-cpu"]
    threads: cluster["prodigal_full_gene_prediction"]["cpus-per-task"]
    shell:""" 
         prodigal -i {input.fa} \
                  -o {output.gtf} \
                  -a {output.fa} \
                  -p meta -f gff    
    """
    
checkpoint cut_prodigal_full_bash:
    """
    Cut prodigal output into smaller contigs (BASH).
    """
    input:
        "%s/PRODIGAL/final.contigs_full/final.contigs_full.prodigal.fa" % (config["project-folder"])
    output:
        directory("%s/PRODIGAL/final.contigs_full/Chunks/" % (config["project-folder"]))
    params:
        out="%s/PRODIGAL/final.contigs_full/Chunks/final.contigs_full.prodigal.chunk" % (config["project-folder"]),
        folder=config["pipeline-folder"],
        split=1000000
    log:
        "%s/logs/cut_prodigal_full_bash.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/cut_prodigal_full_bash.tsv" % (config["project-folder"])
    resources:
        time=cluster["cut_prodigal_full_bash"]["time"],
        mem=cluster["cut_prodigal_full_bash"]["mem-per-cpu"]
    threads: cluster["cut_prodigal_full_bash"]["cpus-per-task"]
    shell:"""
       mkdir -p {output}
       {params.folder}/scripts/cutProdigal.sh {params.split} {params.out} {input}
    """
    
checkpoint cut_prodigal_coas_bash:
    """
    Cut prodigal output into smaller contigs (BASH).
    """
    input:
        "%s/PRODIGAL/final.contigs_group_{cagroup}/final.contigs_group_{cagroup}.prodigal.fa" % (config["project-folder"])
    output:
        directory("%s/PRODIGAL/final.contigs_group_{cagroup}/Chunks/" % (config["project-folder"]))
    params:
        out="%s/PRODIGAL/final.contigs_group_{cagroup}/Chunks/final.contigs_group_{cagroup}.prodigal.chunk" % (config["project-folder"]),
        folder=config["pipeline-folder"],
        split=1000000
    log:
        "%s/logs/cut_prodigal_group_{cagroup}_bash.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/cut_prodigal_group_{cagroup}_bash.tsv" % (config["project-folder"])
    resources:
        time=cluster["cut_prodigal_coas_bash"]["time"],
        mem=cluster["cut_prodigal_coas_bash"]["mem-per-cpu"]
    threads: cluster["cut_prodigal_coas_bash"]["cpus-per-task"]
    shell:"""
       mkdir -p {output}
       {params.folder}/scripts/cutProdigal.sh {params.split} {params.out} {input}
    """

rule prepare_eggnog_database:
    """
    Prepare the database for eggnog (EGGNOG).
    """
    output:
        directory("%s/EGGNOG/DATA/" % (config["project-folder"]))
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

rule eggnog_find_homology_parallel_full:
    """
    Find homolog genes in data (EGGNOG).
    """
    input:
        folder="%s/EGGNOG/DATA/" % (config["project-folder"]),
        files="%s/PRODIGAL/final.contigs_full/Chunks/final.contigs_full.prodigal.chunk.{i}" % (config["project-folder"])
    output:
        file="%s/PRODIGAL/final.contigs_full/Chunks/chunk.{i}.emapper.seed_orthologs" % (config["project-folder"])
    log:
        "%s/logs/eggnog_find_homology_parallel_full.{i}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/eggnog_find_homology_parallel_full.{i}.tsv" % (config["project-folder"])
    params:
        folder="%s/PRODIGAL/final.contigs_full/Chunks/" % (config["project-folder"]),
        tmp=config["local-scratch"],
        out="chunk.{i}",
        fa="%s/EGGNOG/DATA/eggnog*" % (config["project-folder"])
    resources:
        time=cluster["eggnog_find_homology_parallel_full"]["time"],
        mem=cluster["eggnog_find_homology_parallel_full"]["mem-per-cpu"]
    threads: cluster["eggnog_find_homology_parallel_full"]["cpus-per-task"]
    singularity: config["singularity"]["eggnog"]
    shell:""" 
         cp {params.fa} {params.tmp}  &>> {log};
         emapper.py -m diamond --data_dir {params.tmp} --no_annot --no_file_comments --cpu {threads} -i {input.files} --output_dir {params.folder} -o {params.out}  &>> {log};
    """
    
rule eggnog_find_homology_parallel_coas:
    """
    Find homolog genes in data (EGGNOG).
    """
    input:
        folder="%s/EGGNOG/DATA/" % (config["project-folder"]),
        files="%s/PRODIGAL/final.contigs_group_{cagroup}/Chunks/final.contigs_group_{cagroup}.prodigal.chunk.{i}" % (config["project-folder"])
    output:
        file="%s/PRODIGAL/final.contigs_group_{cagroup}/Chunks/chunk.{i}.emapper.seed_orthologs" % (config["project-folder"])
    log:
        "%s/logs/eggnog_find_homology_parallel_group_{cagroup}.{i}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/eggnog_find_homology_parallel_group_{cagroup}.{i}.tsv" % (config["project-folder"])
    params:
        folder="%s/PRODIGAL/final.contigs_group_{cagroup}/Chunks/" % (config["project-folder"]),
        tmp=config["local-scratch"],
        out="chunk.{i}",
        fa="%s/EGGNOG/DATA/eggnog*" % (config["project-folder"])
    resources:
        time=cluster["eggnog_find_homology_parallel_coas"]["time"],
        mem=cluster["eggnog_find_homology_parallel_coas"]["mem-per-cpu"]
    threads: cluster["eggnog_find_homology_parallel_coas"]["cpus-per-task"]
    singularity: config["singularity"]["eggnog"]
    shell:""" 
         cp {params.fa} {params.tmp}  &>> {log};
         emapper.py -m diamond --data_dir {params.tmp} --no_annot --no_file_comments --cpu {threads} -i {input.files} --output_dir {params.folder} -o {params.out}  &>> {log};
    """
    
def aggregate_full_eggnog_search(wildcards):
    """
    Aggregate the input object for the eggnog search
    """
    checkpoint_outputEGG = checkpoints.cut_prodigal_full_bash.get(**wildcards).output[0]
    return expand("%s/PRODIGAL/final.contigs_full/Chunks/chunk.{i}.emapper.seed_orthologs" % (config["project-folder"]),
                  i=glob_wildcards(os.path.join(checkpoint_outputEGG, "final.contigs.prodigal.chunk.{i}")).i)        

rule aggregate_full_eggnog:
    input:
        aggregate_full_eggnog_search
    output:
        "%s/PRODIGAL/final.contigs_full/eggnog.emapper.seed_orthologs" % (config["project-folder"])
    log:
        "%s/logs/aggregate_full_eggnog.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/aggregate_full_eggnog.benchmark.tsv" % (config["project-folder"])
    shell:"""
       cat {input} > {output} 2> {log}
    """
    
def aggregate_coas_eggnog_search(wildcards):
    """
    Aggregate the input object for the eggnog search
    """
    checkpoint_outputEGG = checkpoints.cut_prodigal_coas_bash.get(**wildcards).output[0]
    return expand("%s/PRODIGAL/final.contigs_group_{cagroup}/Chunks/chunk.{i}.emapper.seed_orthologs" % (config["project-folder"]),
                  i=glob_wildcards(os.path.join(checkpoint_outputEGG, "final.contigs.prodigal.chunk.{i}")).i)        

rule aggregate_coas_eggnog:
    input:
        aggregate_coas_eggnog_search,
    output:
        "%s/PRODIGAL/final.contigs_group_{cagroup}/eggnog.emapper.seed_orthologs" % (config["project-folder"])
    log:
        "%s/logs/aggregate_group_{cagroup}_eggnog.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/aggregate_group_{cagroup}_eggnog.benchmark.tsv" % (config["project-folder"])
    shell:"""
       cat {input} > {output} 2> {log}
    """
    
rule eggnog_full_orthology:
    """
    Find orthology and annotate (EGGNOG).
    """
    input:
        orthologs="%s/PRODIGAL/final.contigs_full/eggnog.emapper.seed_orthologs" % (config["project-folder"]),
        folder="%s/EGGNOG/DATA/" % (config["project-folder"])
    output:
        "%s/EGGNOG/final.contigs_full/eggnog_output.emapper.annotations" % (config["project-folder"])
    log:
        "%s/logs/eggnog_full_orthology.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/eggnog_full_orthology.benchmark.tsv" % (config["project-folder"])
    params:
        tmp=config["local-scratch"],
        fa="%s/EGGNOG/DATA/eggnog*" % (config["project-folder"]),
        out="%s/EGGNOG/final.contigs_full/eggnog_output" % (config["project-folder"])
    resources:
        time=cluster["eggnog_full_orthology"]["time"],
        mem=cluster["eggnog_full_orthology"]["mem-per-cpu"]
    threads: cluster["eggnog_full_orthology"]["cpus-per-task"]
    singularity: config["singularity"]["eggnog"]
    shell:"""
    # Actually, not sure if that makes any sense, I think the files should not be concatenated here, but treated still here concatenated and then rather be merged afterwards
    # For now it is alright, as we annotated approx 650 per seconds (=run takes 24h), but if you ever rerun this step again, change it that it will be processed on the chunks
    # and then merge the chunks!
       cp {params.fa} {params.tmp}  &>> {log};
       emapper.py --data_dir {params.tmp} --annotate_hits_table {input.orthologs} --no_file_comments -o {params.out} --cpu {threads} &> {log}
    """
    
rule eggnog_coas_orthology:
    """
    Find orthology and annotate (EGGNOG).
    """
    input:
        orthologs="%s/PRODIGAL/final.contigs_group_{cagroup}/eggnog.emapper.seed_orthologs" % (config["project-folder"]),
        folder="%s/EGGNOG/DATA/" % (config["project-folder"])
    output:
        "%s/EGGNOG/final.contigs_group_{cagroup}/eggnog_output.emapper.annotations" % (config["project-folder"])
    log:
        "%s/logs/eggnog_group_{cagroup}_orthology.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/eggnog_group_{cagroup}_orthology.benchmark.tsv" % (config["project-folder"])
    params:
        tmp=config["local-scratch"],
        fa="%s/EGGNOG/DATA/eggnog*" % (config["project-folder"]),
        out="%s/EGGNOG/final.contigs_group_{cagroup}/eggnog_output" % (config["project-folder"])
    resources:
        time=cluster["eggnog_coas_orthology"]["time"],
        mem=cluster["eggnog_coas_orthology"]["mem-per-cpu"]
    threads: cluster["eggnog_coas_orthology"]["cpus-per-task"]
    singularity: config["singularity"]["eggnog"]
    shell:"""
    # Actually, not sure if that makes any sense, I think the files should not be concatenated here, but treated still here concatenated and then rather be merged afterwards
    # For now it is alright, as we annotated approx 650 per seconds (=run takes 24h), but if you ever rerun this step again, change it that it will be processed on the chunks
    # and then merge the chunks!
       cp {params.fa} {params.tmp}  &>> {log};
       emapper.py --data_dir {params.tmp} --annotate_hits_table {input.orthologs} --no_file_comments -o {params.out} --cpu {threads} &> {log}
    """
    
