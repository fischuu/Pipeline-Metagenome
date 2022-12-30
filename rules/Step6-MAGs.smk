rule cut_filtered_contigs_full_concoct:
    """
    Cut contigs larger than 20kb into smaller subcontigs of min size 10kb (CONCOCT).
    """
    input:
        "%s/MEGAHIT/final.contigs.{fafilter}.fa" % (config["project-folder"])
    output:
        bed="%s/CONCOCT/final.contigs.{fafilter}_10K.bed" % (config["project-folder"]),
        fa="%s/CONCOCT/final.contigs.{fafilter}_10K.fa" % (config["project-folder"])
    log:
        "%s/logs/cut_{fafilter}_filtered_contigs_concoct.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/cut_{fafilter}_filtered_contigs_concoct.tsv" % (config["project-folder"])
    resources:
        time=cluster["cut_filtered_contigs_concoct"]["time"],
        mem=cluster["cut_filtered_contigs_concoct"]["mem-per-cpu"]
    threads: cluster["cut_filtered_contigs_concoct"]["cpus-per-task"]
    singularity: config["singularity"]["concoct"]
    shell:"""
        cut_up_fasta.py {input} -c 10000 -o 0 --merge_last -b {output.bed} > {output.fa} 2> {log}
    """
    
rule cut_filtered_contigs_coas_concoct:
    """
    Cut contigs into smaller contigs (CONCOCT).
    """
    input:
        "%s/MEGAHIT/final.contigs.group_{cagroup}.{fafilter}.fa" % (config["project-folder"])
    output:
        bed="%s/CONCOCT/final.contigs.group_{cagroup}.{fafilter}_10K.bed" % (config["project-folder"]),
        fa="%s/CONCOCT/final.contigs.group_{cagroup}.{fafilter}_10K.fa" % (config["project-folder"])
    log:
        "%s/logs/cut_{fafilter}_filtered_contigs_group_{cagroup}_concoct.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/cut_{fafilter}_filtered_contigs_group_{cagroup}_concoct.tsv" % (config["project-folder"])
    resources:
        time=cluster["cut_filtered_contigs_coas_concoct"]["time"],
        mem=cluster["cut_filtered_contigs_coas_concoct"]["mem-per-cpu"]
    threads: cluster["cut_filtered_contigs_coas_concoct"]["cpus-per-task"]
    singularity: config["singularity"]["concoct"]
    shell:"""
        cut_up_fasta.py {input} -c 10000 -o 0 --merge_last -b {output.bed} > {output.fa} 2> {log}
    """
    
checkpoint cut_input_bed:
    """
    Divide the input bed-file for parallel processing
    """
    input:
        "%s/CONCOCT/final.contigs.{fafilter}_10K.bed" % (config["project-folder"])
    output:
        directory("%s/CONCOCT/bed_subsampled_{fafilter}" % (config["project-folder"])),
    params:
        out="%s/CONCOCT/bed_subsampled_{fafilter}/final.contigs.{fafilter}_10K." % (config["project-folder"]),
        bedcut=config["params"]["concoct"]["bedcut"]
    shell:"""
        mkdir -p {output}

        split -l {params.bedcut} --numeric-suffixes {input} {params.out}
    """
    
checkpoint cut_input_coas_bed:
    """
    Divide the input bed-file for parallel processing
    """
    input:
        "%s/CONCOCT/final.contigs.group_{cagroup}.{fafilter}_10K.bed" % (config["project-folder"])
    output:
        directory("%s/CONCOCT/bed_subsampled_group_{cagroup}_{fafilter}" % (config["project-folder"])),
    params:
        out="%s/CONCOCT/bed_subsampled_group_{cagroup}_{fafilter}/final.contigs.group_{cagroup}.{fafilter}_10K." % (config["project-folder"]),
        bedcut=config["params"]["concoct"]["bedcut"]
    shell:"""
        mkdir -p {output}

        split -l {params.bedcut} --numeric-suffixes {input} {params.out}
    """
    

rule process_covTable_bed_parallel:
    """
    Calculate the coverage for the subpart of the splitted bed
    """
    input:
        bed="%s/CONCOCT/bed_subsampled_{fafilter}/final.contigs.{fafilter}_10K.{i}" % (config["project-folder"]),
        bam=expand("%s/BAM/final.contigs_full/{samples}_mega.bam" % (config["project-folder"]), samples=samples)    
    output:
        woheader="%s/CONCOCT/coverage_table_{fafilter}_{i}.tsv" % (config["project-folder"]),
        header="%s/CONCOCT/coverage_table_{fafilter}_header_{i}.tsv" % (config["project-folder"]),
    params:
       bam="%s/BAM/final.contigs_full/*.bam" % (config["project-folder"])
    resources:
        time=cluster["process_covTable_bed_parallel"]["time"],
        mem=cluster["process_covTable_bed_parallel"]["mem-per-cpu"]
    threads: cluster["process_covTable_bed_parallel"]["cpus-per-task"]
    singularity: config["singularity"]["concoct"]
    shell:"""
      concoct_coverage_table.py {input.bed} {params.bam} > {output.woheader}
      sed -n 1p {output.woheader} > {output.header}
      
      sed -i '1d' {output.woheader}
     """
     
rule process_covTable_bed_coas_parallel:
    """
    Calculate the coverage for the subpart of the splitted bed
    """
    input:
        bed="%s/CONCOCT/bed_subsampled_group_{cagroup}_{fafilter}/final.contigs.group_{cagroup}.{fafilter}_10K.{i}" % (config["project-folder"]),
        bam=expand("%s/BAM/final.contigs_coas{cagroup}/{samples}_mega.bam" % (config["project-folder"]), samples=samples, cagroup=assemblyGroups)    
    output:
        woheader="%s/CONCOCT/coverage_table_group_{cagroup}_{fafilter}_{i}.tsv" % (config["project-folder"]),
        header="%s/CONCOCT/coverage_table_group_{cagroup}_{fafilter}_header_{i}.tsv" % (config["project-folder"]),
    params:
       bam="%s/BAM/final.contigs_group_{cagroup}/*.bam" % (config["project-folder"])
    resources:
        time=cluster["process_covTable_bed_coas_parallel"]["time"],
        mem=cluster["process_covTable_bed_coas_parallel"]["mem-per-cpu"]
    threads: cluster["process_covTable_bed_coas_parallel"]["cpus-per-task"]
    singularity: config["singularity"]["concoct"]
    shell:"""
      concoct_coverage_table.py {input.bed} {params.bam} > {output.woheader}
      sed -n 1p {output.woheader} > {output.header}
      
      sed -i '1d' {output.woheader}
     """
     
def aggregate_input_cov(wildcards):
    """
    Aggregate the input object for the final coverage table
    """
    checkpoint_output = checkpoints.cut_input_bed.get(**wildcards).output[0]
    returnValue = expand("%s/CONCOCT/coverage_table_{ff}_{i}.tsv" % (config["project-folder"]),
                  i=glob_wildcards(os.path.join(checkpoint_output, "final.contigs.{ff}_10K.{i}")).i,
                  ff={wildcards.fafilter})
    return returnValue  
    
def aggregate_input_coas_cov(wildcards):
    """
    Aggregate the input object for the final coverage table
    """
    checkpoint_output = checkpoints.cut_input_coas_bed.get(**wildcards).output[0]
    returnValue = expand("%s/CONCOCT/coverage_table_group_{gr}_{ff}_{i}.tsv" % (config["project-folder"]),
                  i=glob_wildcards(os.path.join(checkpoint_output, "final.contigs.{ff}_10K.{i}")).i,
                  ff={wildcards.fafilter},
                  gr={wildcards.cagroup})
    return returnValue  
      
rule aggregate_coverage_table:
    input:
        aggregate_input_cov
    output:
        all="%s/CONCOCT/coverage_table_{fafilter}.tsv" % (config["project-folder"]),
        tmp="%s/CONCOCT/coverage_table_{fafilter}.tsv.tmp" % (config["project-folder"]),
        tmp2="%s/CONCOCT/coverage_table_{fafilter}.tsv.tmp2" % (config["project-folder"])
    params:
        header="%s/CONCOCT/coverage_table_{fafilter}_header_00.tsv" % (config["project-folder"])
    shell:"""
        cat {input} > {output.tmp}
        sort -t_ -k2,2n {output.tmp} > {output.tmp2}
        cat {params.header} {output.tmp2} > {output.all}
    """
    
rule aggregate_coverage_table_coas:
    input:
        aggregate_input_coas_cov
    output:
        all="%s/CONCOCT/coverage_table_group_{cagroup}_{fafilter}.tsv" % (config["project-folder"]),
        tmp="%s/CONCOCT/coverage_table_group_{cagroup}_{fafilter}.tsv.tmp" % (config["project-folder"]),
        tmp2="%s/CONCOCT/coverage_table_group_{cagroup}_{fafilter}.tsv.tmp2" % (config["project-folder"])
    params:
        header="%s/CONCOCT/coverage_table_group_{cagroup}_{fafilter}_header_00.tsv" % (config["project-folder"])
    shell:"""
        cat {input} > {output.tmp}
        sort -t_ -k2,2n {output.tmp} > {output.tmp2}
        cat {params.header} {output.tmp2} > {output.all}
    """
    
rule run_concoct:
    """
    Apply concoct to the data (CONCOCT).
    """
    input:
        fa="%s/CONCOCT/final.contigs.{fafilter}_10K.fa" % (config["project-folder"]),
        coverage="%s/CONCOCT/coverage_table_{fafilter}.tsv" % (config["project-folder"])
    output:
        "%s/CONCOCT/{fafilter}_clustering_gt1000.csv" % (config["project-folder"])
    log:
        "%s/logs/run_concoct_{fafilter}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/run_concoct_{fafilter}.tsv" % (config["project-folder"])
    resources:
        time=cluster["run_concoct"]["time"],
        mem=cluster["run_concoct"]["mem-per-cpu"]
    threads: cluster["run_concoct"]["cpus-per-task"]
    params:
       outFolder=directory("%s/CONCOCT/{fafilter}" % (config["project-folder"])),
       maxClusters=config["params"]["concoct"]["maxClusters"]
    singularity: config["singularity"]["concoct"]
    shell:"""
      concoct --composition_file {input.fa} --coverage_file {input.coverage} -c {params.maxClusters} -t {threads} -b {params.outFolder} &> {log}
    """
    
rule merge_cutup_concoct:
    """
    Merge cutup (CONCOCT).
    """
    input:
        "%s/CONCOCT/{fafilter}_clustering_gt1000.csv" % (config["project-folder"])
    output:
        "%s/CONCOCT/{fafilter}_clustering_merged.csv" % (config["project-folder"])
    log:
        "%s/logs/merge_cutup_concoct_{fafilter}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/merge_cutup_concoct_{fafilter}.tsv" % (config["project-folder"])
    singularity: config["singularity"]["concoct"]
    shell:"""
      merge_cutup_clustering.py {input} > {output} 2> {log}
    """

rule extract_fasta_concoct:
    """
    Extract fasta (CONCOCT).
    """
    input:
        fa="%s/MEGAHIT/final.contigs.fa" % (config["project-folder"]),
        csv="%s/CONCOCT/{fafilter}_clustering_merged.csv" % (config["project-folder"])
    output:
        directory("%s/CONCOCT/fasta_bins_{fafilter}" % (config["project-folder"]))
    log:
        "%s/logs/extract_fasta_concoct_{fafilter}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/extract_fasta_concoct_{fafilter}.tsv" % (config["project-folder"])
    params:
        out=directory("%s/CONCOCT/fasta_bins_{fafilter}" % (config["project-folder"]))
    singularity: config["singularity"]["concoct"]
    shell:"""
      mkdir -p {params.out}
      
      extract_fasta_bins.py {input.fa} {input.csv} --output_path {params.out}
    """