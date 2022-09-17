rule cut_filtered_contigs_concoct:
    """
    Cut contigs into smaller contigs (CONCOCT).
    """
    input:
        con1k="%s/MEGAHIT/final.contigs.1k.fa" % (config["project-folder"]),
        con2k="%s/MEGAHIT/final.contigs.1k.fa" % (config["project-folder"])
    output:
        bed1k="%s/CONCOCT/final.contigs.1k_10K.bed" % (config["project-folder"]),
        fa1k="%s/CONCOCT/final.contigs.1k_10K.fa" % (config["project-folder"]),
        bed2k="%s/CONCOCT/final.contigs.2k_10K.bed" % (config["project-folder"]),
        fa2k="%s/CONCOCT/final.contigs.2k_10K.fa" % (config["project-folder"])
    log:
        "%s/logs/cut_filtered_contigs_concoct.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/cut_filtered_contigs_concoct.tsv" % (config["project-folder"])
    resources:
        time=cluster["cut_filtered_contigs_concoct"]["time"],
        mem=cluster["cut_filtered_contigs_concoct"]["mem-per-cpu"]
    threads: cluster["cut_filtered_contigs_concoct"]["cpus-per-task"]
    singularity: config["singularity"]["concoct"]
    shell:"""
        cut_up_fasta.py {input.con1k} -c 10000 -o 0 --merge_last -b {output.bed1k} > {output.fa1k} 2> {log}
        cut_up_fasta.py {input.con2k} -c 10000 -o 0 --merge_last -b {output.bed2k} > {output.fa2k} 2> {log}
    """
    
checkpoint cut_input_bed_1k:
    """
    Divide the input bed-file for parallel processing
    """
    input:
        "%s/CONCOCT/final.contigs.1k_10K.bed" % (config["project-folder"]),
    output:
        directory("%s/CONCOCT/bed_subsampled_1k" % (config["project-folder"])),
    params:
        out="%s/CONCOCT/bed_subsampled_1k/final.contigs.1k_10K." % (config["project-folder"]),
    shell:"""
        mkdir -p {output}

        split -l 10000000 --numeric-suffixes {input} {params.out}
    """
    
checkpoint cut_input_bed_2k:
    """
    Divide the input bed-file for parallel processing
    """
    input:
        "%s/CONCOCT/final.contigs.2k_10K.bed" % (config["project-folder"]),
    output:
        directory("%s/CONCOCT/bed_subsampled_2k" % (config["project-folder"])),
    params:
        out="%s/CONCOCT/bed_subsampled_2k/final.contigs.2k_10K." % (config["project-folder"]),
    shell:"""
        mkdir -p {output}

        split -l 10000000 --numeric-suffixes {input} {params.out}
    """
    
rule process_covTable_bed_parallel_1k:
    """
    Calculate the coverage for the subpart of the splitted bed
    """
    input:
        bed="%s/CONCOCT/bed_subsampled_1k/final.contigs.1k_10K.{i}" % (config["project-folder"]),
        bam=expand("%s/BAM/megahit/{samples}_mega.bam" % (config["project-folder"]), samples=samples)    
    output:
        woheader="%s/CONCOCT/coverage_table_1k_{i}.tsv" % (config["project-folder"]),
        header="%s/CONCOCT/coverage_table_1k_header_{i}.tsv" % (config["project-folder"]),
    params:
       bam="%s/BAM/megahit/*.bam" % (config["project-folder"])
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

     
rule process_covTable_bed_parallel_2k:
    """
    Calculate the coverage for the subpart of the splitted bed
    """
    input:
        bed="%s/CONCOCT/bed_subsampled_2k/final.contigs.2k_10K.{i}" % (config["project-folder"]),
        bam=expand("%s/BAM/megahit/{samples}_mega.bam" % (config["project-folder"]), samples=samples)    
    output:
        woheader="%s/CONCOCT/coverage_table_2k_{i}.tsv" % (config["project-folder"]),
        header="%s/CONCOCT/coverage_table_2k_header_{i}.tsv" % (config["project-folder"]),
    params:
       bam="%s/BAM/megahit/*.bam" % (config["project-folder"])
    resources:
        time=cluster["process_covTable_bed_parallel"]["time"],
        mem=cluster["process_covTable_bed_parallel"]["mem-per-cpu"]
    threads: cluster["process_covTable_bed_parallel"]["cpus-per-task"]
    singularity: config["singularity"]["concoct"]
    shell:"""
      concoct_coverage_table.py {input.bed} {params.bam} > {output.woheader}
      
      # Extract the header  
      sed -n 1p {output.woheader} > {output.header}
      
      # Remove the header
      sed -i '1d' {output.woheader}
     """
     
def aggregate_input_1k(wildcards):
    """
    Aggregate the input object for the final coverage table
    """
    checkpoint_output = checkpoints.cut_input_bed_1k.get(**wildcards).output[0]
    return expand("%s/CONCOCT/coverage_table_1k_{i}.tsv" % (config["project-folder"]),
                  i=glob_wildcards(os.path.join(checkpoint_output, "final.contigs.1k_10K.{i}")).i)      
      
def aggregate_input_2k(wildcards):
    """
    Aggregate the input object for the final coverage table
    """
    checkpoint_output = checkpoints.cut_input_bed_2k.get(**wildcards).output[0]
    return expand("%s/CONCOCT/coverage_table_2k_{i}.tsv" % (config["project-folder"]),
                  i=glob_wildcards(os.path.join(checkpoint_output, "final.contigs.2k_10K.{i}")).i)      
      
rule aggregate_coverage_table_1k:
    input:
        aggregate_input_1k
    output:
        all="%s/CONCOCT/coverage_table_1k.tsv" % (config["project-folder"]),
        tmp="%s/CONCOCT/coverage_table_1k.tsv.tmp" % (config["project-folder"]),
        tmp2="%s/CONCOCT/coverage_table_1k.tsv.tmp2" % (config["project-folder"])
    params:
        header="%s/CONCOCT/coverage_table_1k_header_00.tsv" % (config["project-folder"])
    shell:"""
        # This was my manual approach to it... Keep this as reminder until below approach is validated in larger test run
        # head -n 1 coverage_table_00.tsv > header
        # tail -n +2 "coverage_table_00.tsv" > ct_00.tsv
        # tail -n +2 "coverage_table_01.tsv" > ct_01.tsv
        # tail -n +2 "coverage_table_02.tsv" > ct_02.tsv
        # tail -n +2 "coverage_table_03.tsv" > ct_03.tsv
        # tail -n +2 "coverage_table_04.tsv" > ct_04.tsv
        # tail -n +2 "coverage_table_05.tsv" > ct_05.tsv
        # tail -n +2 "coverage_table_06.tsv" > ct_06.tsv
        # tail -n +2 "coverage_table_07.tsv" > ct_07.tsv
        
        # cat ct_00.tsv ct_01.tsv ct_02.tsv ct_03.tsv ct_04.tsv ct_05.tsv ct_06.tsv ct_07.tsv > ct.tsv
        # sort -t_ -k2,2n ct.tsv > ct_sorted.tsv
        # cat header ct_sorted.tsv > coverage_table.tsv
        # rm -rf ct*

        cat {input} > {output.tmp}
        sort -t_ -k2,2n {output.tmp} > {output.tmp2}
        cat {params.header} {output.tmp2} > {output.all}
    """
    
rule aggregate_coverage_table_2k:
    input:
        aggregate_input_2k
    output:
        all="%s/CONCOCT/coverage_table_2k.tsv" % (config["project-folder"]),
        tmp="%s/CONCOCT/coverage_table_2k.tsv.tmp" % (config["project-folder"]),
        tmp2="%s/CONCOCT/coverage_table_2k.tsv.tmp2" % (config["project-folder"])
    params:
        header="%s/CONCOCT/coverage_table_2k_header_00.tsv" % (config["project-folder"])
    shell:"""
        # This was my manual approach to it... Keep this as reminder until below approach is validated in larger test run
        # head -n 1 coverage_table_00.tsv > header
        # tail -n +2 "coverage_table_00.tsv" > ct_00.tsv
        # tail -n +2 "coverage_table_01.tsv" > ct_01.tsv
        # tail -n +2 "coverage_table_02.tsv" > ct_02.tsv
        # tail -n +2 "coverage_table_03.tsv" > ct_03.tsv
        # tail -n +2 "coverage_table_04.tsv" > ct_04.tsv
        # tail -n +2 "coverage_table_05.tsv" > ct_05.tsv
        # tail -n +2 "coverage_table_06.tsv" > ct_06.tsv
        # tail -n +2 "coverage_table_07.tsv" > ct_07.tsv
        
        # cat ct_00.tsv ct_01.tsv ct_02.tsv ct_03.tsv ct_04.tsv ct_05.tsv ct_06.tsv ct_07.tsv > ct.tsv
        # sort -t_ -k2,2n ct.tsv > ct_sorted.tsv
        # cat header ct_sorted.tsv > coverage_table.tsv
        # rm -rf ct*
        
        cat {input} > {output.tmp}
        sort -t_ -k2,2n {output.tmp} > {output.tmp2}
        cat {params.header} {output.tmp2} > {output.all}
    """
    
rule run_concoct_1k:
    """
    Apply concoct to the data (CONCOCT).
    """
    input:
        fa="%s/CONCOCT/final.contigs.1k_10K.fa" % (config["project-folder"]),
        coverage="%s/CONCOCT/coverage_table_1k.tsv" % (config["project-folder"])
    output:
        "%s/CONCOCT/1k_clustering_gt1000.csv" % (config["project-folder"])
    log:
        "%s/logs/run_concoct_1k.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/run_concoct_1k.tsv" % (config["project-folder"])
    resources:
        time=cluster["run_concoct_1k"]["time"],
        mem=cluster["run_concoct_1k"]["mem-per-cpu"]
    threads: cluster["run_concoct_1k"]["cpus-per-task"]
    params:
       outFolder=directory("%s/CONCOCT/1k" % (config["project-folder"])),
       maxClusters=config["params"]["concoct"]["maxClusters"]
    singularity: config["singularity"]["concoct"]
    shell:"""
      concoct --composition_file {input.fa} --coverage_file {input.coverage} -c {params.maxClusters} -t {threads} -b {params.outFolder} &> {log}
    """
    
rule run_concoct_2k:
    """
    Apply concoct to the data (CONCOCT).
    """
    input:
        fa="%s/CONCOCT/final.contigs.2k_10K.fa" % (config["project-folder"]),
        coverage="%s/CONCOCT/coverage_table_2k.tsv" % (config["project-folder"])
    output:
        "%s/CONCOCT/2k_clustering_gt1000.csv" % (config["project-folder"])
    log:
        "%s/logs/run_concoct_2k.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/run_concoct_2k.tsv" % (config["project-folder"])
    resources:
        time=cluster["run_concoct_2k"]["time"],
        mem=cluster["run_concoct_2k"]["mem-per-cpu"]
    threads: cluster["run_concoct_2k"]["cpus-per-task"]
    params:
       outFolder=directory("%s/CONCOCT/2k" % (config["project-folder"])),
       maxClusters=config["params"]["concoct"]["maxClusters"]
    singularity: config["singularity"]["concoct"]
    shell:"""
      concoct --composition_file {input.fa} --coverage_file {input.coverage} -c {params.maxClusters} -t {threads} -b {params.outFolder} &> {log}
    """
    
rule merge_cutup_concoct_1k:
    """
    Merge cutup (CONCOCT).
    """
    input:
        "%s/CONCOCT/1k_clustering_gt1000.csv" % (config["project-folder"])
    output:
        "%s/CONCOCT/1k_clustering_merged.csv" % (config["project-folder"])
    log:
        "%s/logs/merge_cutup_concoct_1k.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/merge_cutup_concoct_1k.tsv" % (config["project-folder"])
    singularity: config["singularity"]["concoct"]
    shell:"""
      merge_cutup_clustering.py {input} > {output} 2> {log}
    """
    
rule merge_cutup_concoct_2k:
    """
    Merge cutup (CONCOCT).
    """
    input:
        "%s/CONCOCT/2k_clustering_gt1000.csv" % (config["project-folder"])
    output:
        "%s/CONCOCT/2k_clustering_merged.csv" % (config["project-folder"])
    log:
        "%s/logs/merge_cutup_concoct_2k.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/merge_cutup_concoct_2k.tsv" % (config["project-folder"])
    singularity: config["singularity"]["concoct"]
    shell:"""
      merge_cutup_clustering.py {input} > {output} 2> {log}
    """
    
rule extract_fasta_concoct_1k:
    """
    Extract fasta (CONCOCT).
    """
    input:
        fa="%s/MEGAHIT/final.contigs.fa" % (config["project-folder"]),
        csv="%s/CONCOCT/1k_clustering_merged.csv" % (config["project-folder"])
    output:
        directory("%s/CONCOCT/fasta_bins_1k" % (config["project-folder"]))
    log:
        "%s/logs/extract_fasta_concoct_1k.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/extract_fasta_concoct_1k.tsv" % (config["project-folder"])
    params:
        out=directory("%s/CONCOCT/fasta_bins_1k" % (config["project-folder"]))
    singularity: config["singularity"]["concoct"]
    shell:"""
      mkdir -p {params.out}
      
      extract_fasta_bins.py {input.fa} {input.csv} --output_path {params.out}
    """
    
    
rule extract_fasta_concoct_2k:
    """
    Extract fasta (CONCOCT).
    """
    input:
        fa="%s/MEGAHIT/final.contigs.fa" % (config["project-folder"]),
        csv="%s/CONCOCT/2k_clustering_merged.csv" % (config["project-folder"])
    output:
        directory("%s/CONCOCT/fasta_bins_2k" % (config["project-folder"]))
    log:
        "%s/logs/extract_fasta_concoct_2k.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/extract_fasta_concoct_2k.tsv" % (config["project-folder"])
    params:
        out=directory("%s/CONCOCT/fasta_bins_2k" % (config["project-folder"]))
    singularity: config["singularity"]["concoct"]
    shell:"""
      mkdir -p {params.out}
      
      extract_fasta_bins.py {input.fa} {input.csv} --output_path {params.out}
    """