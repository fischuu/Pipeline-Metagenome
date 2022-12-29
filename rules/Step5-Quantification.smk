rule quantify_full_predictedGenes_featureCounts:
    """
    Quantify the mapped reads after merging (featureCounts).
    """
    input:
        bam="%s/BAM/final.contigs_full/{samples}_mega.bam" % (config["project-folder"]),
        gtf="%s/PRODIGAL/final.contigs_full/final.contigs_full.prodigal.gtf" % (config["project-folder"])
    output:
        file="%s/QUANTIFICATION/PRODIGAL_FC/final.contigs_full/{samples}_full_fc.txt" % (config["project-folder"])
    log:
        "%s/logs/quantify_full_predictedGenes_featureCounts.{samples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/quantify_full_predictedGenes_featureCounts.{samples}.benchmark.tsv" % (config["project-folder"])
    resources:
        time=cluster["quantify_full_predictedGenes_featureCounts"]["time"],
        mem=cluster["quantify_full_predictedGenes_featureCounts"]["mem-per-cpu"]
    threads: cluster["quantify_full_predictedGenes_featureCounts"]["cpus-per-task"]
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
    
rule quantify_coas_predictedGenes_featureCounts:
    """
    Quantify the mapped reads after merging (featureCounts).
    """
    input:
        bam="%s/BAM/final.contigs_coas{cagroup}/{samples}_mega.bam" % (config["project-folder"]),
        gtf="%s/PRODIGAL/final.contigs_group_{cagroup}/final.contigs_group_{cagroup}.prodigal.gtf" % (config["project-folder"])
    output:
        file="%s/QUANTIFICATION/PRODIGAL_FC/final.contigs_group_{cagroup}/{samples}_group_{cagroup}_fc.txt" % (config["project-folder"])
    log:
        "%s/logs/quantify_group_{cagroup}_predictedGenes_featureCounts.{samples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/quantify_group_{cagroup}_predictedGenes_featureCounts.{samples}.benchmark.tsv" % (config["project-folder"])
    resources:
        time=cluster["quantify_coas_predictedGenes_featureCounts"]["time"],
        mem=cluster["quantify_coas_predictedGenes_featureCounts"]["mem-per-cpu"]
    threads: cluster["quantify_coas_predictedGenes_featureCounts"]["cpus-per-task"]
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
  
rule Full_AlignmentStats:
    """
    Get the mapping stats for the full assembly alignments.
    """
    input:
        "%s/BAM/final.contigs_full/{samples}_mega.bam" % (config["project-folder"])
    output:
        fs="%s/BAM/final.contigs_full/{samples}_mega.bam.flagstat" % (config["project-folder"]),
        cv="%s/BAM/final.contigs_full/{samples}_mega.bam.coverage" % (config["project-folder"]),
    log:
        "%s/logs/Full_AlignmentStats.{samples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/Full_AlignmentStats.{samples}.tsv" % (config["project-folder"])
    threads: cluster["Full_AlignmentStats"]["cpus-per-task"]
    resources:
        time=cluster["Full_AlignmentStats"]["time"],
        mem=cluster["Full_AlignmentStats"]["mem-per-cpu"]
    singularity: config["singularity"]["samtools"]
    shell:"""
        samtools flagstat {input} > {output.fs}
        samtools idxstats {input} | awk '{{print $1\" \"$3}}' > {output.cv}
  	"""  
  	
rule Coas_AlignmentStats:
    """
    Get the mapping stats for the different co-assemblies.
    """
    input:
        "%s/BAM/final.contigs_coas{cagroup}/{samples}_mega.bam" % (config["project-folder"])
    output:
        fs="%s/BAM/final.contigs_coas{cagroup}/{samples}_mega.bam.flagstat" % (config["project-folder"]),
        cv="%s/BAM/final.contigs_coas{cagroup}/{samples}_mega.bam.coverage" % (config["project-folder"]),
    log:
        "%s/logs/coas{cagroup}_AlignmentStats.{samples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/coas{cagroup}_AlignmentStats.{samples}.tsv" % (config["project-folder"])
    threads: cluster["Coas_AlignmentStats"]["cpus-per-task"]
    resources:
        time=cluster["Coas_AlignmentStats"]["time"],
        mem=cluster["Coas_AlignmentStats"]["mem-per-cpu"]
    singularity: config["singularity"]["samtools"]
    shell:"""
        samtools flagstat {input} > {output.fs}
        samtools idxstats {input} | awk '{{print $1\" \"$3}}' > {output.cv}
  	"""  