rule R_createSummaries_file:
    """
    Create the file summaries.tsv (R).
    """
    input:
        script=config["summary-script"],
        samplesheet=config["samplesheet-file"]
    output:
         "%s/RESULTS/summary.tsv" % (config["project-folder"])
    log:
        "%s/logs/R_createSummaries_file.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/R_createSummaries.tsv" % (config["project-folder"])
    singularity:
        config["singularity"]["r-gbs"]
    threads: cluster["R_createSummaries_file"]["cpus-per-task"]
    resources:
        time=cluster["R_createSummaries_file"]["time"],
        mem=cluster["R_createSummaries_file"]["mem-per-cpu"]
    params:
       projFolder=config["project-folder"],
       pipeFolder=config["pipeline-folder"],
       pipeConfig=config["pipeline-config"]
    shell:"""
       Rscript {input.script} {input.samplesheet} {output} &> {log}
    """