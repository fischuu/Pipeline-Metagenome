rule R_finalReport:
    """
    Create the final report (R).
    """
    input:
        script=config["report-script"],
    output:
        "%s/finalReport.html" % (config["project-folder"])
    log:
        "%s/logs/R/finalReport.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/R_finalReport.tsv" % (config["project-folder"])
    singularity:
        config["singularity"]["r-gbs"]
    params:
       projFolder=config["project-folder"],
       pipeFolder=config["pipeline-folder"],
       pipeConfig=config["pipeline-config"],
       serverConfig=config["server-config"],
    shell:"""
       R -e "projFolder <- '{params.projFolder}'; \
             pipelineFolder <- '{params.pipeFolder}'; \
             pipelineConfig.file <- '{params.pipeConfig}'; \
             serverConfig.file <- '{params.serverConfig}'; \
             snakemake <- TRUE;\
             options(knitr.duplicate.label = 'allow');\
             source('{input.script}')" &> {log}
    """