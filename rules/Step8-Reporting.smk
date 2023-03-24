rule R_finalReport:
    """
    Create the final report (R).
    """
    input:
        "%s/QC/RAW/multiqc_R1/" % (config["project-folder"]),
        "%s/QC/CONCATENATED/multiqc_R1/" % (config["project-folder"]),
        "%s/QC/TRIMMED/multiqc_R1/" % (config["project-folder"]),
        "%s/QC/DECONTAMINATED/multiqc_R1/" % (config["project-folder"]),
        expand("%s/QC/RAW/{rawsamples}_R1_qualdist.txt" % (config["project-folder"]), rawsamples=rawsamples),
        expand("%s/QC/CONCATENATED/{samples}_R1_qualdist.txt" % (config["project-folder"]), samples=samples),
        expand("%s/QC/TRIMMED/{samples}_R1_qualdist.txt" % (config["project-folder"]), samples=samples),
        expand("%s/QC/DECONTAMINATED/{samples}_R1_qualdist.txt" % (config["project-folder"]), samples=samples),
        script=config["report-script"]
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