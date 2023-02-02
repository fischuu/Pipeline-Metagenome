rule check_assemblies:
    """
    Check the different assemblies (METAQUAST).
    """
    input:
        full="%s/MEGAHIT/final.contigs.fa" % (config["project-folder"]),
        coas=expand("%s/MEGAHIT/final.contigs.group_{cagroup}.fa" % (config["project-folder"]), cagroup=assemblyGroups),
        fullfilter1k="%s/MEGAHIT/final.contigs.1k.fa" % (config["project-folder"]),
        fullfilter2k="%s/MEGAHIT/final.contigs.2k.fa" % (config["project-folder"]),
        coasfilter1k=expand("%s/MEGAHIT/final.contigs.group_{cagroup}.1k.fa" % (config["project-folder"]), cagroup=assemblyGroups),
        coasfilter2k=expand("%s/MEGAHIT/final.contigs.group_{cagroup}.2k.fa" % (config["project-folder"]), cagroup=assemblyGroups)
    output:
        "%s/QUAST/report.html" % (config["project-folder"])
    log:
        "%s/logs/check_assemblies.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/check_assemblies.benchmark.tsv" % (config["project-folder"])
    singularity: config["singularity"]["quast"]
    resources:
        time=cluster["check_assemblies"]["time"],
        mem=cluster["check_assemblies"]["mem-per-cpu"]
    threads: cluster["check_assemblies"]["cpus-per-task"]
    params: out="%s/QUAST" % (config["project-folder"])
    shell:"""
        metaquast.py \
            -o {params.out} \
	          {input.full} {input.coas} {input.fullfilter1k} {input.fullfilter2k} {input.coasfilter1k} {input.coasfilter2k}\
            --max-ref-number 0 \
            --threads {threads}
    """