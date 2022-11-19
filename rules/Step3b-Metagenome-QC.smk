rule check_full_assembly:
    """
    Check full assembly (METAQUAST).
    """
    input:
        "%s/MEGAHIT/final.contigs.fa" % (config["project-folder"])
    output:
        "%s/METAQUAST/final.contigs.meta" % (config["project-folder"])
    log:
        "%s/logs/check_full_assembly.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/check_full_assembly.benchmark.tsv" % (config["project-folder"])
    singularity: config["singularity"]["metaquast"]
    resources:
        time=cluster["check_full_assembly"]["time"],
        mem=cluster["check_full_assembly"]["mem-per-cpu"]
    threads: cluster["check_full_assembly"]["cpus-per-task"]
    shell:"""
        metaquast.py -t {threads} --no-plots -o {output} -1 {input}
    """

rule check_group_assembly:
    """
    Check group assembly (METAQUAST).
    """
    input:
        "%s/MEGAHIT/final.contigs.group_{cagroup}.fa" % (config["project-folder"])
    output:
        "%s/METAQUAST/final.contigs.group_{cagroup}.meta" % (config["project-folder"])
    log:
        "%s/logs/check_full_assembly.group_{cagroup}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/check_full_assembly.group_{cagroup}.benchmark.tsv" % (config["project-folder"])
    singularity: config["singularity"]["metaquast"]
    resources:
        time=cluster["check_full_assembly"]["time"],
        mem=cluster["check_full_assembly"]["mem-per-cpu"]
    threads: cluster["check_full_assembly"]["cpus-per-task"]
    shell:"""
        metaquast.py -t {threads} --no-plots -o {output} -1 {input}
    """