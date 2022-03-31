# vim: set filetype=sh :
import tempfile
import pandas as pd
from snakemake.utils import validate, min_version
##### set minimum snakemake version #####
min_version("6.0")

##### Snakemake Metagenomics pipeline #####
##### Daniel Fischer (daniel.fischer@luke.fi)
##### Natural Resources Institute Finland (Luke)
##### Version: 0.1
version = "0.1"


##### sample sheets #####

# Put those later into an own config file
config["project-folder"] = "/home/ejo138/git/Pipeline-Metagenome/"
config["samplesheet"] = "toyData.tsv"
config["pipeline-folder"] = "/home/ejo138/git/Pipeline-Metagenome/"

# These are the derived config values
config["samplesheet-full"] = config["project-folder"] + config["samplesheet"]

### Import the data ###
samples = pd.read_table(config["samplesheet-full"], delimiter='\s+', lineterminator='\n').set_index("sample", drop=False)

### setup report #####

report: "report/workflow.rst"

### Define functions ###
def get_fq(wildcards):
    sample_name = wildcards.sample
    r1= config["project-folder"]+"/FASTQ/DECONTAMINATED/" + sample_name + "_R1.decontaminated.fastq.gz"
    r2= config["project-folder"]+"/FASTQ/DECONTAMINATED/" + sample_name + "_R2.decontaminated.fastq.gz"
    result=[r1, r2]
    return result


##### Print the welcome screen #####
print("#################################################################################")
print("##### Welcome to the Snakemake Metagenome pipeline")
print("##### Version: "+version)
print("#####")
print("##### Pipeline configuration")
print("##### --------------------------------")
print("##### project-folder    : "+config["project-folder"])
print("##### pipeline-folder   : "+config["pipeline-folder"])
print("##### samplesheet       : "+config["samplesheet-full"])
print("##### Number of samples : "+str(samples.shape[0]))

##### run complete pipeline #####

rule all:
    input:
        expand("test.{sample}.txt", sample=samples.index)


rule a:
    output:
        "test.{sample}.txt"
    shell:
        "touch {output}"

#rule all:
#    input:
#     # PRODIGAL OUTPUT
#      "%s/prodigal_out/final.contigs.prodigal.fa" % (config["project-folder"]),
#      "%s/prodigal_out/final.contigs.prodigal.chunk.txt" % (config["project-folder"]),
#      expand("%s/prodigal_out/{prodigalChunks}.emapper.seed_orthologs" % (config["project-folder"]), prodigalChunks=prodigalChunks),
#      "%s/prodigal_out/input_file.emapper.seed_orthologs" % (config["project-folder"]),
#      "%s/prodigal_out/eggnog_output.emapper.annotations" % (config["project-folder"]),
#      expand("%s/BAM/megahit/{samples}_mega.bam" % (config["project-folder"]), samples=samples),
#      "%s/megahit_out/final.contigs.1k.fa" % (config["project-folder"]),
#      "%s/megahit_out/final.contigs.2k.fa" % (config["project-folder"]),
#      expand("%s/QUANTIFICATION/prodigal_fc/{samples}_fc.txt" % (config["project-folder"]), samples=samples),
#      "%s/concoct_out/coverage_table.tsv" % (config["project-folder"]),
#      "%s/concoct_out/clustering_gt1000.csv" % (config["project-folder"]),
#      "%s/concoct_out/clustering_merged.csv" % (config["project-folder"]),
#      "%s/concoct_out/fasta_bins" % (config["project-folder"]),
#      "%s/checkm_out/lineage.ms" % (config["project-folder"]),
# Those inputs are temporry and can be deleted later
#        "%s/concoct_out/final.contigs.1k_10K.bed" % (config["project-folder"]),
#        "%s/concoct_out/final.contigs.1k_10K.fa" % (config["project-folder"]),
#        "%s/concoct_out/final.contigs.2k_10K.bed" % (config["project-folder"]),
#        "%s/concoct_out/final.contigs.2k_10K.fa" % (config["project-folder"]),
#        "%s/concoct_out/coverage_table.tsv" % (config["project-folder"]),
#        "%s/concoct_out/coverage_table_1k.tsv" % (config["project-folder"]),
#        "%s/concoct_out/coverage_table_2k.tsv" % (config["project-folder"]),
#        "%s/concoct_out/clustering_merged_1k.csv" % (config["project-folder"]),
#        "%s/concoct_out/clustering_merged_2k.csv" % (config["project-folder"]),
#        "%s/concoct_out/fasta_bins_1k" % (config["project-folder"]),
#        "%s/concoct_out/fasta_bins_2k" % (config["project-folder"]),
#        "%s/checkm_out/1k/lineage.ms" % (config["project-folder"]),
#        "%s/checkm_out/2k/lineage.ms" % (config["project-folder"])
#        

##### load rules #####

#include: "rules/00-control_quality_fastqc.smk"
#include: "rules/01-trim_inputfiles_trimmomatic.smk"
#include: "rules/02-concatenate_inputfiles_cat.smk"
#include: "rules/03-assemble_contigs_megahit.smk"
#include: "rules/12-prodigal_gene_prediction.smk"
#include: "rules/04-create_index_bowtie2.smk"
#include: "rules/05-map_data_metagenome_bowtie2.smk"
#include: "rules/18-quantify_predictedGenes_featureCounts.smk"
#include: "rules/13-cut_prodigal_bash.smk"
#include: "rules/15-copy_database_to_scratch.smk"
#include: "rules/14-eggnog_find_homology.smk"
#include: "rules/16-eggnog_concatenate_homologs_cat.smk"
#include: "rules/17-eggnog_orthology.smk"
#include: "rules/06-cut_contigs_concoct.smk"
#include: "rules/07-coverage_table_concoct.smk"
#include: "rules/08-run_concoct.smk"
#include: "rules/09-merge_cutup_concoct.smk"
#include: "rules/10-extract_fasta_concoct.smk"
#include: "rules/11-qc_checkm.smk"
