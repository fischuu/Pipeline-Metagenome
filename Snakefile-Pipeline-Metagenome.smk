import pandas as pd
from snakemake.utils import validate, min_version
from multiprocessing import cpu_count
import glob
import re
import os, sys
import yaml

from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()
shell.executable("bash")

##### Metagenome Snakemake Pipeline              #####
##### Daniel Fischer (daniel.fischer@luke.fi)    #####
##### Natural Resources Institute Finland (Luke) #####

##### Version: 0.1.15
version = "0.1.15"

##### set minimum snakemake version #####
min_version("6.0")

##### load config and sample sheets #####

samplesheet = pd.read_table(config["samplesheet"]).set_index("rawsample", drop=False)
rawsamples=list(samplesheet.rawsample)
samples=list(set(list(samplesheet.sample_name)))  # Unique list
lane=list(samplesheet.lane)

workdir: config["project-folder"]

##### Complete the input configuration
config["report-script"] = config["pipeline-folder"]+"/scripts/workflow-report.Rmd"

wildcard_constraints:
    rawsamples="|".join(rawsamples),
    samples="|".join(samples)

##### Extract the cluster resource requests from the server config #####
cluster=dict()
if os.path.exists(config["server-config"]):
    with open(config["server-config"]) as yml:
        cluster = yaml.load(yml, Loader=yaml.FullLoader)

##### input checks #####

# rawdata-folder needs to end with "/", add it if missing:
if config["rawdata-folder"][-1] != '/':
   config["rawdata-folder"]=config["rawdata-folder"]+'/'

# project-folder should not end with "/", so remove it
if config["project-folder"][-1] == '/':
   config["project-folder"]=config["project-folder"][:-1]

if config["contamination-folder"] != "":
  if config["contamination-folder"][-1] == '/':
     config["contamination-folder"]=config["contamination-folder"][:-1]

##### input function definitions ######

def get_fastq_for_concatenating_read1(wildcards):
    r1 = samplesheet.loc[samplesheet["sample_name"] == wildcards.samples]["read1"]
    path = config["rawdata-folder"]
    output = [path + x for x in r1]
    return output   

def get_fastq_for_concatenating_read2(wildcards):
    r1 = samplesheet.loc[samplesheet["sample_name"] == wildcards.samples]["read2"]
    path = config["rawdata-folder"]
    output = [path + x for x in r1]
    return output  

def merge_r1_reads(wildcards):
    samples = set(samplesheet["sample_name"])
    if config["contamination-folder"] == "":
      out = expand("%s/FASTQ/TRIMMED/{samples}_R1.trimmed.fastq.gz" % (config["project-folder"]), samples=samples)
    else:
      out = expand("%s/FASTQ/DECONTAMINATED/{samples}_R1.decontaminated.fastq.gz" % (config["project-folder"]), samples=samples)
    return out

def merge_r2_reads(wildcards):
    samples = set(samplesheet["sample_name"])
    if config["contamination-folder"] == "":
      out = expand("%s/FASTQ/TRIMMED/{samples}_R2.trimmed.fastq.gz" % (config["project-folder"]), samples=samples)
    else:
      out = expand("%s/FASTQ/DECONTAMINATED/{samples}_R2.decontaminated.fastq.gz" % (config["project-folder"]), samples=samples)
    return out

##### Define the required docker images #####
config["singularity"] = {}
config["singularity"]["samtools"] = "docker://fischuu/samtools:1.9-0.2"
config["singularity"]["fastp"] = "docker://fischuu/fastp:0.20.1"
config["singularity"]["star"] = "docker://fischuu/star:2.7.5a"
config["singularity"]["wgs"] = "docker://fischuu/wgs:0.2"
config["singularity"]["r-gbs"] = "docker://fischuu/r-gbs:4.2.1-0.1"
config["singularity"]["megahit"] = "docker://fischuu/megahit:1.2.9-0.2"
config["singularity"]["bowtie2"] = "docker://fischuu/bowtie2:2.4.4-0.1"
config["singularity"]["prodigal"] = "docker://fischuu/prodigal:2.6.3-0.1"
config["singularity"]["eggnog"] = "docker://fischuu/eggnog:latest"
config["singularity"]["subread"] = "docker://fischuu/subread:2.0.1-0.1"
config["singularity"]["concoct"] = "docker://nanozoo/concoct:latest"

##### Deriving runtime paramteres ######

##### Print some welcoming summary #####

##### Print the welcome screen #####
print("#################################################################################")
print("##### Welcome to the Metagenome pipeline")
print("##### version: "+version)
print("##### Number of rawsamples : "+str(len(rawsamples)))
print("##### Number of samples    : "+str(len(samples)))
print("##### Rawdata folder       : "+config["rawdata-folder"])
print("##### Project folder       : "+config["project-folder"])
print("##### Contamination folder : " +config["contamination-folder"])
#print("##### Contamination refs   :" +config["contamination-refs"])
print("#####")
print("##### Singularity configuration")
print("##### --------------------------------")
print("##### concoct        : "+config["singularity"]["concoct"])
print("##### samtools       : "+config["singularity"]["samtools"])
print("##### fastp          : "+config["singularity"]["fastp"])
print("##### star           : "+config["singularity"]["star"])
print("##### wgs            : "+config["singularity"]["wgs"])
print("##### r-gbs          : "+config["singularity"]["r-gbs"])
print("##### megahit        : "+config["singularity"]["megahit"])
print("##### bowtie2        : "+config["singularity"]["bowtie2"])
print("##### prodigal       : "+config["singularity"]["prodigal"])
print("##### eggnog         : "+config["singularity"]["eggnog"])
print("##### subread        : "+config["singularity"]["subread"])

##### run complete pipeline #####

rule all:
    input:
      "%s/FASTQ/MERGED/all_merged_R1.fastq.gz" % (config["project-folder"]),
      "%s/FASTQ/MERGED/all_merged_R2.fastq.gz" % (config["project-folder"]),
      "%s/MEGAHIT/final.contigs.fa" % (config["project-folder"]),
      expand("%s/BAM/megahit/{samples}_mega.bam" % (config["project-folder"]), samples=samples),
      "%s/PRODIGAL/final.contigs.prodigal.gtf" % (config["project-folder"]),
      "%s/PRODIGAL/EGGNOG-DATA/" % (config["project-folder"]),
      "%s/PRODIGAL/eggnog.emapper.seed_orthologs" % (config["project-folder"]),
      "%s/EGGNOG/eggnog_output.emapper.annotations" % (config["project-folder"]),
      expand("%s/QUANTIFICATION/PRODIGAL_FC/{samples}_fc.txt" % (config["project-folder"]), samples=samples),
      "%s/CONCOCT/final.contigs.1k_10K.fa" % (config["project-folder"]),
      "%s/CONCOCT/final.contigs.2k_10K.fa" % (config["project-folder"]),
      "%s/CONCOCT/coverage_table_1k.tsv" % (config["project-folder"]),
      "%s/CONCOCT/coverage_table_2k.tsv" % (config["project-folder"]),
      "%s/CONCOCT/1k_clustering_gt1000.csv" % (config["project-folder"]),
      "%s/CONCOCT/2k_clustering_gt1000.csv" % (config["project-folder"]),
      "%s/CONCOCT/fasta_bins_1k" % (config["project-folder"]),
      "%s/CONCOCT/fasta_bins_2k" % (config["project-folder"])

rule preparations:
    input:

rule trimming:
    input:
      expand("%s/BAM/{samples}.bam" % (config["project-folder"]), samples=samples)
      
rule qc:
    input:
      expand("%s/QC/RAW/{rawsamples}_R1_001_fastqc.zip" % (config["project-folder"]), rawsamples=rawsamples),
      expand("%s/QC/RAW/{rawsamples}_R2_001_fastqc.zip" % (config["project-folder"]), rawsamples=rawsamples),
      expand("%s/QC/TRIMMED/{rawsamples}_R1_fastqc.zip" % (config["project-folder"]), rawsamples=rawsamples),
      expand("%s/QC/TRIMMED/{rawsamples}_R2_fastqc.zip" % (config["project-folder"]), rawsamples=rawsamples)

rule alignment:
    input:

### setup report #####
report: "report/workflow.rst"

##### load rules #####
include: "rules/Step1-Preparations.smk"
include: "rules/Step2-ReadProcessing.smk"
include: "rules/Step2b-Decontamination.smk"
include: "rules/Step3-CreateMetagenome.smk"
include: "rules/Step4-GenePrediction.smk"
include: "rules/Step5-MAGs.smk"
