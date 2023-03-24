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

##### Metagenomics Snakemake Pipeline              #####
##### Daniel Fischer (daniel.fischer@luke.fi)    #####
##### Natural Resources Institute Finland (Luke) #####

##### Version: 0.4
version = "0.4"

##### set minimum snakemake version #####
min_version("6.0")


##### Fill the configuration lines for relative paths
# project-folder should not end with "/", so remove it

if config["project-folder"][-1] == '/':
   config["project-folder"]=config["project-folder"][:-1]
   
if(config["pipeline-config"][0]!='/'):
    config["pipeline-config"] = config["project-folder"] + '/' + config["pipeline-config"]

if(config["server-config"][0]!='/'):
    config["server-config"] = config["project-folder"] + '/' + config["server-config"]

if(config["rawdata-folder"][0]!='/'):
    config["rawdata-folder"] = config["project-folder"] + '/' + config["rawdata-folder"]

if(config["contamination-folder"][0]!='/'):
    config["contamination-folder"] = config["project-folder"] + '/' + config["contamination-folder"]
      
if(config["samplesheet-file"][0]!='/'):
    config["samplesheet-file"] = config["project-folder"] + '/' + config["samplesheet-file"]
    
if(config["tmp"][0]!='/'):
    config["tmp"] = config["project-folder"] + '/' + config["tmp"]

##### load config and sample sheets #####

samplesheet = pd.read_table(config["samplesheet-file"]).set_index("rawsample", drop=False)
rawsamples=list(samplesheet.rawsample)
samples=list(set(list(samplesheet.sample_name)))  # Unique list
lane=list(samplesheet.lane)

assemblyGroups=list(set(list(samplesheet.assemblyGroup)))  # Unique list
assemblyGroups = [str(i) for i in assemblyGroups]          # Needed if more than one group per sample is given (via ',')
assemblyGroups = [i.split(',') for i in assemblyGroups]
assemblyGroups = list(set(sum(assemblyGroups, [])))
assemblyGroups = [int(i) for i in assemblyGroups]
samplesheet["assemblyGroup"] = samplesheet["assemblyGroup"].astype(str)

fafilter = ["1k", "2k"]

workdir: config["project-folder"]

##### Complete the input configuration
config["report-script"] = config["pipeline-folder"]+"/scripts/workflow-report.Rmd"

wildcard_constraints:
    rawsamples="|".join(rawsamples),
    samples="|".join(samples),
    cagroup="|".join(str(assemblyGroups)),
    fafilter="|".join(fafilter)

##### Extract the cluster resource requests from the server config #####
cluster=dict()
if os.path.exists(config["server-config"]):
    with open(config["server-config"]) as yml:
        cluster = yaml.load(yml, Loader=yaml.FullLoader)

##### input checks #####

# rawdata-folder needs to end with "/", add it if missing:
if config["rawdata-folder"][-1] != '/':
   config["rawdata-folder"]=config["rawdata-folder"]+'/'

if config["contamination-folder"] != "":
  if config["contamination-folder"][-1] == '/':
     config["contamination-folder"]=config["contamination-folder"][:-1]

# Get the basename fastq inputs
possible_ext = [".fastq", ".fq.gz", ".fastq.gz", ".fasta", ".fa", ".fa.gz", ".fasta.gz"]
ext = ".null"

reads1_tmp = list(samplesheet.read1)
reads1_trim = []
for r in reads1_tmp:
    for e in possible_ext:
        if r.endswith(e):
            addThis = r[:-len(e)]
            reads1_trim += [addThis]
            ext=e

reads2_tmp = list(samplesheet.read2)
reads2_trim = []
for r in reads2_tmp:
    for e in possible_ext:
        if r.endswith(e):
            addThis = r[:-len(e)]
            reads2_trim += [addThis] 
            ext=e

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

def merge_r1_ca_reads(wildcards):
    samples = samplesheet[samplesheet["assemblyGroup"].str.contains(str(wildcards.cagroup))]["sample_name"].drop_duplicates()
    if config["contamination-folder"] == "":
      out = expand("%s/FASTQ/TRIMMED/{samples}_R1.trimmed.fastq.gz" % (config["project-folder"]), samples=samples)
    else:
      out = expand("%s/FASTQ/DECONTAMINATED/{samples}_R1.decontaminated.fastq.gz" % (config["project-folder"]), samples=samples)
    return out

def merge_r2_ca_reads(wildcards):
    samples = samplesheet[samplesheet["assemblyGroup"].str.contains(str(wildcards.cagroup))]["sample_name"].drop_duplicates()
    if config["contamination-folder"] == "":
      out = expand("%s/FASTQ/TRIMMED/{samples}_R2.trimmed.fastq.gz" % (config["project-folder"]), samples=samples)
    else:
      out = expand("%s/FASTQ/DECONTAMINATED/{samples}_R2.decontaminated.fastq.gz" % (config["project-folder"]), samples=samples)
    return out

def get_raw_input_read_bs1(wildcards):
    reads = wildcards.reads1 + ext
    output = config["rawdata-folder"] + reads
    return output

def get_raw_input_read_bs2(wildcards):
    reads = wildcards.reads2 + ext
    output = config["rawdata-folder"] + reads
    return output

def get_raw_input_read1(wildcards):
    reads = samplesheet.loc[wildcards.rawsamples][["read1"]]
    path = config["rawdata-folder"]
    output = [path + x for x in reads]
    return output

def get_raw_input_read2(wildcards):
    reads = samplesheet.loc[wildcards.rawsamples][["read2"]]
    path = config["rawdata-folder"]
    output = [path + x for x in reads]
    return output

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
config["singularity"]["gbs"] = "docker://fischuu/gbs:0.2"
config["singularity"]["quast"] = "docker://fischuu/quast:5.2.0-0.2"

##### Deriving runtime paramteres ######

##### Print some welcoming summary #####

##### Print the welcome screen #####
print("#################################################################################")
print("##### Welcome to the Snakebite - Metagenomics pipeline")
print("##### version: "+version)
print("##### Number of rawsamples : "+str(len(rawsamples)))
print("##### Number of samples    : "+str(len(samples)))
print("##### Rawdata folder       : "+config["rawdata-folder"])
print("##### Project folder       : "+config["project-folder"])
print("##### Contamination folder : " +config["contamination-folder"])
print("##### No. of ass. groups   : "+str(len(assemblyGroups)))
print("##### Assembly groups      : "+str(assemblyGroups))
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
print("##### gbs            : "+config["singularity"]["gbs"])
print("##### quast          : "+config["singularity"]["quast"])


##### run complete pipeline #####

rule all:
    input:
      "%s/FASTQ/MERGED/all_merged_R1.fastq.gz" % (config["project-folder"]),
      "%s/FASTQ/MERGED/all_merged_R2.fastq.gz" % (config["project-folder"]),
      "%s/MEGAHIT/final.contigs.fa" % (config["project-folder"]),
      "%s/QUAST/report.html" % (config["project-folder"]),
      expand("%s/BAM/final.contigs_full/{samples}_mega.bam" % (config["project-folder"]), samples=samples),
      expand("%s/BAM/final.contigs_coas{cagroup}/{samples}_mega.bam" % (config["project-folder"]), cagroup=assemblyGroups, samples=samples),
      "%s/PRODIGAL/final.contigs_full/final.contigs_full.prodigal.fa" % (config["project-folder"]),
      expand("%s/PRODIGAL/final.contigs_group_{cagroup}/final.contigs_group_{cagroup}.prodigal.fa" % (config["project-folder"]), cagroup=assemblyGroups),
      "%s/EGGNOG/final.contigs_full/eggnog_output.emapper.annotations" % (config["project-folder"]),
      expand("%s/EGGNOG/final.contigs_group_{cagroup}/eggnog_output.emapper.annotations" % (config["project-folder"]), cagroup=assemblyGroups),
      expand("%s/QUANTIFICATION/PRODIGAL_FC/final.contigs_full/{samples}_full_fc.txt" % (config["project-folder"]), samples=samples),
      expand("%s/QUANTIFICATION/PRODIGAL_FC/final.contigs_group_{cagroup}/{samples}_group_{cagroup}_fc.txt" % (config["project-folder"]), samples=samples, cagroup=assemblyGroups),
      expand("%s/BAM/final.contigs_full/{samples}_mega.bam.flagstat" % (config["project-folder"]), samples=samples),
      expand("%s/BAM/final.contigs_full/{samples}_mega.bam.coverage" % (config["project-folder"]), samples=samples),
      expand("%s/BAM/final.contigs_coas{cagroup}/{samples}_mega.bam.flagstat" % (config["project-folder"]), samples=samples, cagroup=assemblyGroups),
      expand("%s/BAM/final.contigs_coas{cagroup}/{samples}_mega.bam.coverage" % (config["project-folder"]), samples=samples, cagroup=assemblyGroups),
      expand("%s/CONCOCT/COVERAGE_TABLE/coverage_table_{fafilter}.tsv" % (config["project-folder"]), fafilter=fafilter),
      expand("%s/CONCOCT/COVERAGE_TABLE/coverage_table_group_{cagroup}_{fafilter}.tsv" % (config["project-folder"]), cagroup=assemblyGroups, fafilter=fafilter),
      expand("%s/CONCOCT/FASTA_BINS/fasta_bins_{fafilter}" % (config["project-folder"]), fafilter=fafilter),
      expand("%s/CONCOCT/FASTA_BINS/fasta_bins_group_{cagroup}.{fafilter}" % (config["project-folder"]), cagroup=assemblyGroups, fafilter=fafilter)


#      "%s/CONCOCT/final.contigs.1k_10K.fa" % (config["project-folder"]),
#      "%s/CONCOCT/final.contigs.2k_10K.fa" % (config["project-folder"]),
#      "%s/CONCOCT/coverage_table_1k.tsv" % (config["project-folder"]),
#      "%s/CONCOCT/coverage_table_2k.tsv" % (config["project-folder"]),
#      "%s/CONCOCT/1k_clustering_gt1000.csv" % (config["project-folder"]),
#      "%s/CONCOCT/2k_clustering_gt1000.csv" % (config["project-folder"]),
#      "%s/CONCOCT/fasta_bins_1k" % (config["project-folder"]),
#      "%s/CONCOCT/fasta_bins_2k" % (config["project-folder"])

rule preparations:
    input:
        expand("%s/FASTQ/TRIMMED/{samples}_R1.trimmed.fastq.gz" % (config["project-folder"]), samples=samples)

rule qc:
    input:
        "%s/QC/RAW/multiqc_R1/" % (config["project-folder"]),
#        "%s/QC/CONCATENATED/multiqc_R1/" % (config["project-folder"]),
        "%s/QC/TRIMMED/multiqc_R1/" % (config["project-folder"]),
        expand("%s/QC/RAW/{rawsamples}_R1_qualdist.txt" % (config["project-folder"]), rawsamples=rawsamples),
#        expand("%s/QC/CONCATENATED/{samples}_R1_qualdist.txt" % (config["project-folder"]), samples=samples),
        expand("%s/QC/TRIMMED/{samples}_R1_qualdist.txt" % (config["project-folder"]), samples=samples)

rule decontaminate:
    input:
        expand("%s/FASTQ/DECONTAMINATED/{samples}_R1.decontaminated.fastq.gz" % (config["project-folder"]), samples=samples)

rule combineFastq:
    input:
        "%s/FASTQ/MERGED/all_merged_R1.fastq.gz" % (config["project-folder"]),
        "%s/FASTQ/MERGED/all_merged_R2.fastq.gz" % (config["project-folder"])

rule createMetagenome:
    input:
        "%s/MEGAHIT/final.contigs.fa" % (config["project-folder"]),
        "%s/QUAST/report.html" % (config["project-folder"]),
        expand("%s/MEGAHIT/final.contigs.group_{cagroup}.fa" % (config["project-folder"]), cagroup=assemblyGroups),
        "%s/MEGAHIT/final.contigs.1k.fa" % (config["project-folder"]),
        "%s/MEGAHIT/final.contigs.2k.fa" % (config["project-folder"]),
        expand("%s/MEGAHIT/final.contigs.group_{cagroup}.1k.fa" % (config["project-folder"]), cagroup=assemblyGroups),
        expand("%s/MEGAHIT/final.contigs.group_{cagroup}.2k.fa" % (config["project-folder"]), cagroup=assemblyGroups)

rule alignment:
    input:
        expand("%s/BAM/final.contigs_full/{samples}_mega.bam" % (config["project-folder"]), samples=samples),
        expand("%s/BAM/final.contigs_coas{cagroup}/{samples}_mega.bam" % (config["project-folder"]), cagroup=assemblyGroups, samples=samples)

rule genePrediction:
    input:
        "%s/PRODIGAL/final.contigs_full/final.contigs_full.prodigal.fa" % (config["project-folder"]),
        expand("%s/PRODIGAL/final.contigs_group_{cagroup}/final.contigs_group_{cagroup}.prodigal.fa" % (config["project-folder"]), cagroup=assemblyGroups),
        "%s/EGGNOG/final.contigs_full/eggnog_output.emapper.annotations" % (config["project-folder"]),
        expand("%s/EGGNOG/final.contigs_group_{cagroup}/eggnog_output.emapper.annotations" % (config["project-folder"]), cagroup=assemblyGroups)
        
rule quantification:
    input:
        expand("%s/QUANTIFICATION/PRODIGAL_FC/final.contigs_full/{samples}_full_fc.txt" % (config["project-folder"]), samples=samples),
        expand("%s/QUANTIFICATION/PRODIGAL_FC/final.contigs_group_{cagroup}/{samples}_group_{cagroup}_fc.txt" % (config["project-folder"]), samples=samples, cagroup=assemblyGroups),
        expand("%s/BAM/final.contigs_full/{samples}_mega.bam.flagstat" % (config["project-folder"]), samples=samples),
        expand("%s/BAM/final.contigs_full/{samples}_mega.bam.coverage" % (config["project-folder"]), samples=samples),
        expand("%s/BAM/final.contigs_coas{cagroup}/{samples}_mega.bam.flagstat" % (config["project-folder"]), samples=samples, cagroup=assemblyGroups),
        expand("%s/BAM/final.contigs_coas{cagroup}/{samples}_mega.bam.coverage" % (config["project-folder"]), samples=samples, cagroup=assemblyGroups)
        
rule mags:
    input:
        expand("%s/CONCOCT/COVERAGE_TABLE/coverage_table_{fafilter}.tsv" % (config["project-folder"]), fafilter=fafilter),
        expand("%s/CONCOCT/COVERAGE_TABLE/coverage_table_group_{cagroup}_{fafilter}.tsv" % (config["project-folder"]), cagroup=assemblyGroups, fafilter=fafilter),
        expand("%s/CONCOCT/FASTA_BINS/fasta_bins_{fafilter}" % (config["project-folder"]), fafilter=fafilter),
        expand("%s/CONCOCT/FASTA_BINS/fasta_bins_group_{cagroup}.{fafilter}" % (config["project-folder"]), cagroup=assemblyGroups, fafilter=fafilter)

### setup report #####
report: "report/workflow.rst"

##### load rules #####
include: "rules/Step1-Preparations.smk"
include: "rules/Step2-ReadProcessing.smk"
include: "rules/Step2b-Decontamination.smk"
include: "rules/Step2c-QC.smk"
include: "rules/Step3-CreateMetagenome.smk"
include: "rules/Step3b-Metagenome-QC.smk"
include: "rules/Step3c-ReadAlignments.smk"
include: "rules/Step4-GenePrediction.smk"
include: "rules/Step5-Quantification.smk"
include: "rules/Step6-MAGs.smk"
