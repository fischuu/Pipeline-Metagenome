# CHANGELOG

## PLANNED UPDATES FOR FUTURE VERSIONS

## UPDATES

### Development version 0.5.*:
--------------------------------------------------------------------------------
2: Added the templete to prepare the summary statistics

1: Added the final report template


### Release version 0.4:
--------------------------------------------------------------------------------

### Development version 0.3.*:
--------------------------------------------------------------------------------
29: Added sed command to remove comment lines from the prodigal gene prediction

28: Added the k-min option for megahit to avoid crashes in case of deep sequenced genomes

27: Adjusted the output paths for the MAG module

26: Fixed a bug in the EGGNOG annotation step, that prevented the annotation to happen...

25: Added the concoct steps for the co-assemblies

24: Some concoct intermediate files are set to temp now

23: Concoct module  for the full assembly finalised

22: Added the --cluster-cancel option to the start script

21: Added the filtered assemblies to the QUAST report

20: Added a length filter for the metagenome assemblies

19: Constrained the cagroup wildcard

18: Report reads per contig for full and coas

17: Added alignment statistics to the quantification

16: Quantification step added

15: Eggnog annotation added

14: Prodigal gene prediction added

13: Samples can now be used in several co-assemblies, separated by comma in samplesheet file

12: Server resource allocation names adjusted

11: Bowtie2 indices are now forced to be large to avoid index naming issues

10: Added rules to align reads against the various assemblies (full and co-assemblies)

9: Metaquast added to the pipeline

8: Group-wise co-assemblies are now also generated

7: Added the QC (FastQC and MultiQC) for the Raw samples and concatenated as well as trimmed ones

6: Example start and config files copy+paste error fixed

5: More adjustments to the default resource allocations (extract_fasta_concoct_1k and extract_fasta_concoct_2k)

4: MEGAHIT_temp folder is now temporary only, to save tons of disc space

3: Adjusted the defaults resources for concoct parts 1k and 2k

2: Critical bugfix, 2k concoct calculations were done with 1k output

1: Inital submission

### Release version 0.2, includes all changes from dev version 0.1.*
--------------------------------------------------------------------------------

### Development version 0.1.*:
--------------------------------------------------------------------------------
16: Default resource allocations adjusted for avoid crashes

15: rule extract_fasta_concoct for 1k and 2k added

14: rule run_concoct added for 1k and 2k

13: rule cut_filtered_contigs_concoct added

12: Step5-MAGs.smk ruleset started

11: rule quantify_predictedGenes_featureCounts added

10: rule filter_length_megahit added

9: Eggnot orthology rule added

8: Eggnog gene annotation added

7: Execution script adjusted to bind local nvme discs

6: Prodigal output is now cut in chunks for parallel processing

5: Eggnog databases are now downloaded automatically

4: Bugfix in bowtie2 index creation step fixed

3: create_index_bowtie2 write now its output into a logfile

2: Added the decontamination step

1: Inital submission
