# CHANGELOG

## PLANNED UPDATES FOR FUTURE VERSIONS

## UPDATES

### Development version 0.3.*:
--------------------------------------------------------------------------------
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
