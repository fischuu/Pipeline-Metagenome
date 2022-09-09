if config["contamination-folder"] == "":
    pass
else:
    rule build_contamination_indices:
        """
        Create genome index for contamination genomes (STAR).
        """
        output:
            "%s/indexCreationReady" % (config["contamination-folder"])
        log:
            "%s/logs/build_contamination_indices.log" % (config["project-folder"])
        benchmark:
            "%s/benchmark/build_contamination_indices.benchmark.tsv" % (config["project-folder"])
        resources:
            time=cluster["build_contamination_indices"]["time"],
            mem=cluster["build_contamination_indices"]["mem-per-cpu"]
        threads: cluster["build_contamination_indices"]["cpus-per-task"]
        params:
           contamination_folder=config["contamination-folder"],
           contamination_fasta=config["contamination-refs"],
           overhang=config["params"]["star"]["overhang"],
           limitGenomeGenerateRAM=config["params"]["star"]["limitGenomeGenerateRAM"]
        singularity: config["singularity"]["star"]
        shell:"""
            refs=({params.contamination_fasta})
            
            for i in "${{refs[@]}}"
            do
            : 
                mkdir -p {params.contamination_folder}/Index-"$i"
                STAR --runThreadN {threads} --runMode genomeGenerate --limitGenomeGenerateRAM={params.limitGenomeGenerateRAM} --genomeDir {params.contamination_folder}/Index-"$i" --genomeFastaFiles {params.contamination_folder}/"$i"
            done
            touch {output};
        """
        
if config["contamination-folder"] == "":
    pass
else:
    rule star_decontaminate_reads:
        """
        Map the samples to the reference genome (STAR).
        """
        input:
            "%s/indexCreationReady" % (config["contamination-folder"]),
            fastq=["%s/FASTQ/TRIMMED/{samples}_R1.trimmed.fastq.gz" % (config["project-folder"]),
                   "%s/FASTQ/TRIMMED/{samples}_R2.trimmed.fastq.gz" % (config["project-folder"])]
        output:
            fastq=["%s/FASTQ/DECONTAMINATED/{samples}_R1.decontaminated.fastq.gz" % (config["project-folder"]),
                   "%s/FASTQ/DECONTAMINATED/{samples}_R2.decontaminated.fastq.gz" % (config["project-folder"])]
        log:
            "%s/logs/star_decontaminate_reads.{samples}.log" % (config["project-folder"])
        benchmark:
            "%s/benchmark/star_decontaminate_reads.{samples}.benchmark.tsv" % (config["project-folder"])
        resources:
            time=cluster["star_decontaminate_reads"]["time"],
            mem=cluster["star_decontaminate_reads"]["mem-per-cpu"]
        threads: cluster["star_decontaminate_reads"]["cpus-per-task"]
        params:
            project_folder=config["project-folder"],
            contamination_folder=config["contamination-folder"],
            contamination_fasta=config["contamination-refs"]
        singularity: config["singularity"]["star"]
        shell:"""
        
            refs=({params.contamination_fasta});
            tmpFastq=({input.fastq});
            currentInput="${{tmpFastq[@]}}"

            for i in "${{refs[@]}}"
            do
            :
                DIRout={params.project_folder}/BAM/DECONTAMINATION/"$i";
                FILEout={params.project_folder}/BAM/DECONTAMINATION/"$i"/{wildcards.samples}.bam;
                currentIndex={params.contamination_folder}/Index-"$i";
                
                mkdir -p $DIRout

                STAR --genomeDir "$currentIndex" \
                    --readFilesIn $currentInput \
                    --readFilesCommand zcat \
                    --outSAMtype BAM SortedByCoordinate \
                    --outSAMunmapped Within \
                    --runThreadN {threads} \
                    --outReadsUnmapped Fastx \
                    --outFileNamePrefix {wildcards.samples}_"$i"_ 2> {log};
        
                mv {wildcards.samples}_"$i"_Aligned.sortedByCoord.out.bam "$FILEout"
                mv {wildcards.samples}_"$i"_Log.final.out {wildcards.samples}_"$i"_Log.progress.out {wildcards.samples}_"$i"_Log.out {wildcards.samples}_"$i"_SJ.out.tab "$DIRout"
         
                mv {wildcards.samples}_"$i"_Unmapped.out.mate1 {wildcards.samples}_"$i"_R1.decontaminated.fastq
                mv {wildcards.samples}_"$i"_Unmapped.out.mate2 {wildcards.samples}_"$i"_R2.decontaminated.fastq
                
                wc -l {wildcards.samples}_"$i"_R1.decontaminated.fastq > {wildcards.samples}_"$i"_R1.decontaminated.fastq.wcl
                wc -l {wildcards.samples}_"$i"_R2.decontaminated.fastq > {wildcards.samples}_"$i"_R2.decontaminated.fastq.wcl
                
                gzip {wildcards.samples}_"$i"_R1.decontaminated.fastq
                gzip {wildcards.samples}_"$i"_R2.decontaminated.fastq
                
                tmpFastq=({params.project_folder}/{wildcards.samples}_"$i"_R1.decontaminated.fastq.gz {params.project_folder}/{wildcards.samples}_"$i"_R2.decontaminated.fastq.gz)
                currentInput="${{tmpFastq[@]}}"
                lastI="$i"
            done
           
            mv {wildcards.samples}_"$lastI"_R1.decontaminated.fastq.gz {wildcards.samples}_R1.decontaminated.fastq.gz;
            mv {wildcards.samples}_"$lastI"_R2.decontaminated.fastq.gz {wildcards.samples}_R2.decontaminated.fastq.gz;

            mkdir -p {params.project_folder}/FASTQ/DECONTAMINATED;
            mv {wildcards.samples}_R1.decontaminated.fastq.gz {params.project_folder}/FASTQ/DECONTAMINATED;
            mv {wildcards.samples}*decontaminated.fastq.wcl {params.project_folder}/FASTQ/DECONTAMINATED;

            rm {wildcards.samples}_*_R1.decontaminated.fastq.gz;
            rm {wildcards.samples}_*_R2.decontaminated.fastq.gz;
        """
