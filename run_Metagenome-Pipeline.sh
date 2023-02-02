# This conda module is just to make snakemake available
module load Snakemake/7.17.1

pipelineFolder="/users/XXX/git/Pipeline-GBS"
projectFolder="/scratch/project_XXX/"

# Create the rulegraph
snakemake -s $pipelineFolder/Snakefile-Pipeline-Metagenome.smk \
          --configfile $projectFolder/Pipeline-Metagenome_config.yaml \
          --rulegraph | dot -T png > $projectFolder/workflow.png

snakemake -s $pipelineFolder/Snakefile-Pipeline-Metagenome.smk \
          -j 300 \
          --use-singularity \
          --singularity-args "-B /scratch:/scratch,/run:/run" \
          --configfile $projectFolder/Pipeline-Metagenome_config.yaml \
          --latency-wait 60 \
          --cluster-config $pipelineFolder/Pipeline-Metagenome_server_config.yaml \
          --cluster "sbatch -t {cluster.time} --account={cluster.account} --gres=nvme:{cluster.nvme} --job-name={cluster.job-name} --tasks-per-node={cluster.ntasks} --cpus-per-task={cluster.cpus-per-task} --mem-per-cpu={cluster.mem-per-cpu} -p {cluster.partition} -D {cluster.working-directory} --parsable" \
          --cluster-cancel scancel $@
