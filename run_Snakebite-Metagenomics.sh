module load snakemake

pipelineFolder="/users/fischerd/git/Snakebite-Metagenomics"
projectFolder="/scratch/project_2001746/Metagenomics_Example"

export APPTAINER_TMPDIR="/scratch/project_2001746/tmp"
export APPTAINER_CACHEDIR="/scratch/project_2001746/tmp"

mkdir -p $APPTAINER_TMPDIR

# Create the rulegraph
snakemake -s $pipelineFolder/Snakebite-Metagenomics.smk \
          --configfile $projectFolder/Snakebite-Metagenomics_config.yaml \
          --rulegraph | dot -T png > $projectFolder/workflow.png

snakemake -s $pipelineFolder/Snakebite-Metagenomics.smk \
          -j 300 \
          --use-singularity \
          --singularity-args "-B /scratch:/scratch,/run:/run" \
          --configfile $projectFolder/Snakebite-Metagenomics_config.yaml \
          --latency-wait 60 \
          --cluster-config $pipelineFolder/Snakebite-Metagenomics_server_config.yaml \
          --cluster "sbatch -t {cluster.time} --account={cluster.account} --gres=nvme:{cluster.nvme} --job-name={cluster.job-name} --tasks-per-node={cluster.ntasks} --cpus-per-task={cluster.cpus-per-task} --mem-per-cpu={cluster.mem-per-cpu} -p {cluster.partition} -D {cluster.working-directory} --parsable" \
          --cluster-cancel scancel $@
