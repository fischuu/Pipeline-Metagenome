# This conda module is just to make snakemake available
module load bioconda/3
source activate Snakemake

export SINGULARITY_TMPDIR="/scratch/project_xxx/tmp"
export SINGULARITY_CACHEDIR="/scratch/project_xxx/tmp"

# Create the rulegraph
snakemake -s ~/git/Pipeline-Metagenome/Snakefile-Pipeline-Metagenome.smk \
          --configfile /scratch/project_XXX/Pipeline-WGS_config.yaml \
          --rulegraph | dot -T png > /scratch/project_XXX/workflow.png

snakemake -s ~/git/Pipeline-Metagenome/Snakefile-Pipeline-Metagenome.smk \
          -j 300 \
          --use-singularity \
          --singularity-args "-B /scratch:/scratch,/run:/run" \
          --configfile /scratch/project_XXX/Pipeline-WGS_config.yaml \
          --latency-wait 60 \
          --cluster-config ~/git/Pipeline-Metagenome/Pipeline-Metagenome_server_config.yaml \
          --cluster "sbatch -t {cluster.time} --account={cluster.account} --gres=nvme:{cluster.nvme} --job-name={cluster.job-name} --tasks-per-node={cluster.ntasks} --cpus-per-task={cluster.cpus-per-task} --mem-per-cpu={cluster.mem-per-cpu} -p {cluster.partition} -D {cluster.working-directory}" \
          $@
