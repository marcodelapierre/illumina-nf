#!/bin/bash -l

#SBATCH --job-name=nxf-small
#SBATCH --account=pawsey0001
#SBATCH --partition=workq
#SBATCH --time=10:00
#SBATCH --no-requeue
#SBATCH --export=none

### NOTE: may note longer runtime above, if container images have yet to be pulled

unset SBATCH_EXPORT

module load singularity
module load nextflow

gunzip -k tinydb.fasta.gz
singularity exec docker://quay.io/biocontainers/blast:2.7.1--h96bfa4b_5 makeblastdb -in tinydb.fasta -dbtype nucl -parse_seqids

nextflow run main.nf \
  --reads='small_R{1,2}.fastq.gz' \
  --blast_db="$(pwd)/tinydb.fasta" \
  -profile zeus --slurm_account='pawsey0001' \
  -name nxf-${SLURM_JOB_ID}
