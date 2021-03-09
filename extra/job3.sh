#!/bin/bash -l

#SBATCH --job-name=Nextflow-master-nanopore
#SBATCH --account=pawsey0281
#SBATCH --partition=longq
#SBATCH --time=4-00:00:00
#SBATCH --no-requeue
#SBATCH --export=none

module load nextflow

nextflow run marcodelapierre/illumina-nf \
  --reads='MK27_D1YHBACXX_GTGAAA_L004_R{1,2}.fastq.gz' \
  --seqs='HG970869.1,JX173278.1,KF632713.1' \
  --contigs='NODE_1,NODE_2,NODE_3' \
  -profile zeus --slurm_account='pawsey0281' \
  -name nxf-${SLURM_JOB_ID} \
  -with-trace trace-${SLURM_JOB_ID}.txt \
  -with-report report-${SLURM_JOB_ID}.html