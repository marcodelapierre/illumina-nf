#!/bin/bash

. use-nextflow.sh

gunzip -k tinydb.fasta.gz
singularity exec docker://quay.io/biocontainers/blast:2.12.0--pl5262h3289130_0 makeblastdb -in tinydb.fasta -dbtype nucl -parse_seqids

for profile in joey1 joey2 joey3; do

  nextflow run main.nf \
    --reads='small_R{1,2}.fastq.gz' \
    --blast_db="$(pwd)/tinydb.fasta" \
    -profile $profile

done
