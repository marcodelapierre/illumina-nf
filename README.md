## UNDER CONSTRUCTION - Illumina pipeline for DPIRD - Nextflow edition

The pipeline requires [Nextflow](https://github.com/nextflow-io/nextflow) to run.  
DSL2 syntax is used, so that Nextflow version `20.07.1` or higher is recommended.


### Pipeline

Merge(+QC) -> Trim(+QC) -> De-novo assemble -> (Map contigs && Blast) -> Map ref. sequences\# -> Align\#

\# Requires additional input in a subsequent run


### Basic usage

```
nextflow run marcodelapierre/illumina-nf \
  --reads='reads_{1,2}.fastq.gz' \
  -profile zeus --slurm_account='pawsey0001'
```

The flag `--reads` specifies the name of the pair of input read files.  
Note some syntax requirements:
- encapsulate the file name specification between single quotes;
- within a file pair, use names that differ only by a character, which distinguishes the two files, in this case `1` or `2`;
- use curly brackets to specify the wild character within the file pair, *e.g.* `{1,2}`;
- the prefix to the wild character serves as the sample ID, *e.g.* `reads_`.

Output files are stored in subdirectory(ies) with name `results_$sampleID`.  Reference sequences that you select for mapping and/or alignment are stored in the subdirectory `refseqs`.

The flag `-profile` allows to select the appropriate profile for the machine in use, Zeus in this case.  The flag `--slurm_account` sets your Pawsey account to run on Zeus.  

After blasting and identifying reference sequences of interest, mapping of input reads against these sequences can be performed, providing the sequence IDs via the flag `--seqs` (comma separated list):

```
nextflow run marcodelapierre/illumina-nf \
  --reads='reads_{1,2}.fastq.gz' \
  --seqs='comma,separated,list,of,ids,from,blast' \
  -profile zeus --slurm_account='pawsey0001'
```

Finally, after selecting contigs of interest from the assembly, alignment with the reference sequences can be performed, by additionally providing the contig IDs via the flag `--contigs` (comma separated list), *e.g.*:

```
nextflow run marcodelapierre/illumina-nf \
  --reads='reads_{1,2}.fastq.gz' \
  --seqs='comma,separated,list,of,ids,from,blast' \
  --contigs='NODE_1,NODE_234,NODE_56' \
  -profile zeus --slurm_account='pawsey0001'
```


### Multiple inputs at once

The pipeline allows to feed in multiple datasets at once.  You can use input file name patterns to this end:

1. multiple input read pairs in the same directory, *e.g.* `sample1_R{1,2}.fq`, `sample2_R{1,2}.fq` and so on, use: `--reads='sample*{1,2}.fq'`;

2. multiple read pairs in distinct directories, *e.g.* `sample1/R{1,2}.fq`, `sample2/R{1,2}.fq` and so on, use: `--reads='sample*/R{1,2}.fq'`.


### Optional parameters

* Change *evalue* for blasting: `--evalue='0.1'`.
* Change minimum length threshold for assembled contigs to be considered for blasting: `--min_len_contig='1000'`.


### Requirements

Setup:
* Singularity or Docker -- if you wish to use containerised software

Software:
* BBmap
* FastQC
* Spades
* SAMtools
* BCFtools
* Blast
* MAFFT

Reference data:
* Database for Blast


### Additional resources

The `extra` directory contains example Slurm scripts, `job1.sh`,  `job2.sh` and `job3.sh`, to run on Zeus.  There is also a sample script `log.sh` that takes a run name as input and displays formatted runtime information.

The `test` directory contains a small input dataset and a launching script for quick testing of the pipeline, with total runtime around five minutes.


