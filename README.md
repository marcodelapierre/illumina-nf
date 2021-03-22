## Illumina pipeline for DPIRD - Nextflow edition

The pipeline requires [Nextflow](https://github.com/nextflow-io/nextflow) to run.  
DSL2 syntax is used, so that Nextflow version `20.07.1` or higher is required.


### Pipeline

This pipeline is the Nextflow translation of all *Illumina* workflows in the project https://github.com/pawseySC/dpird-mk:  

Merge\*(+QC) -> Trim(+QC) -> De-novo assemble -> (Map contigs && Blast) -> Map ref. sequences\+\# -> Align\#

\* Can be *Interleave* instead  
\+ Can be performed in *cascade* mode  
\# Require additional inputs in subsequent runs  

Note how some of these steps map to multiple Nextflow processes, to separate executions that use distinct packages (which is useful when using containerised software).


### Basic usage

```
nextflow run marcodelapierre/illumina-nf \
  --reads='reads_{1,2}.fastq.gz' \
  -profile zeus --slurm_account='<Your Pawsey Project>'
```

The flag `--reads` specifies the name of the pair of input read files.  
Note some syntax requirements:
- encapsulate the file name specification between single quotes;
- within a file pair, use names that differ only by a character, which distinguishes the two files, in this case `1` or `2`;
- use curly brackets to specify the wild character within the file pair, *e.g.* `{1,2}`;
- the prefix to the wild character serves as the sample ID, *e.g.* `reads_`.

Output files are stored in subdirectory(ies) with name `results_$sampleID`.  Reference sequences that you select for mapping and/or alignment are stored in the subdirectory `refseqs`.

The flag `-profile` (note the single dash) allows to select the appropriate profile for the machine in use, Zeus in this case.  The flag `--slurm_account` sets your Pawsey account to run on Zeus.  

After blasting and identifying reference sequences of interest, mapping of input reads against these sequences can be performed, providing the sequence IDs via the flag `--seqs` (case insensitive, comma separated list):

```
nextflow run marcodelapierre/illumina-nf \
  --reads='reads_{1,2}.fastq.gz' \
  --seqs='HG970869.1,JX173278.1,HG970865.1,KF632713.1' \
  -profile zeus --slurm_account='<Your Pawsey Project>'
```

Finally, after selecting contigs of interest from the assembly, alignment with the reference sequences can be performed, by additionally providing the contig IDs via the flag `--contigs` (case insensitive, comma separated list), *e.g.*:

```
nextflow run marcodelapierre/illumina-nf \
  --reads='reads_{1,2}.fastq.gz' \
  --seqs='HG970869.1,JX173278.1,KF632713.1' \
  --contigs='NODE_1,NODE_2,NODE_3' \
  -profile zeus --slurm_account='<Your Pawsey Project>'
```


### Pipeline variants

1. Depending on the nature of the input files, if you need to ***interleave*** rather than to *merge* the read pairs you can achieve so by using the optional flag `--interleave`.

2. If you need to use the ***reverse-complement*** of a reference sequence or contig, just append the suffix `/rc`, or `_rc`, to its ID, as in `HG970869.1/rc` or `NODE_1/rc`.

3. By default, read mapping to reference sequences is performed independently for each sequence from the input list.  In alternative, a ***cascade mode*** is available for read mapping, where mapping is performed in an ordered, progressive way.  Here, the first read mapping is performed against the first sequence in the input list using all input reads;  then the second read mapping is performed against the second sequence, by considering only those reads that went unmapped in the previous mapping step;  the process continues until all input sequences are processed.  
   To activate cascade mapping use the optional flag `--cascade`.  Note how the *ordered set* of input reference sequences is an information that needs to be retained at the alignment step, to uniquely identify the cascade mapping process.  For this reason, you need to provide two sets of reference sequences for the alignment step:  as an example, `--seqs='SEQ1,SEQ2,SEQ3,SEQ4'` has the full set of sequences from the cascade mapping step, whereas the new flag `--cascade_align_seqs='SEQ1,SEQ3'` has the subset of sequences that have been chosen for the alignment.


### Optional parameters

* Change *evalue* for blasting: `--evalue='0.1'`.
* Change minimum length threshold for assembled contigs to be considered for blasting: `--min_len_contig='1000'`.


### Multiple inputs at once

The pipeline allows to feed in multiple datasets at once.  You can use input file name patterns to this end:

1. multiple input read pairs in the same directory, *e.g.* `sample1_R{1,2}.fq`, `sample2_R{1,2}.fq` and so on, use: `--reads='sample*{1,2}.fq'`;

2. multiple read pairs in distinct directories, *e.g.* `sample1/R{1,2}.fq`, `sample2/R{1,2}.fq` and so on, use: `--reads='sample*/R{1,2}.fq'`.


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

The `extra` directory contains example Slurm scripts, `job1.sh`, `job2.sh`, `job3.sh`, `job2_cascade.sh` and `job3_cascade.sh`, to run on Zeus.  There is also a sample script `log.sh` that takes a run name as input and displays formatted runtime information.

The `test` directory contains a small input dataset and two launching scripts (based on Zeus) for quick testing of the pipeline, with total runtime below ten minutes.
