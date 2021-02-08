#!/usr/bin/env nextflow

nextflow.enable.dsl=2


process merge {
  input:
  
  output:
  
  script:
  """
  bbmerge.sh \
    in1=R1.fastq.gz in2=R2.fastq.gz \
    out=merged.fastq.gz
  """
}


process qc_post_merge {
  input:
  
  output:
  
  script:
  """
  fastqc merged.fastq.gz
  """
}


process trim {
  input:
  
  output:
  
  script:
  """
  bbduk.sh \
    in=merged.fastq.gz \
    out=trimmed-partial.fastq.gz \
    ref=adapters ktrim=r k=27 hdist=2 edist=0 mink=4

  bbduk.sh \
    in=trimmed-partial.fastq.gz \
    out=clean.fastq.gz \
    ref=adapters ktrim=l k=27 hdist=2 edist=0 mink=4 \
    qtrim=rl trimq=13 \
    minlength=30
  """
}


process qc_post_trim {
  input:
  
  output:
  
  script:
  """
  fastqc clean.fastq.gz
  """
}


process assemble {
  input:
  
  output:
  
  script:
  """
  spades.py \
    -s clean.fastq.gz \
    --only-assembler \
    -t $OMP_NUM_THREADS -m $((SLURM_MEM_PER_NODE/1024)) \
    -o .

  min_len_contig="1000"

  awk -v min_len_contig=$min_len_contig -F _ '{ if( \$1 == ">NODE" ){ if( \$4 < min_len_contig ) {exit} } ; print }' contigs.fasta >contigs_sub.fasta
  """
}


process map_contigs {
  input:
  
  output:
  
  script:
  """
  bbmap.sh \
    in=clean.fastq.gz ref=contigs_sub.fasta \
    out=mapped_contigs_sub_unsorted.sam \
    k=13 maxindel=16000 ambig=random \
    threads=$OMP_NUM_THREADS
  """
}


process sam_post_map {
  input:
  
  output:
  
  script:
  """
  samtools \
    view -b -o mapped_contigs_sub_unsorted.bam mapped_contigs_sub_unsorted.sam

  samtools \
    sort -o mapped_contigs_sub.bam mapped_contigs_sub_unsorted.bam

  samtools \
    index mapped_contigs_sub.bam

  samtools \
    depth -aa mapped_contigs_sub.bam >depth_contigs_sub.dat
  """
}


process bcf_post_map {
  input:
  
  output:
  
  script:
  """
  bcftools \
    mpileup -Ou -f contigs_sub.fasta mapped_contigs_sub.bam \
    | bcftools \
    call --ploidy 1 -mv -Oz -o calls_contigs_sub.vcf.gz

  bcftools \
    tabix calls_contigs_sub.vcf.gz

  bcftools \
    consensus -f contigs_sub.fasta -o consensus_contigs_sub.fasta calls_contigs_sub.vcf.gz
  """
}


process sam_pre_blast {
  input:
  
  output:
  
  script:
  """
  samtools faidx \
    -i -o contigs_sub_rc.fasta \
    contigs_sub.fasta \$(grep '^>' contigs_sub.fasta | tr -d '>')
  """
}


process blast {
  input:
  
  output:
  
  script:
  """
  blastn \
    -query ${prefix_contig}.fasta -db /group/data/blast/nt \
    -outfmt 11 -out blast_${prefix_contig}.asn \
    -max_hsps 50 \
    -word_size 28 -evalue 0.1 \
    -reward 1 -penalty -2 \
    -num_threads $OMP_NUM_THREADS

  blast_formatter \
    -archive blast_${prefix_contig}.asn \
    -outfmt 5 -out blast_${prefix_contig}.xml

  blast_formatter \
    -archive blast_${prefix_contig}.asn \
    -outfmt 6 -out blast_unsort_default_${prefix_contig}.tsv

blast_formatter \
    -archive blast_${prefix_contig}.asn \
    -outfmt "6 qaccver saccver pident length evalue bitscore stitle" -out blast_unsort_${prefix_contig}.tsv
 
 sort -n -r -k 6 blast_unsort_${prefix_contig}.tsv >blast_${prefix_contig}.tsv
 



  """
}


process map_refseq {
  input:
  
  output:
  
  script:
  """
  
  """
}


process align {
  input:
  
  output:
  
  script:
  """
  
  """
}


workflow {

  

}