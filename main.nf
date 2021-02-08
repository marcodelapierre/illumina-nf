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


process sam_post_map_contigs {
  input:
  
  output:
  
  script:
  """
  samtools \
    view -b -o mapped_contigs_sub_unsorted.bam \
    mapped_contigs_sub_unsorted.sam

  samtools \
    sort -o mapped_contigs_sub.bam \
    mapped_contigs_sub_unsorted.bam

  samtools \
    index mapped_contigs_sub.bam

  samtools \
    depth -aa mapped_contigs_sub.bam \
    >depth_contigs_sub.dat
  """
}


process bcf_post_map_contigs {
  input:
  
  output:
  
  script:
  """
  bcftools \
    mpileup -Ou -f contigs_sub.fasta \
    mapped_contigs_sub.bam \
    | bcftools \
    call --ploidy 1 -mv -Oz \
    -o calls_contigs_sub.vcf.gz

  bcftools \
    tabix calls_contigs_sub.vcf.gz

  bcftools \
    consensus -f contigs_sub.fasta \
    -o consensus_contigs_sub.fasta \
    calls_contigs_sub.vcf.gz
  """
}


process sam_pre_blast {
  input:
  
  output:
  
  script:
  """
  samtools faidx \
    -i -o contigs_sub_revcom.fasta \
    contigs_sub.fasta \$(grep '^>' contigs_sub.fasta | tr -d '>')
  """
}


process blast {
  input:
  
  output:
  
  script:
  """
  blastn \
    -query ${prefix_contig}.fasta -db $blast_db_path \
    -outfmt 11 -out blast_${prefix_contig}.asn \
    -max_hsps 50 \
    -word_size 28 -evalue 0.1 \
    -reward 1 -penalty -2 \
    -num_threads $OMP_NUM_THREADS

  blast_formatter \
    -archive blast_${prefix_contig}.asn \
    -outfmt 5 \
    -out blast_${prefix_contig}.xml

  blast_formatter \
    -archive blast_${prefix_contig}.asn \
    -outfmt "6 qaccver saccver pident length evalue bitscore stitle" \
    -out blast_unsort_${prefix_contig}.tsv
 
  sort -n -r -k 6 blast_unsort_${prefix_contig}.tsv >blast_${prefix_contig}.tsv
 """
}


process seqfile {
  input:
  
  output:
  
  script:
  """
  blastdbcmd \
    -db $blast_db_path -entry ${seqid%/rc} \
    -line_length 60 \
    -out refseq_${MID}.fasta

  sed -i '/^>/ s/ .*//g' refseq_${MID}.fasta
  """
}


process sam_post_seqfile {
  input:
  
  output:
  
  script:
  """
  if [ "${seqid: -3}" == "/rc" ] ; then
  samtools faidx \
	-i -o refseq_${MID}_rc.fasta \
	refseq_${MID}.fasta ${seqid%/rc}
   mv refseq_${MID}_rc.fasta refseq_${MID}.fasta
  fi

  sed -i '/^>/ s/ .*//g' refseq_${MID}.fasta
  """
}


process map_refs {
  input:
  
  output:
  
  script:
  """
  bbmap.sh \
    in=clean.fastq.gz ref=refseq_${MID}.fasta \
    out=mapped_refseq_${MID}_unsorted.sam \
    k=13 maxindel=16000 ambig=random \
    path=ref_${MID} \
    threads=$OMP_NUM_THREADS
  """
}


process sam_post_map_refs {
  input:
  
  output:
  
  script:
  """
  samtools \
    view -b -o mapped_refseq_${MID}_unsorted.bam \
    mapped_refseq_${MID}_unsorted.sam

samtools \
    sort -o mapped_refseq_${MID}.bam \
    mapped_refseq_${MID}_unsorted.bam

samtools \
    index mapped_refseq_${MID}.bam

samtools \
    depth -aa mapped_refseq_${MID}.bam \
    >depth_refseq_${MID}.dat
  """
}


process bcf_post_map_refs {
  input:
  
  output:
  
  script:
  """
  bcftools \
    mpileup -Ou -f refseq_${MID}.fasta \
    mapped_refseq_${MID}.bam \
    | bcftools \
    call --ploidy 1 -mv -Oz \
    -o calls_refseq_${MID}.vcf.gz

bcftools \
    tabix calls_refseq_${MID}.vcf.gz

bcftools \
    consensus -f refseq_${MID}.fasta \
    -o consensus_refseq_${MID}.fasta \
    calls_refseq_${MID}.vcf.gz

  """
}



process align {
  input:
  
  output:
  
  script:
  """
  echo Dummy
  """
}


workflow {

  

}