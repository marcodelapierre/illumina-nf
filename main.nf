#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.reads='R{1,2}.fastq.gz'
params.seqs=''
params.nodes=''

params.min_len_contig='1000'
params.evalue='0.1'
params.outprefix='results_'
params.refdir='refseqs'


process merge_reads {
  tag "${dir}/${name}"
  publishDir "${dir}/${params.outprefix}${name}", mode: 'copy'
  stageInMode 'symlink'

  input:
  tuple val(dir), val(name), path(read1), path(read2)

  output:
  tuple val(dir), val(name), path('merged.fastq.gz')

  script:
  """
  bbmerge.sh \
    in1=$read1 in2=$read2 \
    out=merged.fastq.gz
  """
}


process qc_post_merge {
  tag "${dir}/${name}"
  publishDir "${dir}/${params.outprefix}${name}", mode: 'copy'

  input:
  tuple val(dir), val(name), path('merged.fastq.gz')

  output:
  tuple val(dir), val(name), path('merged_fastqc.html'), path('merged_fastqc.zip')

  script:
  """
  fastqc merged.fastq.gz
  """
}


process trim {
  tag "${dir}/${name}"
  publishDir "${dir}/${params.outprefix}${name}", mode: 'copy'

  input:
  tuple val(dir), val(name), path('merged.fastq.gz')

  output:
  tuple val(dir), val(name), path('clean.fastq.gz')

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
  tag "${dir}/${name}"
  publishDir "${dir}/${params.outprefix}${name}", mode: 'copy'

  input:
  tuple val(dir), val(name), path('clean.fastq.gz')

  output:
  tuple val(dir), val(name), path('clean_fastqc.html'), path('clean_fastqc.zip')

  script:
  """
  fastqc clean.fastq.gz
  """
}


process assemble {
  tag "${dir}/${name}"
  publishDir "${dir}/${params.outprefix}${name}", mode: 'copy'

  input:
  tuple val(dir), val(name), path('clean.fastq.gz')

  output:
  tuple val(dir), val(name), path('contigs_sub.fasta'), emit: sub
  tuple val(dir), val(name), path('contigs_ALL.fasta'), emit: all

  script:
  """
  memory='${task.memory}'
  memory=\${memory% *}

  spades.py \
    -s clean.fastq.gz \
    --only-assembler \
    -t ${task.cpus} -m \${memory} \
    -o .

  awk -v min_len_contig=${params.min_len_contig} -F _ '{ if( \$1 == ">NODE" ){ if( \$4 < min_len_contig ) {exit} } ; print }' contigs.fasta >contigs_sub.fasta

  mv contigs.fasta contigs_ALL.fasta
  """
}


process map_contigs {
  tag "${dir}/${name}"

  input:
  tuple val(dir), val(name), path('clean.fastq.gz'), path('contigs_sub.fasta')

  output:
  tuple val(dir), val(name), path('mapped_contigs_sub_unsorted.sam')

  script:
  """
  bbmap.sh \
    in=clean.fastq.gz \
    ref=contigs_sub.fasta \
    out=mapped_contigs_sub_unsorted.sam \
    k=13 maxindel=16000 ambig=random \
    threads=${task.cpus}
  """
}


process sam_post_map_contigs {
  tag "${dir}/${name}"
  publishDir "${dir}/${params.outprefix}${name}", mode: 'copy'

  input:
  tuple val(dir), val(name), path('mapped_contigs_sub_unsorted.sam')

  output:
  tuple val(dir), val(name), path('mapped_contigs_sub.bam'), path('mapped_contigs_sub.bam.bai'), emit: bam
  tuple val(dir), val(name), path('depth_contigs_sub.dat'), emit: depth

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
  tag "${dir}/${name}"
  publishDir "${dir}/${params.outprefix}${name}", mode: 'copy'

  input:
  tuple val(dir), val(name), path('mapped_contigs_sub.bam'), path('mapped_contigs_sub.bam.bai'), path('contigs_sub.fasta')

  output:
  tuple val(dir), val(name), path('calls_contigs_sub.vcf.gz'), emit: call
  tuple val(dir), val(name), path('consensus_contigs_sub.fasta'), emit: cons

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


process blast {
  tag "${dir}/${name}"
  publishDir "${dir}/${params.outprefix}${name}", mode: 'copy'

  input:
  tuple val(dir), val(name), path('contigs_sub.fasta')

  output:
  tuple val(dir), val(name), path('blast_contigs_sub.tsv'), emit: tsv
  tuple val(dir), val(name), path('blast_contigs_sub.xml'), emit: xml

  script:
  """
  blastn \
    -query contigs_sub.fasta -db ${params.blast_db} \
    -outfmt 11 -out blast_contigs_sub.asn \
    -max_hsps 50 \
    -word_size 28 -evalue ${params.evalue} \
    -reward 1 -penalty -2 \
    -num_threads ${task.cpus}

  blast_formatter \
    -archive blast_contigs_sub.asn \
    -outfmt 5 \
    -out blast_contigs_sub.xml

  blast_formatter \
    -archive blast_contigs_sub.asn \
    -outfmt "6 qaccver saccver pident length evalue bitscore stitle" \
    -out blast_unsort_contigs_sub.tsv
 
  sort -n -r -k 6 blast_unsort_contigs_sub.tsv >blast_contigs_sub.tsv
 """
}


process seqfile {
  tag "${seqid}"
  publishDir "${params.refdir}/", mode: 'copy', saveAs: { filename -> "refseq_${seqid}.fasta" }

  input:
  val(seqid)

  output:
  tuple val(seqid), path('refseq.fasta')

  script:
  """
  seqid="${seqid}"
  blastdbcmd \
    -db ${params.blast_db} -entry \${seqid} \
    -line_length 60 \
    -out refseq.fasta

  sed -i '/^>/ s/ .*//g' refseq.fasta
  """
}


// process sam_post_seqfile {
//   tag "${seqid}"
//   publishDir "${params.refdir}/", mode: 'copy', saveAs: { filename -> "refseq_${seqid}.fasta" }
//   input:
  
//   output:
  
//   script:
//   """
//   seqid="${seqid}"
//   seqid="\${seqid//_rc/\/rc}"
//   if [ "\${seqid: -3}" == "/rc" ] ; then
//     samtools faidx \
//       -i -o refseq_revcom.fasta \
//       refseq.fasta \${seqid%/rc}
//     mv refseq_revcom.fasta refseq.fasta
//     sed -i '/^>/ s/ .*//g' refseq.fasta
//   fi
//   """
// }


process map_refs {
  tag "${dir}/${name}_${seqid}"

  input:
  tuple val(dir), val(name), path('clean.fastq.gz'), val(seqid), path('refseq.fasta')

  output:
  tuple val(dir), val(name), val(seqid), path('mapped_refseq_unsorted.sam')
  script:
  """
  bbmap.sh \
    in=clean.fastq.gz \
    ref=refseq.fasta \
    out=mapped_refseq_unsorted.sam \
    k=13 maxindel=16000 ambig=random \
    threads=${task.cpus}
  """
}


process sam_post_map_refs {
  tag "${dir}/${name}_${seqid}"
  publishDir "${dir}/${params.outprefix}${name}/${seqid}/", mode: 'copy'
//  publishDir "${dir}/${params.outprefix}${name}/", mode: 'copy', saveAs: { filename -> filename.replaceFirst(/_refseq/,"_refseq_$seqid") }

  input:
  tuple val(dir), val(name), val(seqid), path('mapped_refseq_unsorted.sam')

  output:
// NOTE: val(seqid) first here, as needed by the cross operator
  tuple val(seqid), val(dir), val(name), path('mapped_refseq.bam'), path('mapped_refseq.bam.bai'), emit: bam
  tuple val(dir), val(name), val(seqid), path('depth_refseq.dat'), emit: depth

  script:
  """
  samtools \
    view -b -o mapped_refseq_unsorted.bam \
    mapped_refseq_unsorted.sam

samtools \
    sort -o mapped_refseq.bam \
    mapped_refseq_unsorted.bam

samtools \
    index mapped_refseq.bam

samtools \
    depth -aa mapped_refseq.bam \
    >depth_refseq.dat
  """
}


process bcf_post_map_refs {
  tag "${dir}/${name}_${seqid}"
  publishDir "${dir}/${params.outprefix}${name}/${seqid}/", mode: 'copy'

  input:
  tuple val(dir), val(name), val(seqid), path('mapped_refseq.bam'), path('mapped_refseq.bam.bai'), path('refseq.fasta')
  
  output:
  tuple val(dir), val(name), val(seqid), path('calls_refseq.vcf.gz'), emit: call
  tuple val(dir), val(name), val(seqid), path('consensus_refseq.fasta'), emit: cons

  script:
  """
  bcftools \
    mpileup -Ou -f refseq.fasta \
    mapped_refseq.bam \
    | bcftools \
    call --ploidy 1 -mv -Oz \
    -o calls_refseq.vcf.gz

bcftools \
    tabix calls_refseq.vcf.gz

bcftools \
    consensus -f refseq.fasta \
    -o consensus_refseq.fasta \
    calls_refseq.vcf.gz

  """
}


// process nodefile {
//   input:
  
//   output:
  
//   script:
//   """
//   echo Dummy
//   """
// }


// process align {
//   input:
  
//   output:
  
//   script:
//   """
//   echo Dummy
//   """
// }



workflow {

// inputs
  read_ch = channel.fromFilePairs( params.reads )
                   .map{ it -> [ it[1][0].parent, it[0], it[1][0], it[1][1] ] }

  seqs_list = params.seqs?.tokenize(',')
  seqs_ch = seqs_list ? channel.fromList( seqs_list ) : channel.empty()


// upstream
  merge_reads(read_ch)
  qc_post_merge(merge_reads.out)

  trim(merge_reads.out)
  qc_post_trim(trim.out)

  assemble(trim.out)

  map_contigs(trim.out.join(assemble.out.sub, by: [0,1]))
  sam_post_map_contigs(map_contigs.out)
  bcf_post_map_contigs(sam_post_map_contigs.out.bam.join(assemble.out.sub, by: [0,1]))

  blast(assemble.out.sub)


// downstream
  seqfile(seqs_ch)

  map_refs(trim.out.combine(seqfile.out))
  sam_post_map_refs(map_refs.out)
  bcf_post_map_refs( seqfile.out
    .cross(sam_post_map_refs.out.bam)
    .map{ zit -> [ zit[1][1], zit[1][2], zit[1][0], zit[1][3], zit[1][4], zit[0][1] ] } )

// to be addded : nodefile, align
}
