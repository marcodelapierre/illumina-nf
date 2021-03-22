#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.reads = 'R{1,2}.fastq.gz'
params.seqs = ''
params.contigs = ''
params.interleave = false
params.cascade = false
params.cascade_align_seqs = ''

params.min_len_contig = '1000'
params.evalue = '0.1'
params.outprefix = 'results_'
params.refdir = 'refseqs'


process merge_interl {
  tag "${dir}/${name}"
  publishDir "${dir}/${params.outprefix}${name}", mode: 'copy'
  stageInMode { ( workflow.profile == 'zeus' ) ? 'copy' : 'symlink' }

  input:
  tuple val(dir), val(name), path(read1), path(read2)

  output:
  tuple val(dir), val(name), path{ params.interleave ? 'interleaved.fastq.gz' : 'merged.fastq.gz' }

  script:
  if ( params.interleave )
    """
    reformat.sh \
      in1=$read1 in2=$read2 \
      out=interleaved.fastq.gz
    """
  else
    """
    bbmerge.sh \
      in1=$read1 in2=$read2 \
      out=merged.fastq.gz
    """
}


process qc_post_merge_interl {
  tag "${dir}/${name}"
  publishDir "${dir}/${params.outprefix}${name}", mode: 'copy'

  input:
  tuple val(dir), val(name), path(processed_fastq_gz)

  output:
  tuple val(dir), val(name), path{ params.interleave ? 'interleaved_fastqc.html' : 'merged_fastqc.html' }, path{ params.interleave ? 'interleaved_fastqc.zip' : 'merged_fastqc.zip' }

  script:
  """
  fastqc $processed_fastq_gz
  """
}


process trim {
  tag "${dir}/${name}"
  publishDir "${dir}/${params.outprefix}${name}", mode: 'copy'

  input:
  tuple val(dir), val(name), path(processed_fastq_gz)

  output:
  tuple val(dir), val(name), path('clean.fastq.gz')

  script:
  """
  if [ "${params.interleave}" == "true" ] ; then
    INTERL="interleaved=t"
  else
    INTERL=""
  fi
  bbduk.sh \
    in=$processed_fastq_gz \
    out=trimmed-partial.fastq.gz \
    \$INTERL ref=adapters ktrim=r k=27 hdist=2 edist=0 mink=4

  bbduk.sh \
    in=trimmed-partial.fastq.gz \
    out=clean.fastq.gz \
    \$INTERL ref=adapters ktrim=l k=27 hdist=2 edist=0 mink=4 \
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
  if [ "${params.interleave}" == "true" ] ; then
    INTERL="--12"
  else
    INTERL="-s"
  fi
  memory='${task.memory}'
  memory=\${memory% *}

  spades.py \
    \$INTERL clean.fastq.gz \
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
  if [ "${params.interleave}" == "true" ] ; then
    INTERL="interleaved=t"
  else
    INTERL=""
  fi
  bbmap.sh \
    in=clean.fastq.gz \
    ref=contigs_sub.fasta \
    out=mapped_contigs_sub_unsorted.sam \
    \$INTERL k=13 maxindel=16000 ambig=random \
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

  input:
  tuple val(order), val(seqid)

  output:
  tuple val(order), val(seqid), path('refseq.fasta')

  script:
  """
  seqid="${seqid}"
  blastdbcmd \
    -db ${params.blast_db} -entry \${seqid%_RC} \
    -line_length 60 \
    -out refseq.fasta

  sed -i '/^>/ s/ .*//g' refseq.fasta
  """
}


process sam_post_seqfile {
  tag "${seqid}"
  publishDir "${params.refdir}", mode: 'copy', saveAs: { filename -> "refseq_${seqid}.fasta" }

  input:
  tuple val(order), val(seqid), path('refseq.fasta')

  output:
  tuple val(order), val(seqid), path('refseq.fasta')

  script:
  """
  seqid="${seqid}"
  if [ "\${seqid: -3}" == "_RC" ] ; then
    samtools faidx \
      -i -o refseq_revcom.fasta \
      refseq.fasta \${seqid%_RC}

    mv refseq_revcom.fasta refseq.fasta

    sed -i '/^>/ s/ .*//g' refseq.fasta
  fi
  """
}


process map_refs {
  tag "${dir}/${name}/${seqid}"

  input:
  tuple val(dir), val(name), path('clean.fastq.gz'), val(seqid), path('refseq.fasta')

  output:
  tuple val(dir), val(name), val(seqid), path('mapped_refseq_unsorted.sam')

  script:
  """
  if [ "${params.interleave}" == "true" ] ; then
    INTERL="interleaved=t"
  else
    INTERL=""
  fi
  bbmap.sh \
    in=clean.fastq.gz \
    ref=refseq.fasta \
    out=mapped_refseq_unsorted.sam \
    \$INTERL k=13 maxindel=16000 ambig=random \
    threads=${task.cpus}
  """
}


// All cascade maps in one task, 
// because DSL2 does not allow feedback loop processes
process map_refs_cascade_bulk {
  tag "${dir}/${name}${params.hash_cascade}"
  publishDir "${dir}/${params.outprefix}${name}${params.hash_cascade}", mode: 'copy', pattern: 'labels_map_cascade.txt'

  input:
  tuple val(dir), val(name), path('clean.fastq.gz'), val(seqids), path('refseq*.fasta')

  output:
  tuple val(dir), val(name), val(seqids), path('mapped_refseq_unsorted*.sam'), emit: sam
  tuple val(dir), val(name), val(seqids), path('labels_map_cascade.txt'), emit: labels

  script:
  //new_seqids = seqids.collect{ params.hash_cascade + '_' + it }
  """
  if [ "${params.interleave}" == "true" ] ; then
    INTERL="interleaved=t"
  else
    INTERL=""
  fi
  echo $seqids >labels_map_cascade.txt
  num=\$( ls refseq*.fasta | wc -w )

  bbmap.sh \
    in=clean.fastq.gz \
    ref=refseq1.fasta \
    out=mapped_refseq_unsorted1.sam \
    outu=unmapped_refseq1.fastq.gz \
    path=ref1 \
    \$INTERL k=13 maxindel=16000 ambig=random \
    threads=${task.cpus}

  for i in \$( seq 2 \$num ) ; do
    bbmap.sh \
      in=unmapped_refseq\$((i-1)).fastq.gz \
      ref=refseq\${i}.fasta \
      out=mapped_refseq_unsorted\${i}.sam \
      outu=unmapped_refseq\${i}.fastq.gz \
      path=ref\${i} \
      \$INTERL k=13 maxindel=16000 ambig=random \
      threads=${task.cpus}
  done
  """
}


process sam_post_map_refs {
  tag "${dir}/${name}/${seqid}"
  publishDir "${dir}/${params.outprefix}${name}${params.hash_cascade}/${seqid}", mode: 'copy'
//  publishDir "${dir}/${params.outprefix}${name}", mode: 'copy', saveAs: { filename -> filename.replaceFirst(/_refseq/,"_refseq_$seqid") }

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
  tag "${dir}/${name}/${seqid}"
  publishDir "${dir}/${params.outprefix}${name}${params.hash_cascade}/${seqid}", mode: 'copy'

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


process contigfile {
  tag "${dir}/${name}_${contigid}"
  publishDir "${dir}/${params.outprefix}${name}", mode: 'copy', saveAs: { filename -> "consensus_contig_${contigid}.fasta" }

  input:
  tuple val(dir), val(name), path('consensus_contigs_sub.fasta'), val(contigid)
  
  output:
  tuple val(dir), val(name), val(contigid), path('consensus_contig.fasta')
  
  script:
  """
  contigid="${contigid}"
  idawk=\${contigid#NODE_}
  idawk=\${idawk%_RC}
  awk -F _ -v id=\$idawk '{ if(ok==1){if(\$1==">NODE"){exit}; print} ; if(ok!=1 && \$1==">NODE" && \$2==id){ok=1; print} }' consensus_contigs_sub.fasta >consensus_contig.fasta

  if [ "\${contigid: -3}" == "_RC" ] ; then
    samtools faidx \
       -i -o consensus_contig_revcom.fasta \
       consensus_contig.fasta \$(grep "^>\${contigid%_RC}_" consensus_contig.fasta | tr -d '>')
    mv consensus_contig_revcom.fasta consensus_contig.fasta
  fi
  """
}


process align {
  tag "${dir}/${name}/${hash_align}"
  publishDir "${dir}/${params.outprefix}${name}${params.hash_cascade}", mode: 'copy', saveAs: { filename -> file(filename).getSimpleName()+"_${hash_align}."+file(filename).getExtension() }

  input:
  tuple val(dir), val(name), val(seqids), path("consensus_refseq_*.fasta"), val(contigids), path("consensus_contig_*.fasta")
  
  output:
  tuple val(dir), val(name), path('labels_refseqs_contigs.txt'), path('aligned.fasta')
  
  script:
  hash_align = (seqids+contigids).toString().digest('SHA-1').substring(0,8)  // an alternative is .md5()
  """
  echo $seqids >labels_refseqs_contigs.txt
  echo $contigids >>labels_refseqs_contigs.txt
  cat consensus_refseq_*.fasta consensus_contig_*.fasta >input_align.fasta
  
  mafft-linsi \
    --thread ${task.cpus} \
    input_align.fasta >aligned.fasta
  """
}



workflow {

// inputs
  read_ch = channel.fromFilePairs( params.reads )
                   .map{ it -> [ it[1][0].parent, it[0], it[1][0], it[1][1] ] }

  seqs_list = params.seqs?.toUpperCase()
  seqs_list = seqs_list?.replaceAll(/\/RC/, "_RC")
  seqs_list = seqs_list?.tokenize(',')
  seqs_list_ord = []
  if ( seqs_list ) {
    for (i = 0 ; i < seqs_list.size() ; i++) {
      seqs_list_ord.push([ i, seqs_list[i] ])
    }
  }
  if ( params.cascade ) {
    seqs_ch = seqs_list_ord ? ( channel.fromList( seqs_list_ord ).map{ uit -> [ uit[0], uit[1]] } ) : channel.empty()
    params.hash_cascade = '/cascade_map_' + seqs_list?.toString().digest('SHA-1').substring(0,8)  // an alternative is .md5()
  } else {
    seqs_ch = seqs_list_ord ? ( channel.fromList( seqs_list_ord ).map{ uit -> [ '0', uit[1] ] } ) : channel.empty()
    params.hash_cascade = ''
  }

  contigs_list = params.contigs?.toUpperCase()
  contigs_list = contigs_list?.replaceAll(/\/RC/, "_RC")
  contigs_list = contigs_list?.tokenize(',')
  contigs_ch = contigs_list ? channel.fromList( contigs_list ) : channel.empty()


// upstream
  merge_interl(read_ch)
  qc_post_merge_interl(merge_interl.out)

  trim(merge_interl.out)
  qc_post_trim(trim.out)

  assemble(trim.out)

  map_contigs(trim.out.join(assemble.out.sub, by: [0,1]))
  sam_post_map_contigs(map_contigs.out)
  bcf_post_map_contigs( sam_post_map_contigs.out.bam
    .join(assemble.out.sub, by: [0,1]) )

  blast(assemble.out.sub)


// downstream
  seqfile(seqs_ch)
  sam_post_seqfile(seqfile.out)

  contigfile(bcf_post_map_contigs.out.cons.combine(contigs_ch))

  if ( params.cascade ) {
    map_refs_cascade_bulk( trim.out
      .combine( sam_post_seqfile.out
        .toSortedList({ a, b -> a[0] <=> b[0] })
        .transpose()
        .collate(3)
        .map{ vit -> [ vit[1], vit[2] ] } ) )
    sam_post_map_refs(map_refs_cascade_bulk.out.sam.transpose())
  } else {
    map_refs(trim.out.combine(sam_post_seqfile.out.map{ wit -> [ wit[1], wit[2] ] }))
    sam_post_map_refs(map_refs.out)
  }
  bcf_post_map_refs( sam_post_seqfile.out.map{ xit -> [ xit[1], xit[2] ] }
    .cross(sam_post_map_refs.out.bam)
    .map{ zit -> [ zit[1][1], zit[1][2], zit[1][0], zit[1][3], zit[1][4], zit[0][1] ] } )

  if ( params.cascade ) {
    cas_list = params.cascade_align_seqs?.toUpperCase()
    cas_list = cas_list?.replaceAll(/\/RC/, "_RC")
    cas_list = cas_list?.tokenize(',')
    cas_ch = cas_list ? channel.fromList( cas_list ) : channel.empty()

    align( cas_ch.cross( bcf_post_map_refs.out.cons
      .map { jit -> [ jit[2], jit[0], jit[1], jit[3] ] } )
      .map { kit -> [ kit[1][1], kit[1][2], kit[1][0], kit[1][3] ] }
      .groupTuple(by: [0,1])
      .join(contigfile.out.groupTuple(by: [0,1]), by: [0,1])
      .map{ yit -> [ yit[0], yit[1], yit[2].sort(), yit[3].sort(), yit[4].sort(), yit[5].sort() ] } )
  } else {
    align( bcf_post_map_refs.out.cons.groupTuple(by: [0,1])
      .join(contigfile.out.groupTuple(by: [0,1]), by: [0,1])
      .map{ yit -> [ yit[0], yit[1], yit[2].sort(), yit[3].sort(), yit[4].sort(), yit[5].sort() ] } )
  }

}
