#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


process prep_seqs {

  input:
  tuple val(num), val(id)

  output:
  tuple val(num), val(id), val(what)

  script:
  what = "seq_$id"
  """
  if [ "$id" == "id2"] ; then sleep 2 ; fi
  """

}


// testing purposes
process echo_seq {

  input:
  val(bulk_seqs)

  output:
  val(bulk_seqs)

  exec:
//  active0 = bulk_seqs.pop()
//  active0 = bulk_seqs.pop()
  active = bulk_seqs.pop()
  id = active[1]
  seq = active[2]
  println("Out ID: "+id+" Seq: "+seq)
  assert "Out ID: "+id+" Seq: "+seq == 'Out ID: id1 Seq: seq_id1'
  //println(bulk_seqs)
}


process run_refs {

  input:
  tuple val(read), val(id), val(seq)

  output:
  tuple val('whatever'), path('good'), emit: out
  tuple path('minus'), val(seq), emit: loop
  
  script:
  """
  echo "this is $read after $seq" >>good
  echo "this is $read minus $seq" >>minus
  """
}


workflow {

  reads = channel.fromList(['sam1.fq', 'sam2.fq'])

  list = ['id1', 'id2', 'id3']
  list2 = []
  for (i = 0 ; i < list.size() ; i++) {
      list2.push([ i, list[i] ])
  }
  seqs = channel.fromList(list2).map{ zit -> [ zit[0], zit[1]] }

  prep_seqs(seqs)
//  prep_seqs.out.view()


  bulk = prep_seqs.out.toSortedList({ a, b -> a[0] <=> b[0] })
//  bulk.view()

// Nice test for reduce operator
  //bulk_seqs = prep_seqs.out
  //  .reduce( [ ''.toList(), ''.toList(), ''.toList() ] ){ a,b -> return [ a[0]+b[0], a[1]+b[1], a[2]+b[2] ] }


//
// Testing feedback loop -- apparently missing Nextflow capability
//
  echo_seq(bulk)  // one loop iteration only
//  feedback = channel.of()
//  input = bulk.mix(feedback)
//  input | echo_seq | tap(feedback)

// Ideally aiming to test prototype process
  //run_refs(reads.combine(bulk_seqs))

}
