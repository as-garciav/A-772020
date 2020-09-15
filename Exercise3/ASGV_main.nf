#!/usr/bin/env nextflow
params.input_dir = false
params.output_dir = false
params.fasta_pattern = false


if (params.input_dir) {
  input_dir = params.input_dir - ~/\/$/
  output_dir = params.output_dir - ~/\/$/
  fasta_pattern = params.fasta_pattern
  fasta_files = input_dir + '/' + fasta_pattern
  Channel
    .fromPath(fasta_files)
    .ifEmpty { error "Cannot find any fastas matching: ${fasta_files}" }
    .set { fastas }
}

//seqsero
process serotyping {
   memory '2 GB'
   publishDir "${output_dir}/seqsero",
   mode:'copy', 
   saveAs: { file -> "SeqSero_result_${fasta}_dir"}
 

  input:
  file (fasta) from fastas

  output:
  file('SeqSero_result_*') 

  """
  SeqSero2_package.py -m k -p 2 -t 4 -i ${fasta} 
  """
}

//SISTR
process insilico_serotyping {
   memory '2 GB'
   publishDir "${output_dir}/sistr",
   mode:'copy', 
   saveAs: { file -> "SISTR_result_${fasta}_dir"}
 

  input:
  file (fasta) from fastas

  output:
  file('*sistr-output.tab') into tables

  """
  sistr --qc -vv --alleles-output allele-results.json --novel-alleles novel-alleles.fasta --cgmlst-profiles cmgmlst-profiles.csv -f tab -o sistr-output.tab ${fasta} 
  """
}

SeqSero2_tabs = Channel.fromPath("SeqSero_result_${fasta}_dir")

// Concatenate files
process cat_files {
	input:
	file query_file from SeqSero2_tabs
	file sistr_tab from tables
	
	output:
	file "antimicrobial_genes_serotyping.txt"
	
	"""
	cat *SeqSero_result_* > SeqSer2_results.txt
	cat sistr_tab > SISTR_results.txt
	// the idea would be to use awk to concatenate both results
  // awk   | query_file > antimicrobial_genes_serotyping.txt
	"""
	
}
