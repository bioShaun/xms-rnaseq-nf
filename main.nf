#!/usr/bin/env nextflow

def helpMessage() {
    log.info """

    A Pipeline to run rnaseq analysis

    Usage:

    ==================================================================

    nextflow run xms-rnaseq-nf --fasta transcriptome.fa

    ==================================================================

    Mandatory arguments:
      --reads                       Path to input data (must be surrounded with quotes)
      --outdir                      Path to output directory
    
    Options:
      --singleEnd                   Specifies that the input is single end reads
      --gene2tr                     Gene_id vs Transcript_id map file

    References
      --fasta                       Path to Fasta reference
      --kallisto_index              Path to Kallisto index

    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

 // Show help emssage
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

params.outdir = false
params.singleEnd = false
params.gene2tr = false
params.kallisto_index = false

gene2tr_file = file(params.gene2tr)

kallisto_to_table = "${baseDir}/script/kallisto_to_table.R"

Channel
    .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\nIf this is single-end data, please specify --singleEnd on the command line." }
    .set { raw_reads }

if (params.kallisto_index) {
    kallisto_idx = file(params.kallisto_index)
} else if (params.fasta) {
    transcriptome_file = file(params.fasta)
    kallisto_idx = file("${params.fasta}.kallisto_idx")
} else {
    exit 1, "No reference fasta/index specified!"
}

if (kallisto_idx.exists()) {
    kallisto_idx_file = kallisto_idx
} else {
    fasta_path = transcriptome_file.getParent() 
    process kallisto_index {
        tag "KALLISTO INDEX on ${transcriptome_file.getName()}"

        publishDir fasta_path, mode: 'copy'

        input:
        file transcriptome from transcriptome_file

        output:
        file "${transcriptome}.kallisto_idx" into kallisto_idx_file

        script:
        """
        kallisto index -i "${transcriptome}.kallisto_idx" ${transcriptome}
        """
    }
}



process kallisto_quant {
    tag "KALLISTO QUANT on $sample_id"

    publishDir "${params.outdir}/kallisto/", mode: 'copy'

    input:
    set val(sample_id), file(reads) from raw_reads
    file kallisto_idx from kallisto_idx_file

    output:
    file sample_id into kallisto_out

    cpus 4

    script:
    """
    kallisto quant \
        --threads ${task.cpus} \
        -i ${kallisto_idx} \
        --output-dir=${sample_id} \
        ${reads}
    """
}

process merge_kallisto_count {
    tag "MERGE KALLISTO QUANT"

    publishDir "${params.outdir}/expression", mode: 'copy'

    input:
    file ('kallisto/*') from kallisto_out.collect()
    file gene2tr from gene2tr_file

    output:
    file "*tpm.txt" into tpm_file
    file "*count.txt" into count_file    

    script:
    """
    Rscript ${kallisto_to_table} \
        --kallisto_dir ./kallisto \
        --out_dir . \
        --gene2tr ${gene2tr}
    """
}