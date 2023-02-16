#!/usr/bin/env nextflow
nextflow.enable.dsl=1
// If the user uses the --help flag, print the help text below
params.help = false

// Function which prints help message text
def helpMessage() {
    log.info"""
    Run customized DADA2 pipeline to process MiSeq reads for strain purity check

    Required Arguments:
      --project       string      MiSeq project name
      --config        string      input parameter file
      --input_path    string      input:  s3 path
      --output_path   string      output: s3 path
      --db            string      db:  silva efs path

    Options:
      -profile        docker run locally

    """.stripIndent()
}

// Show help message if the user specifies the --help flag at runtime
if (params.help){
    // Invoke the function above which prints the help message
    helpMessage()
    // Exit out and do not run anything else
    exit 0
}

/*
 * Defines the pipeline inputs parameters (giving a default value for each for them)
 * Each of the following parameters can be specified as command line options
 */

def output_path = "${params.output_path}"
def input_path = "${params.input_path}"
def config_file = "${params.config}"

/*
 * Run the pipeline
 * 16 v4 region and Genus level only
 */
process run_dada2 {

    container params.container
    cpus 8
    memory 32.GB
    publishDir "${output_path}", mode:'copy'

    input:

    output:

    script:
    """
    export S3INPUTPATH="${input_path}/fastqs"
    export S3OUTPUTPATH="${output_path}/${params.project}"
    export CONFIG="${config_file}"
    export DB="${params.db}"
    16s_wrapper.sh
    """
}
