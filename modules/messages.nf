def help_message() {
  log.info """
        Usage:
        The typical command for running the pipeline is as follows:
        nextflow run main.nf --input "PathToCSVfile" --outdir "PathToOutputDir"
	
	or for hybrid assembly (short read assembly, then long reads to cover gap)
	
	nextflow run main.nf --input "PathToCSVfile" --outdir "PathToOutputDir" --assembler 'unicycler' --assembler_type 'hybrid'
	
	
	or for hybrid assembly (long read assembly, then short read polishing)
	
	nextflow run main.nf --gff_input "PathToGFFfiles" --output_dir "PathToOutputDir" --assembler 'dragonflye' --assembler_type 'long'
	
	
	or for short-read assembly alone
	
	nextflow run main.nf --input "PathToCSVfile" --outdir "PathToOutputDir" --assembler 'unicycler' --assembler_type 'short'
	
	
	or for long-read assembly alone
	
	nextflow run main.nf --input "PathToCSVfile" --outdir "PathToOutputDir" --assembler 'dragonflye' --assembler_type 'long' --skip_polish
	

        Mandatory arguments:
         --input                   	Path to sample sheet, either tab-separated (.tsv), comma-separated (.csv), or in YAML format (.yml/.yaml), that points to compressed fastq files.\n\nThe sample sheet must have six tab-separated columns/entries with the following headers: \n- `ID` (required): Unique sample IDs, must start with a letter, and can only contain letters, numbers or underscores\n- `R1` (optional): Paths to (forward) reads zipped FastQ files\n- `R2` (optional): Paths to reverse reads zipped FastQ files, required if the data is paired-end\n- `LongFastQ` (optional): Paths to long reads zipped FastQ files\n- `Fast5` (optional): Paths to the directory containing FAST5 files\n- `GenomeSize` (optional): A number (including decimals) ending with 'm', representing genome size.\n\n Please be aware that files will be required based on the chosen assembly type specified with the '--assembly_type' option, which can be set to one of the following values: ['short', 'long', 'hybrid'].`
         --outdir                   	Output directory to place output reads 
	 --kraken2db			Path to kraken2 database directory
	          
        Optional arguments:
	 --assembler			Allowed: ['unicycler', 'canu', 'miniasm', 'dragonflye']. Default: 'unicycler' 
	 --assembly_type		Allowed: ['short', 'long', 'hybrid'] (hybrid works only with Unicycler). Default: 'short'
	 --polish_method		Allowed: ['medaka', 'nanopolish']. Default: 'medaka' 
	 --canu_mode                    Allowed: ['-pacbio', '-nanopore', '-pacbio-hifi']. Default: '-nanopore' 
	 --iqtree			if interested in building a maximum-likelihood phylogenetic tree, this option has to be supplied.
         --help                         This usage statement.
         --version                      Version statement
        """
}



def version_message(String version) {
      println(
            """
            ============================================================================
             BACASS GENOME ASSEMBLY PIPELINE: TAPIR Pipeline version ${version}
            ============================================================================
            """.stripIndent()
        )

}

def pipeline_start_message(String version, Map params){
    log.info "================================================================================"
    log.info "      BACASS GENOME ASSEMBLY PIPELINE: TAPIR Pipeline version ${version}            "
    log.info "================================================================================"
    log.info "Running version            : ${version}"
    log.info "Input file                 : ${params.input}"
    log.info ""
    log.info "-------------------------- Other parameters --------------------------"
    params.sort{ it.key }.each{ k, v ->
        if (v){
            log.info "${k}: ${v}"
        }
    }
    log.info "================================================================================"
    log.info "Outputs written to path '${params.outdir}'"
    log.info "================================================================================"

    log.info ""
}


def complete_message(Map params, nextflow.script.WorkflowMetadata workflow, String version){
    // Display complete message
    log.info ""
    log.info "Ran the workflow: ${workflow.scriptName} ${version}"
    log.info "Command line    : ${workflow.commandLine}"
    log.info "Completed at    : ${workflow.complete}"
    log.info "Duration        : ${workflow.duration}"
    log.info "Success         : ${workflow.success}"
    log.info "Work directory  : ${workflow.workDir}"
    log.info "Exit status     : ${workflow.exitStatus}"
    log.info ""
}

def error_message(nextflow.script.WorkflowMetadata workflow){
    // Display error message
    log.info ""
    log.info "Workflow execution stopped with the following message:"
    log.info "  " + workflow.errorMessage
}
