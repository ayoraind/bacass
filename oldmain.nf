#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// include non-process modules
include { help_message; version_message; complete_message; error_message; pipeline_start_message } from './modules/messages.nf'
include { default_params; check_params } from './modules/params_parser.nf'
include { help_or_version } from './modules/params_utilities.nf'

version = '1.0dev'

// setup default params
default_params = default_params()

// merge defaults with user params
merged_params = default_params + params

// help and version messages
help_or_version(merged_params, version)

final_params = check_params(merged_params)

// starting pipeline
pipeline_start_message(version, final_params)


// include processes
include { PROKKA; BAKTA; PANAROO_RUN; IQTREE } from './modules/processes.nf' addParams(final_params)

workflow {
	
	if (final_params.panaroo) {
	
	// check if user supplied the filepath to gff files
	
	if (final_params.gff_input) {
	gffinput_ch = Channel
			     .fromPath(final_params.gff_input)
        		     .collect()
        		     .ifEmpty { error "Cannot find any gff file matching: ${final_params.gff_input}" }
			     
	PANAROO_RUN(gffinput_ch)
	
	} else {
	
	// Message to user
	if (!final_params.bakta && !final_params.prokka) {
	exit("""
	ERROR!
	A major error has occured!
	
	==> User forgot to specify --prokka or --bakta arguments. Check nextflow run /path/to/main.nf --help for further details."
	
	""")
	
	} else {
	
	// if gff files are absent, run genome annotation afresh and run pangenome analyses afterwards, provided that the user supplied the filepath to genome assemblies
	
	if (final_params.assemblies && final_params.output_dir) {
        assemblies_ch = Channel
        			.fromPath(final_params.assemblies)
        			.map { file -> tuple(file.simpleName, file) }
        			.ifEmpty { error "Cannot find any assembly fasta file matching: ${final_params.assemblies}" }
				
	// if assemblies filepaths are specified, check if the argument '--prokka' or '--bakta' is provided.
        if (final_params.prokka) {
            PROKKA(assemblies_ch)
	    PANAROO_RUN(PROKKA.out.gff_ch.collect())
    }
    
    	if (final_params.bakta && !final_params.bakta_db) {
            exit("""
	ERROR!
	A major error has occured!
	
	==> User forgot to specify the --bakta_db argument. Check nextflow run /path/to/main.nf --help for further details."
	
	""")
    }
        if (final_params.bakta && final_params.bakta_db) {
            BAKTA(assemblies_ch, final_params.bakta_db)
	    PANAROO_RUN(BAKTA.out.gff_ch.collect())
    }
        
    }
	}
	
	 
	}
	
	} else {
	
	// Message to user
	if (!final_params.prokka && !final_params.bakta) {
	exit("""
	ERROR!
	A major error has occured!
	
	==> User forgot to specify --prokka or --bakta arguments. Check nextflow run /path/to/main.nf --help for further details."
	
	""")
	
	} else {
	
	// this means that only the genome annotation part would be run since the --panaroo option was not supplied
		// check if user supplied the file path to the genome assemblies
	
	if (final_params.assemblies && final_params.output_dir) {
        assemblies_ch = Channel
        			.fromPath(final_params.assemblies)
        			.map { file -> tuple(file.simpleName, file) }
        			.ifEmpty { error "Cannot find any assembly fasta file matching: ${final_params.assemblies}" }
				
	// if assemblies filepaths are specified, check if the argument '--prokka' or '--bakta' is given.
        if (final_params.prokka) {
            PROKKA(assemblies_ch)
	   
    }
    
   	 if (final_params.bakta && !final_params.bakta_db) {
            exit("""
	ERROR!
	A major error has occured!
	
	==> User forgot to specify the --bakta_db argument. Check nextflow run /path/to/main.nf --help for further details."
	
	""")
    }
    
        if (final_params.bakta && final_params.bakta_db) {
            BAKTA(assemblies_ch, final_params.bakta_db)
	    
    }
	}
	
	}
	
}

}

workflow.onComplete {
    complete_message(final_params, workflow, version)
}

workflow.onError {
    error_message(workflow)
}
