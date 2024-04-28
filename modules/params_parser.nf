include { check_mandatory_parameter; check_optional_parameters; check_parameter_value } from './params_utilities.nf'
// includeConfig './nextflow.config'

def default_params(){
    /***************** Setup inputs and channels ************************/
    def params = [:] as nextflow.script.ScriptBinding$ParamsMap
    // Defaults for configurable variables
    
    
    // Input options
    params.input                    		= null

    // QC and trimming options
    params.save_trimmed_fail        		= false
    params.save_merged              		= false
    params.single_end				= false
    params.val_skip_fastqc			= false

    // Contamination_screening
    params.kraken2db                		= ""
    params.adapter_file				= ""

    // Assembly parameters
    params.assembler                		= 'unicycler'   // Allowed: ['unicycler', 'canu', 'miniasm', 'dragonflye']
    params.assembly_type            		= 'short'       // Allowed: ['short', 'long', 'hybrid'] (hybrid works only with Unicycler)
    params.unicycler_args           		= ""
    params.canu_mode                		= '-nanopore'   // Allowed: ['-pacbio', '-nanopore', '-pacbio-hifi']
    params.canu_args                		= ''            // Default no extra options, can be adjusted by the user
    params.dragonflye_args          		= ''

    // Assembly polishing
    params.polish_method            		= 'medaka'      // Allowed: ['medaka', 'nanopolish']

    
    // Skipping options
    params.skip_fastqc              		= false
    params.skip_fastp               		= false
    params.skip_kraken2             		= false
    params.skip_pycoqc              		= false
    params.skip_polish             	 	= false
    params.skip_multiqc             		= false

    // MultiQC options
    params.multiqc_config           		= null
    params.multiqc_title            		= null
    params.multiqc_logo             		= null
    params.max_multiqc_email_size   		= '25.MB'
    params.multiqc_methods_description     	= null

    // Boilerplate options
    params.outdir                          	= null
    params.publish_dir_mode                	= 'copy'
    params.email                           	= null
    params.email_on_fail                   	= null
    params.plaintext_email                 	= false
    params.monochrome_logs                 	= false
    params.hook_url                        	= null
    params.help                            	= false
    params.validate_params                 	= true
    params.schema_ignore_params            	= 'modules,igenomes_base'
    params.version                         	= false

    // Config options
    params.config_profile_name             	= null
    params.config_profile_description      	= null
    // custom_config_version           = 'master'
    params.custom_config_base              	= "$projectDir/configs"
    params.config_profile_contact          	= null
    params.config_profile_url              	= null

    // Max resource options
    // Defaults only, expecting to be overwritten
    params.max_memory                      	= '128.GB'
    params.max_cpus                        	= 16
    params.max_time                        	= '240.h'

    // Schema validation default options
    params.validationFailUnrecognisedParams	= false
    params.validationLenientMode           	= false
    params.validationSchemaIgnoreParams    	= 'genomes'
    params.validationShowHiddenParams      	= false
    
    return params
}

def check_params(Map params) { 
    final_params = params
     
    // set up input file
    final_params.input = check_mandatory_parameter(params, 'input') - ~/\/$/
    
    // set up output directory
    final_params.outdir = check_mandatory_parameter(params, 'outdir') - ~/\/$/
    
                 
    return final_params
}
