include { check_mandatory_parameter; check_optional_parameters; check_parameter_value } from './params_utilities.nf'
includeConfig './nextflow.config'

def default_params(){
    /***************** Setup inputs and channels ************************/
    def params = [:] as nextflow.script.ScriptBinding$ParamsMap
    // Defaults for configurable variables
    
    params {
     
     // Input options
    input                           = null

    // QC and trimming options
    save_trimmed_fail               = false
    save_merged                     = false

    // Contamination_screening
    kraken2db                       = ""

    // Assembly parameters
    assembler                       = 'unicycler'   // Allowed: ['unicycler', 'canu', 'miniasm', 'dragonflye']
    assembly_type                   = 'short'       // Allowed: ['short', 'long', 'hybrid'] (hybrid works only with Unicycler)
    unicycler_args                  = ""
    canu_mode                       = '-nanopore'   // Allowed: ['-pacbio', '-nanopore', '-pacbio-hifi']
    canu_args                       = ''            // Default no extra options, can be adjusted by the user
    dragonflye_args                 = ''

    // Assembly polishing
    polish_method                   = 'medaka'      // Allowed: ['medaka', 'nanopolish']

    
    // Skipping options
    skip_fastqc                     = false
    skip_fastp                      = false
    skip_kraken2                    = false
    skip_pycoqc                     = false
    skip_polish                     = false
    skip_multiqc                    = false

    // MultiQC options
    multiqc_config                  = null
    multiqc_title                   = null
    multiqc_logo                    = null
    max_multiqc_email_size          = '25.MB'
    multiqc_methods_description     = null

    // Boilerplate options
    outdir                          = null
    publish_dir_mode                = 'copy'
    email                           = null
    email_on_fail                   = null
    plaintext_email                 = false
    monochrome_logs                 = false
    hook_url                        = null
    help                            = false
    validate_params                 = true
    schema_ignore_params            = 'modules,igenomes_base'
    version                         = false

    // Config options
    config_profile_name             = null
    config_profile_description      = null
    // custom_config_version           = 'master'
    custom_config_base              = "$projectDir/configs"
    config_profile_contact          = null
    config_profile_url              = null

    // Max resource options
    // Defaults only, expecting to be overwritten
    max_memory                      = '128.GB'
    max_cpus                        = 16
    max_time                        = '240.h'

    // Schema validation default options
    validationFailUnrecognisedParams= false
    validationLenientMode           = false
    validationSchemaIgnoreParams    = 'genomes'
    validationShowHiddenParams      = false
    validate_params                 = true
    
    }
    
    return params
}

def check_params(Map params) { 
    final_params = params
     
    // set up input file
    final_params.input = check_mandatory_parameter(params, 'input') - ~/\/$/
    
    // set up output directory
    final_params.outdir = check_mandatory_parameter(params, 'outdir') - ~/\/$/
    
    final_params.kraken2db = if(! final_params.skip_kraken2){
    							check_mandatory_parameter(params, 'kraken2db') - ~/\/$/
    								}
              
    return final_params
}


// Load modules.config for DSL2 module specific options
includeConfig '$projectDir/conf/modules.config'

// Function to ensure that resource requirements don't go beyond
// a maximum limit
def check_max(obj, type) {
    if (type == 'memory') {
        try {
            if (obj.compareTo(params.max_memory as nextflow.util.MemoryUnit) == 1)
                return params.max_memory as nextflow.util.MemoryUnit
            else
                return obj
        } catch (all) {
            println "   ### ERROR ###   Max memory '${params.max_memory}' is not valid! Using default value: $obj"
            return obj
        }
    } else if (type == 'time') {
        try {
            if (obj.compareTo(params.max_time as nextflow.util.Duration) == 1)
                return params.max_time as nextflow.util.Duration
            else
                return obj
        } catch (all) {
            println "   ### ERROR ###   Max time '${params.max_time}' is not valid! Using default value: $obj"
            return obj
        }
    } else if (type == 'cpus') {
        try {
            return Math.min( obj, params.max_cpus as int )
        } catch (all) {
            println "   ### ERROR ###   Max cpus '${params.max_cpus}' is not valid! Using default value: $obj"
            return obj
        }
    }
}
