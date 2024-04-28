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
include { FASTQC as FASTQC_RAW; FASTP; FASTQC as FASTQC_TRIM; NANOPLOT; PYCOQC; PORECHOP_PORECHOP; UNICYCLER; CANU; MINIMAP2_ALIGN as MINIMAP2_POLISH; MINIMAP2_ALIGN as MINIMAP2_CONSENSUS; MINIASM; RACON; DRAGONFLYE; MEDAKA; KRAKEN2_DB_PREPARATION; KRAKEN2; KRAKEN2_LONG; QUAST; CUSTOM_DUMPSOFTWAREVERSIONS } from './modules/processes.nf' addParams(final_params)

//
// SUBWORKFLOWS: Consisting of a mix of local and nf-core/modules
//
include { paramsSummaryMultiqc                  } from './subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                } from './subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                } from './subworkflows/local/utils_nfcore_bacass_pipeline'


//
// Function that parses fastp json output file to get total number of reads after trimming
//
import groovy.json.JsonSlurper

def getFastpReadsAfterFiltering(json_file) {
    def Map json = (Map) new JsonSlurper().parseText(json_file.text).get('summary')
    return json['after_filtering']['total_reads'].toLong()
}


workflow {

		ch_versions = Channel.empty()
    		ch_multiqc_files = Channel.empty()
    		
    		
		
	//	ch_samplesheet = Channel
		//			.fromPath(final_params.input, checkIfExists: true)
			//		.splitCsv(header:true, sep: '\t')
				//	.map { row -> [row.ID, row.R1, row.R2, row.LongFastQ, row.Fast5, row.GenomeSize] }
				
		
		ch_samplesheet = Channel
					.fromPath(final_params.input, checkIfExists: true)
					.splitCsv(header:true, sep: '\t')
					
		ch_adapter = Channel
				    .fromPath(final_params.adapter_file, checkIfExists: true)
				    
		// Check krakendb
			if(! params.skip_kraken2){
    				if(params.kraken2db){
        				kraken2db = file(params.kraken2db)
    				} else {
        				exit 1, "Missing Kraken2 DB arg"
    				}
			}


	//	ch_kraken2_db = Channel
		//		       .fromPath(final_params.kraken2db, checkIfExists: true)
					
		
		// channel for short reads
		ch_shortreads = ch_samplesheet
					     .map{  row ->
							def meta = row['ID']
							def fastq_1 = row['R1']
							def fastq_2 = row['R2']
						
						// Return a tuple with ID and R1, R2 if available, else only ID and R1
						// return fastq_2 ? [meta, fastq_1, fastq_2] : [id, fastq_1]
						
						// Combine R1 and R2 in a list; if only R1 is provided, create a list with one element
						def reads = fastq_2 ? [fastq_1, fastq_2] : [fastq_1]
						return [meta, reads]
						}
		
					
    		// channel for long reads
		ch_longreads = ch_samplesheet
					     .map{  row ->
							def meta = row['ID']
							def long_fastq = row['LongFastQ']
							return [meta, long_fastq]
					    
					     }
		
		// channel for fast5
		ch_fast5 = ch_samplesheet
					 .map{  row ->
							def meta = row['ID']
							def fast5 = row['Fast5']
							return [meta, fast5]
					    
					     }
		
	//	ch_shortreads.view()
	//	ch_longreads.view()
	//	ch_fast5.view()
		
              // QC
	      
	      if (final_params.assembly_type != 'long') {
	      
	      ch_fastqc_raw_html = Channel.empty()
              ch_fastqc_raw_zip  = Channel.empty()
	      
	      FASTQC_RAW(ch_shortreads)
	      
	      ch_fastqc_raw_html = FASTQC_RAW.out.html
       	      ch_fastqc_raw_zip  = FASTQC_RAW.out.zip
              ch_versions     = ch_versions.mix(FASTQC_RAW.out.versions.first())
	
		
	      // trimming step
	      ch_trim_json         = Channel.empty()
    	      ch_trim_html         = Channel.empty()
              ch_trim_log          = Channel.empty()
              ch_trim_reads_fail   = Channel.empty()
              ch_trim_reads_merged = Channel.empty()
              ch_fastqc_trim_html  = Channel.empty()
              ch_fastqc_trim_zip   = Channel.empty()
	      
	     	
	      
	      FASTP(ch_shortreads.combine(ch_adapter), final_params.save_trimmed_fail, final_params.save_merged)
	      
	       	ch_trim_reads        = FASTP.out.reads
        	ch_trim_json         = FASTP.out.json
        	ch_trim_html         = FASTP.out.html
        	ch_trim_log          = FASTP.out.log
        	ch_trim_reads_fail   = FASTP.out.reads_fail
        	ch_trim_reads_merged = FASTP.out.reads_merged
        	ch_versions       = ch_versions.mix(FASTP.out.versions.first())
		
	       //
               // Filter empty FastQ files after adapter trimming so FastQC doesn't fail
               //
            ch_trim_reads
            	.join(ch_trim_json)
            	.map {
                	meta, reads, json ->
                    	if (getFastpReadsAfterFiltering(json) > 0) {
                        	[ meta, reads ]
                    	}
            	}
            	.set { ch_trim_reads }
		
	   if (!params.val_skip_fastqc) {
            FASTQC_TRIM (
                ch_trim_reads
            )
            ch_fastqc_trim_html = FASTQC_TRIM.out.html
            ch_fastqc_trim_zip  = FASTQC_TRIM.out.zip
            ch_versions         = ch_versions.mix(FASTQC_TRIM.out.versions.first())
        }
	      }
	      
	      NANOPLOT (
        	ch_longreads
    		)
    		ch_nanoplot_txt_multiqc = NANOPLOT.out.txt
    		ch_versions = ch_versions.mix(NANOPLOT.out.versions)
		
		//
    		// MODULE: PYCOQC, quality check for nanopore reads and Quality/Length Plots
    		//
    		// TODO: Couldn't be tested. No configuration test available (lack of fast5 file or params.skip_pycoqc=false).
    	//	ch_pycoqc_multiqc = Channel.empty()
    	//	if ( !params.skip_pycoqc ) {
        	//	PYCOQC (
            	//	ch_fast5.dump(tag: 'fast5')
        	//	)
        	//	ch_pycoqc_multiqc = PYCOQC.out.json
        	//	ch_versions       = ch_versions.mix(PYCOQC.out.versions)
  //  }
  
  		//
   	 // MODULE: PORECHOP, quality check for nanopore reads and Quality/Length Plots
    		//
    		ch_porechop_log_multiqc = Channel.empty()
    		if ( params.assembly_type == 'hybrid' || params.assembly_type == 'long' && !('short' in params.assembly_type) ) {
        		PORECHOP_PORECHOP (
            		ch_longreads.dump(tag: 'longreads')
        		)
        		ch_porechop_log_multiqc = PORECHOP_PORECHOP.out.log
        		ch_versions = ch_versions.mix( PORECHOP_PORECHOP.out.versions )
    		}
		
		
		//
    		// Join channels for assemblers. As samples have the same meta data, we can simply use join() to merge the channels based on this. If we only have one of the channels we insert 'NAs' which are not used in the unicycler process then subsequently, in case of short or long read only assembly.
    		// Prepare channel for Kraken2
    		//
    		if(final_params.assembly_type == 'hybrid'){
        		ch_for_kraken2_short    = ch_trim_reads
        		ch_for_kraken2_long     = PORECHOP_PORECHOP.out.reads
        		ch_trim_reads
            			.dump(tag: 'fastp')
            			.join(PORECHOP_PORECHOP.out.reads)
            			.dump(tag: 'ch_for_assembly')
            			.set { ch_for_assembly }
    		} else if ( final_params.assembly_type == 'short' ) {
        		ch_for_kraken2_short    = ch_trim_reads
        		ch_for_kraken2_long     = Channel.empty()
        		ch_trim_reads
            			.dump(tag: 'fastp')
            			.map{ meta,reads -> tuple(meta,reads,[]) }
            			.dump(tag: 'ch_for_assembly')
            			.set { ch_for_assembly }
    		} else if ( final_params.assembly_type == 'long' ) {
        		ch_for_kraken2_short    = Channel.empty()
        		ch_for_kraken2_long     = PORECHOP_PORECHOP.out.reads
        		PORECHOP_PORECHOP.out.reads
            					.dump(tag: 'porechop')
            					.map{ meta,lr -> tuple(meta,[],lr) }
            					.dump(tag: 'ch_for_assembly')
            					.set { ch_for_assembly }
    }
    
    		//
    		// ASSEMBLY: Unicycler, Canu, Miniasm, Dragonflye
    		//
    		ch_assembly = Channel.empty()

    		//
    		// MODULE: Unicycler, genome assembly, nf-core module allows only short, long and hybrid assembly
    		//
    if ( final_params.assembler == 'unicycler' ) {
        UNICYCLER (
            ch_for_assembly
        )
        ch_assembly = ch_assembly.mix( UNICYCLER.out.scaffolds.dump(tag: 'unicycler') )
        ch_versions = ch_versions.mix( UNICYCLER.out.versions )
    }


    //
    // MODULE: Canu, genome assembly, long reads
    //
    if ( final_params.assembler == 'canu' ) {
        CANU (
            ch_for_assembly.map { meta, reads, lr -> tuple( meta, lr ) },
            params.canu_mode,
            ch_for_assembly.map { meta, reads, lr -> meta.genome_size }
        )
        ch_assembly = ch_assembly.mix( CANU.out.assembly.dump(tag: 'canu') )
        ch_versions = ch_versions.mix(CANU.out.versions)
    }

    //
    // MODULE: Miniasm, genome assembly, long reads
    if ( final_params.assembler == 'miniasm' ) {
        MINIMAP2_ALIGN (
            ch_for_assembly.map{ meta,sr,lr -> tuple(meta,lr) },
            [[:],[]],
            false,
            false,
            false
        )
        ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions)

        ch_for_assembly
            .join(MINIMAP2_ALIGN.out.paf)
            .map { meta, sr, lr, paf-> tuple(meta, lr, paf) }
            .set { ch_for_miniasm }

        MINIASM (
            ch_for_miniasm
        )
        ch_versions = ch_versions.mix(MINIASM.out.versions)

        MINIMAP2_CONSENSUS (
            ch_for_assembly.map{ meta,sr,lr -> tuple(meta,lr) },
            MINIASM.out.assembly,
            false,
            false,
            false
        )
        ch_versions = ch_versions.mix(MINIMAP2_CONSENSUS.out.versions)

        ch_for_assembly
            .join(MINIASM.out.assembly)
            .join(MINIMAP2_CONSENSUS.out.paf)
            .map { meta, sr, lr, assembly, paf -> tuple(meta, lr, assembly, paf) }
            .set{ ch_for_racon }

        RACON (
            ch_for_racon
        )
        ch_assembly = ch_assembly.mix( RACON.out.improved_assembly.dump(tag: 'miniasm') )
        ch_versions = ch_versions.mix( RACON.out.versions )
    }

    //
    // MODULE: Dragonflye, genome assembly of long reads. Moreover, it provides the option for polishing the draft genome using short reads when both short and long reads are available.
    //
    if( final_params.assembler == 'dragonflye' ){
        DRAGONFLYE(
            ch_for_assembly
        )
        ch_assembly = ch_assembly.mix( DRAGONFLYE.out.contigs.dump(tag: 'dragonflye') )
        ch_versions = ch_versions.mix( DRAGONFLYE.out.versions )
    }
    
     //
    // MODULE: Nanopolish, polishes assembly using FAST5 files - should take either miniasm, canu, or unicycler consensus sequence
    //
    if ( !params.skip_polish && params.assembly_type == 'long' && params.polish_method != 'medaka' ) {
        ch_for_assembly
            .join( ch_assembly )
            .set { ch_for_polish }

        MINIMAP2_POLISH (
            ch_for_polish.map { meta, sr, lr, fasta -> tuple(meta, lr)  },
            ch_for_polish.map { meta, sr, lr, fasta -> fasta  },
            true,
            false,
            false
        )
        ch_versions = ch_versions.mix(MINIMAP2_POLISH.out.versions)

        SAMTOOLS_INDEX (
            MINIMAP2_POLISH.out.bam.dump(tag: 'samtools_sort')
        )
        ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

      //  ch_for_polish    // tuple val(meta), val(reads), file(longreads), file(assembly)
        //    .join( MINIMAP2_POLISH.out.bam )    // tuple val(meta), file(bam)
         //   .join( SAMTOOLS_INDEX.out.bai )     // tuple  val(meta), file(bai)
          //  .join( ch_fast5 )             // tuple val(meta), file(fast5)
          //  .set { ch_for_nanopolish }          // tuple val(meta), val(reads), file(longreads), file(assembly), file(bam), file(bai), file(fast5)

        // TODO: 'nanopolish index' couldn't be tested. No fast5 provided in test datasets.
       // NANOPOLISH (
         //   ch_for_nanopolish.dump(tag: 'into_nanopolish')
       // )
       // ch_versions = ch_versions.mix(NANOPOLISH.out.versions)
    }
    
    
     //
    // MODULE: Medaka, polishes assembly - should take either miniasm, canu, or unicycler consensus sequence
    //
    if ( !params.skip_polish && params.assembly_type == 'long' && params.polish_method == 'medaka' ) {
        ch_for_assembly
            .join( ch_assembly )
            .map { meta, sr, lr, assembly -> tuple(meta, lr, assembly) }
            .set { ch_for_medaka }

        MEDAKA ( ch_for_medaka.dump(tag: 'into_medaka') )
        ch_versions = ch_versions.mix(MEDAKA.out.versions)
    }

   
   //
    // MODULE: Kraken2, QC for sample purity
    //
    ch_kraken_short_multiqc = Channel.empty()
    ch_kraken_long_multiqc  = Channel.empty()
    if ( !params.skip_kraken2 ) {
        KRAKEN2_DB_PREPARATION (
            kraken2db
        )
        ch_versions = ch_versions.mix(KRAKEN2_DB_PREPARATION.out.versions)
        KRAKEN2 (
            ch_for_kraken2_short.dump(tag: 'kraken2_short'),
            KRAKEN2_DB_PREPARATION.out.db.map { info, db -> db }.dump(tag: 'kraken2_db_preparation'),
            false,
            false
        )
        ch_kraken_short_multiqc = KRAKEN2.out.report
        ch_versions = ch_versions.mix(KRAKEN2.out.versions)

        KRAKEN2_LONG (
            ch_for_kraken2_long
              //  .map { meta, reads ->
                //    info = [:]
                 //   info.id = meta
                 //   params.single_end = true
                 //   [ info, reads ]
               // }
                .dump(tag: 'kraken2_long'),
            KRAKEN2_DB_PREPARATION.out.db.map { ch_for_kraken2_long, db -> db }.dump(tag: 'kraken2_db_preparation'),
            false,
            false
        )
        ch_kraken_long_multiqc = KRAKEN2_LONG.out.report
        ch_versions = ch_versions.mix(KRAKEN2_LONG.out.versions)
    }
    
     //
    // MODULE: QUAST, assembly QC
    //
  //  ch_assembly
    //    .collect{ it[1] }
     //   .map { consensus_collect -> tuple("report", consensus_collect) }
     //   .set { ch_to_quast }

    QUAST (
        ch_assembly,
        [[:],[]],
        [[:],[]]
    )
    ch_quast_multiqc = QUAST.out.tsv
    ch_versions      = ch_versions.mix(QUAST.out.versions)
    
    //
    // Collate and save software versions
    //
    CUSTOM_DUMPSOFTWAREVERSIONS (
                                ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )


}


workflow.onComplete {
    complete_message(final_params, workflow, version)
}

workflow.onError {
    error_message(workflow)
}
