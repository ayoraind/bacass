params.gff_input = "/MIGE/04_PROJECTS/STRAIN_RESOLUTION_AND_TRANSMISSION_DYNAMICS/R10_genome_annotation/*_bakta/*.gff3"
params.output_dir = "test_out"

process PANAROO_RUN {
    tag "pangenome analyses"
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    path(gff)

    output:
    path("results/*")                                      		, emit: results_ch
    path("results/core_gene_alignment.aln"), optional: true		, emit: aln_ch
    path "versions.yml"                                                 , emit: versions_ch

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    
    """
    panaroo \\
        $args \\
        -t $task.cpus \\
        -o results \\
	-a core \\
	--clean-mode strict \\
	--remove-invalid-genes \\
        -i $gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        panaroo: \$(echo \$(panaroo --version 2>&1) | sed 's/^.*panaroo //' ))
    END_VERSIONS
    """
}

workflow {

	gffinput_ch = Channel
                             .fromPath(params.gff_input)
                             .collect()
                             .ifEmpty { error "Cannot find any gff file matching: ${params.gff_input}" }

     // gffinput_ch.view()  
       PANAROO_RUN(gffinput_ch)
}
