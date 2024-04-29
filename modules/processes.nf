process FASTQC {
    tag "$meta"
    label 'process_medium'
    cache 'lenient'
    
    conda "${projectDir}/conda_environments/fastp.yml"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip") , emit: zip
    path  "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta}"
    // Make list of old name and new name pairs to use for renaming in the bash while loop
    def old_new_pairs = reads instanceof Path || reads.size() == 1 ? [[ reads, "${prefix}.${reads.extension}" ]] : reads.withIndex().collect { entry, index -> [ entry, "${prefix}_${index + 1}.${entry.extension}" ] }
    def rename_to = old_new_pairs*.join(' ').join(' ')
    def renamed_files = old_new_pairs.collect{ old_name, new_name -> new_name }.join(' ')
    """
    printf "%s %s\\n" $rename_to | while read old_name new_name; do
        [ -f "\${new_name}" ] || ln -s \$old_name \$new_name
    done

    fastqc \\
        $args \\
        --threads $task.cpus \\
        $renamed_files

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$( fastqc --version | sed '/FastQC v/!d; s/.*v//' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta}"
    """
    touch ${prefix}.html
    touch ${prefix}.zip

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$( fastqc --version | sed '/FastQC v/!d; s/.*v//' )
    END_VERSIONS
    """
}

process FASTP {
    tag "$meta"
    label 'process_medium'
    cache 'lenient'
    
    conda "${projectDir}/conda_environments/fastp.yml"
    
    input:
    tuple val(meta), path(reads), path(adapter_fasta)
    val   save_trimmed_fail
    val   save_merged

    output:
    tuple val(meta), path('*.fastp.fastq.gz') , optional:true, emit: reads
    tuple val(meta), path('*.json')           , emit: json
    tuple val(meta), path('*.html')           , emit: html
    tuple val(meta), path('*.log')            , emit: log
    path "versions.yml"                       , emit: versions
    tuple val(meta), path('*.fail.fastq.gz')  , optional:true, emit: reads_fail
    tuple val(meta), path('*.merged.fastq.gz'), optional:true, emit: reads_merged

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta}"
    def adapter_list = adapter_fasta ? "--adapter_fasta ${adapter_fasta}" : ""
    def fail_fastq = save_trimmed_fail && params.single_end ? "--failed_out ${prefix}.fail.fastq.gz" : save_trimmed_fail && !params.single_end ? "--failed_out ${prefix}.paired.fail.fastq.gz --unpaired1 ${prefix}_1.fail.fastq.gz --unpaired2 ${prefix}_2.fail.fastq.gz" : ''
    // Added soft-links to original fastqs for consistent naming in MultiQC
    // Use single ended for interleaved. Add --interleaved_in in config.
    if ( task.ext.args?.contains('--interleaved_in') ) {
        """
        [ ! -f  ${prefix}.fastq.gz ] && ln -sf $reads ${prefix}.fastq.gz

        fastp \\
            --stdout \\
            --in1 ${prefix}.fastq.gz \\
            --thread $task.cpus \\
            --json ${prefix}.fastp.json \\
            --html ${prefix}.fastp.html \\
            $adapter_list \\
            $fail_fastq \\
            $args \\
            2> >(tee ${prefix}.fastp.log >&2) \\
        | gzip -c > ${prefix}.fastp.fastq.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
        END_VERSIONS
        """
    } else if (params.single_end) {
        """
        [ ! -f  ${prefix}.fastq.gz ] && ln -sf $reads ${prefix}.fastq.gz

        fastp \\
            --in1 ${prefix}.fastq.gz \\
            --out1  ${prefix}.fastp.fastq.gz \\
            --thread $task.cpus \\
            --json ${prefix}.fastp.json \\
            --html ${prefix}.fastp.html \\
            $adapter_list \\
            $fail_fastq \\
            $args \\
            2> >(tee ${prefix}.fastp.log >&2)

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
        END_VERSIONS
        """
    } else {
        def merge_fastq = save_merged ? "-m --merged_out ${prefix}.merged.fastq.gz" : ''
        """
        [ ! -f  ${prefix}_1.fastq.gz ] && ln -sf ${reads[0]} ${prefix}_1.fastq.gz
        [ ! -f  ${prefix}_2.fastq.gz ] && ln -sf ${reads[1]} ${prefix}_2.fastq.gz
        fastp \\
            --in1 ${prefix}_1.fastq.gz \\
            --in2 ${prefix}_2.fastq.gz \\
            --out1 ${prefix}_1.fastp.fastq.gz \\
            --out2 ${prefix}_2.fastp.fastq.gz \\
            --json ${prefix}.fastp.json \\
            --html ${prefix}.fastp.html \\
            $adapter_list \\
            $fail_fastq \\
            $merge_fastq \\
            --thread $task.cpus \\
            --detect_adapter_for_pe \\
            $args \\
            2> >(tee ${prefix}.fastp.log >&2)

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
        END_VERSIONS
        """
    }
    
     stub:
    def prefix = task.ext.prefix ?: "${meta}"

    """
    touch $touch_reads
    touch "${prefix}.fastp.json"
    touch "${prefix}.fastp.html"
    touch "${prefix}.fastp.log"
    $touch_merged

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
    END_VERSIONS
    """
}


process NANOPLOT {
    tag "$meta"
    label 'process_low'
    cache 'lenient'
    
    conda "${projectDir}/conda_environments/nanoplot.yml"

    input:
    tuple val(meta), path(ontfile)

    output:
    tuple val(meta), path("*.html")                , emit: html
    tuple val(meta), path("*.png") , optional: true, emit: png
    tuple val(meta), path("*.txt")                 , emit: txt
    tuple val(meta), path("*.log")                 , emit: log
    path  "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta}"
    def input_file = ("$ontfile".endsWith(".fastq.gz") || "$ontfile".endsWith(".fq.gz")) ? "--fastq ${ontfile}" :  ("$ontfile".endsWith(".txt")) ? "--summary ${ontfile}" : ''
    """
    NanoPlot \\
        $args \\
        -t $task.cpus \\
        $input_file

    mv NanoStats.txt ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanoplot: \$(echo \$(NanoPlot --version 2>&1) | sed 's/^.*NanoPlot //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch LengthvsQualityScatterPlot_dot.html
    touch LengthvsQualityScatterPlot_kde.html
    touch NanoPlot-report.html
    touch NanoPlot_20240301_1130.log
    touch NanoStats.txt
    touch Non_weightedHistogramReadlength.html
    touch Non_weightedLogTransformed_HistogramReadlength.html
    touch WeightedHistogramReadlength.html
    touch WeightedLogTransformed_HistogramReadlength.html
    touch Yield_By_Length.html


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanoplot: \$(echo \$(NanoPlot --version 2>&1) | sed 's/^.*NanoPlot //; s/ .*\$//')
    END_VERSIONS
    """
}


process PYCOQC {
    label 'process_medium'
    cache 'lenient'
    
    conda "${projectDir}/conda_environments/pycoqc.yml"
    
    input:
    tuple val(meta), path(fast5)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.json"), emit: json
    path  "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta}"
    def run_summary = file("${fast5}/sequencing_summary.txt").exists() ? "cp ${fast5}/sequencing_summary.txt ./sequencing_summary.txt" : "Fast5_to_seq_summary -f $fast5 -t ${task.cpus} -s './sequencing_summary.txt' --verbose_level 2"
    def barcode_me  = file("${fast5}/barcoding_sequencing.txt").exists() ? "-b ${fast5}/barcoding_sequencing.txt" : ''

    """
    $run_summary

    pycoQC \\
        $args \\
        -f "sequencing_summary.txt" \\
        $barcode_me \\
        -o ${prefix}.html \\
        -j ${prefix}.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pycoqc: \$(pycoQC --version 2>&1 | sed 's/^.*pycoQC v//; s/ .*\$//')
    END_VERSIONS
    """
}

process PORECHOP_PORECHOP {
    tag "$meta"
    label 'process_medium'
    cache 'lenient'
    
    conda "${projectDir}/conda_environments/porechop.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/porechop:0.2.4--py39h7cff6ad_2' :
        'biocontainers/porechop:0.2.4--py39h7cff6ad_2' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: reads
    tuple val(meta), path("*.log")     , emit: log
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta}"
    """
    porechop \\
        -i $reads \\
        -t $task.cpus \\
        $args \\
        -o ${prefix}.PORECHOP.fastq.gz \\
        > ${prefix}.log
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        porechop: \$( porechop --version )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta}"
    """
    touch ${prefix}.fastq
    gzip ${prefix}.fastq
    touch ${prefix}.log
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        porechop: \$( porechop --version )
    END_VERSIONS
    """
}

process UNICYCLER {
    tag "$meta"
    label 'process_high'
    cache 'lenient'

    conda "${projectDir}/conda_environments/unicycler.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/unicycler:0.4.8--py38h8162308_3' :
        'biocontainers/unicycler:0.4.8--py38h8162308_3' }"

    input:
    tuple val(meta), path(shortreads), path(longreads)

    output:
    tuple val(meta), path('*.scaffolds.fa.gz'), emit: scaffolds
    tuple val(meta), path('*.assembly.gfa.gz'), emit: gfa
    tuple val(meta), path('*.log')            , emit: log
    path  "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta}"
    if(params.assembly_type == 'long'){
        input_reads = "-l $longreads"
    } else if (params.assembly_type == 'short'){
        input_reads = "-1 ${shortreads[0]} -2 ${shortreads[1]}"
    } else if (params.assembly_type == 'hybrid'){
        input_reads = "-1 ${shortreads[0]} -2 ${shortreads[1]} -l $longreads"
    }
    """
    unicycler \\
        --threads $task.cpus \\
        $args \\
        $input_reads \\
        --out ./

    mv assembly.fasta ${prefix}.scaffolds.fa
    gzip -n ${prefix}.scaffolds.fa
    mv assembly.gfa ${prefix}.assembly.gfa
    gzip -n ${prefix}.assembly.gfa
    mv unicycler.log ${prefix}.unicycler.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        unicycler: \$(echo \$(unicycler --version 2>&1) | sed 's/^.*Unicycler v//; s/ .*\$//')
    END_VERSIONS
    """
}

process CANU {
    tag "$meta"
    label 'process_high'
    cache 'lenient'

    conda "${projectDir}/conda_environments/canu.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/canu:2.2--ha47f30e_0':
        'biocontainers/canu:2.2--ha47f30e_0' }"

    input:
    tuple val(meta), path(reads)
    val mode
    val genomesize

    output:
    tuple val(meta), path("*.report")                   , emit: report
    tuple val(meta), path("*.contigs.fasta.gz")         , emit: assembly                , optional: true
    tuple val(meta), path("*.unassembled.fasta.gz")     , emit: contigs
    tuple val(meta), path("*.correctedReads.fasta.gz")	, emit: corrected_reads         , optional: true
    tuple val(meta), path("*.trimmedReads.fasta.gz")	, emit: corrected_trimmed_reads , optional: true
    tuple val(meta), path("*.contigs.layout")           , emit: metadata                , optional: true
    tuple val(meta), path("*.contigs.layout.readToTig") , emit: contig_position         , optional: true
    tuple val(meta), path("*.contigs.layout.tigInfo")   , emit: contig_info             , optional: true
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta}"
    def valid_mode = ["-pacbio", "-nanopore", "-pacbio-hifi"]
    if ( !valid_mode.contains(mode) )  { error "Unrecognised mode to run Canu. Options: ${valid_mode.join(', ')}" }
    """
    canu \\
        -p ${prefix} \\
        $mode \\
        genomeSize=${genomesize} \\
        $args \\
        maxThreads=$task.cpus \\
        $reads

    gzip *.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        canu: \$(echo \$(canu --version 2>&1) | sed 's/^.*canu //; s/Using.*\$//' )
    END_VERSIONS
    """
}

process MINIMAP2_ALIGN {
    tag "$meta"
    label 'process_medium'
    cache 'lenient'

    // Note: the versions here need to match the versions used in the mulled container below and minimap2/index
    conda "${projectDir}/conda_environments/minimap2_align.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:365b17b986c1a60c1b82c6066a9345f38317b763-0' :
        'biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:365b17b986c1a60c1b82c6066a9345f38317b763-0' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(reference)
    val bam_format
    val cigar_paf_format
    val cigar_bam

    output:
    tuple val(meta), path("*.paf"), optional: true, emit: paf
    tuple val(meta), path("*.bam"), optional: true, emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta}"
    def bam_output = bam_format ? "-a | samtools sort -@ ${task.cpus} -o ${prefix}.bam ${args2}" : "-o ${prefix}.paf"
    def cigar_paf = cigar_paf_format && !bam_format ? "-c" : ''
    def set_cigar_bam = cigar_bam && bam_format ? "-L" : ''
    """
    minimap2 \\
        $args \\
        -t $task.cpus \\
        ${reference ?: reads} \\
        $reads \\
        $cigar_paf \\
        $set_cigar_bam \\
        $bam_output


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta}"
    def output_file = bam_format ? "${prefix}.bam" : "${prefix}.paf"
    """
    touch $output_file

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
    END_VERSIONS
    """
}


process MINIASM {
    tag "$meta"
    label 'process_high'
    cache 'lenient'

    conda "${projectDir}/conda_environments/miniasm.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/miniasm:0.3_r179--h5bf99c6_2' :
        'biocontainers/miniasm:0.3_r179--h5bf99c6_2' }"

    input:
    tuple val(meta), path(reads), path(paf)

    output:
    tuple val(meta), path("*.gfa.gz")  , emit: gfa
    tuple val(meta), path("*.fasta.gz"), emit: assembly
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta}"
    """
    miniasm \\
        $args \\
        -f $reads \\
        $paf > \\
        ${prefix}.gfa

    awk '/^S/{print ">"\$2"\\n"\$3}' "${prefix}.gfa" | fold > ${prefix}.fasta

    gzip -n ${prefix}.gfa
    gzip -n ${prefix}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        miniasm: \$( miniasm -V 2>&1 )
    END_VERSIONS
    """
}

process DRAGONFLYE {
    tag "$meta"
    label 'process_medium'
    cache 'lenient'

    conda "${projectDir}/conda_environments/dragonflye.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/dragonflye:1.1.2--hdfd78af_0' :
        'biocontainers/dragonflye:1.1.2--hdfd78af_0' }"

    input:
    tuple val(meta), path(shortreads), path(longreads)

    output:
    tuple val(meta), path("*.fa")                                               , emit: contigs
    tuple val(meta), path("dragonflye.log")                                     , emit: log
    tuple val(meta), path("{flye,miniasm,raven}.fasta")                         , emit: raw_contigs
    tuple val(meta), path("{flye,miniasm,raven}-unpolished.gfa"), optional:true , emit: gfa
    tuple val(meta), path("flye-info.txt")                      , optional:true , emit: txt
    path "versions.yml"                                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta}"
    def memory  = task.memory.toGiga()
    def shortreads_polishing = shortreads ? "--R1 ${shortreads[0]} --R2 ${shortreads[1]}" : ''
    """
    dragonflye \\
        --reads ${longreads} \\
        $shortreads_polishing \\
        $args \\
        --prefix ${prefix} \\
        --cpus $task.cpus \\
        --ram $memory \\
        --outdir ./ \\
        --force

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragonflye: \$(dragonflye --version 2>&1 | sed 's/^.*dragonflye //' )
    END_VERSIONS
    """
}

process RACON {
    tag "$meta"
    label 'process_high'
    cache 'lenient'

    conda "${projectDir}/conda_environments/racon.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/racon:1.4.20--h9a82719_1' :
        'biocontainers/racon:1.4.20--h9a82719_1' }"

    input:
    tuple val(meta), path(reads), path(assembly), path(paf)

    output:
    tuple val(meta), path('*_assembly_consensus.fasta.gz') , emit: improved_assembly
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta}"
    """
    racon -t "$task.cpus" \\
        "${reads}" \\
        "${paf}" \\
        $args \\
        "${assembly}" > \\
        ${prefix}_assembly_consensus.fasta

    gzip -n ${prefix}_assembly_consensus.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        racon: \$( racon --version 2>&1 | sed 's/^.*v//' )
    END_VERSIONS
    """
}


process SAMTOOLS_SORT {
    tag "$meta"
    label 'process_medium'
    cache 'lenient'

    conda "${projectDir}/conda_environments/samtools_sort.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.19.2--h50ea8bc_0' :
        'biocontainers/samtools:1.19.2--h50ea8bc_0' }"

    input:
    tuple val(meta) , path(bam)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("*.bam"),     emit: bam,  optional: true
    tuple val(meta), path("*.cram"),    emit: cram, optional: true
    tuple val(meta), path("*.crai"),    emit: crai, optional: true
    tuple val(meta), path("*.csi"),     emit: csi,  optional: true
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta}"
    def extension = args.contains("--output-fmt sam") ? "sam" :
                    args.contains("--output-fmt cram") ? "cram" :
                    "bam"
    def reference = fasta ? "--reference ${fasta}" : ""
    if ("$bam" == "${prefix}.bam") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"

    """
    samtools cat \\
        --threads $task.cpus \\
        ${bam} \\
    | \\
    samtools sort \\
        $args \\
        -T ${prefix} \\
        --threads $task.cpus \\
        ${reference} \\
        -o ${prefix}.${extension} \\
        -

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta}"
    """
    touch ${prefix}.bam
    touch ${prefix}.bam.csi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}


process SAMTOOLS_INDEX {
    tag "$meta"
    label 'process_low'
    cache 'lenient'

    conda "${projectDir}/conda_environments/samtools_index.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.19.2--h50ea8bc_0' :
        'biocontainers/samtools:1.19.2--h50ea8bc_0' }"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*.bai") , optional:true, emit: bai
    tuple val(meta), path("*.csi") , optional:true, emit: csi
    tuple val(meta), path("*.crai"), optional:true, emit: crai
    path  "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    samtools \\
        index \\
        -@ ${task.cpus-1} \\
        $args \\
        $input

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch ${input}.bai
    touch ${input}.crai
    touch ${input}.csi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}

process KRAKEN2_DB_PREPARATION {
    tag "${db.simpleName}"
    label 'process_low'
    cache 'lenient'

    conda "${projectDir}/conda_environments/kraken2_db_preparation.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    path db

    output:
    tuple val("${db.simpleName}"), path("database"), emit: db
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir db_tmp
    tar -xf "${db}" -C db_tmp
    mkdir database
    mv `find db_tmp/ -name "*.k2d"` database/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tar: \$(tar --version | sed -n 's/^tar (GNU tar) \$([0-9.]*\$).*/\$1/p')
    END_VERSIONS
    """
}

process KRAKEN2 {
    tag "$meta"
    label 'process_high'
    cache 'lenient'

    conda "${projectDir}/conda_environments/kraken2.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-5799ab18b5fc681e75923b2450abaa969907ec98:87fc08d11968d081f3e8a37131c1f1f6715b6542-0' :
        'biocontainers/mulled-v2-5799ab18b5fc681e75923b2450abaa969907ec98:87fc08d11968d081f3e8a37131c1f1f6715b6542-0' }"

    input:
    tuple val(meta), path(reads)
    path  db
    val save_output_fastqs
    val save_reads_assignment

    output:
    tuple val(meta), path('*.classified{.,_}*')     , optional:true, emit: classified_reads_fastq
    tuple val(meta), path('*.unclassified{.,_}*')   , optional:true, emit: unclassified_reads_fastq
    tuple val(meta), path('*classifiedreads.txt')   , optional:true, emit: classified_reads_assignment
    tuple val(meta), path('*report.txt')                           , emit: report
    path "versions.yml"                                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta}"
    def paired       = params.single_end ? "" : "--paired"
    def classified   = params.single_end ? "${prefix}.classified.fastq"   : "${prefix}.classified#.fastq"
    def unclassified = params.single_end ? "${prefix}.unclassified.fastq" : "${prefix}.unclassified#.fastq"
    def classified_option = save_output_fastqs ? "--classified-out ${classified}" : ""
    def unclassified_option = save_output_fastqs ? "--unclassified-out ${unclassified}" : ""
    def readclassification_option = save_reads_assignment ? "--output ${prefix}.kraken2.classifiedreads.txt" : "--output /dev/null"
    def compress_reads_command = save_output_fastqs ? "pigz -p $task.cpus *.fastq" : ""

    """
    kraken2 \\
        --db $db \\
        --threads $task.cpus \\
        --report ${prefix}.kraken2.report.txt \\
        --gzip-compressed \\
        $unclassified_option \\
        $classified_option \\
        $readclassification_option \\
        $paired \\
        $args \\
        $reads

    $compress_reads_command

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken2: \$(echo \$(kraken2 --version 2>&1) | sed 's/^.*Kraken version //; s/ .*\$//')
        pigz: \$( pigz --version 2>&1 | sed 's/pigz //g' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta}"
    def paired       = params.single_end ? "" : "--paired"
    def classified   = params.single_end ? "${prefix}.classified.fastq.gz"   : "${prefix}.classified_1.fastq.gz ${prefix}.classified_2.fastq.gz"
    def unclassified = params.single_end ? "${prefix}.unclassified.fastq.gz" : "${prefix}.unclassified_1.fastq.gz ${prefix}.unclassified_2.fastq.gz"
    def readclassification_option = save_reads_assignment ? "--output ${prefix}.kraken2.classifiedreads.txt" : "--output /dev/null"
    def compress_reads_command = save_output_fastqs ? "pigz -p $task.cpus *.fastq" : ""

    """
    touch ${prefix}.kraken2.report.txt
    if [ "$save_output_fastqs" == "true" ]; then
        touch $classified
        touch $unclassified
    fi
    if [ "$save_reads_assignment" == "true" ]; then
        touch ${prefix}.kraken2.classifiedreads.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken2: \$(echo \$(kraken2 --version 2>&1) | sed 's/^.*Kraken version //; s/ .*\$//')
        pigz: \$( pigz --version 2>&1 | sed 's/pigz //g' )
    END_VERSIONS
    """

}


process KRAKEN2_LONG {
    tag "$meta"
    label 'process_high'
    cache 'lenient'

    conda "${projectDir}/conda_environments/kraken2.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-5799ab18b5fc681e75923b2450abaa969907ec98:87fc08d11968d081f3e8a37131c1f1f6715b6542-0' :
        'biocontainers/mulled-v2-5799ab18b5fc681e75923b2450abaa969907ec98:87fc08d11968d081f3e8a37131c1f1f6715b6542-0' }"

    input:
    tuple val(meta), path(reads)
    path  db
    val save_output_fastqs
    val save_reads_assignment

    output:
    tuple val(meta), path('*.classified{.,_}*')     , optional:true, emit: classified_reads_fastq
    tuple val(meta), path('*.unclassified{.,_}*')   , optional:true, emit: unclassified_reads_fastq
    tuple val(meta), path('*classifiedreads.txt')   , optional:true, emit: classified_reads_assignment
    tuple val(meta), path('*report.txt')                           , emit: report
    path "versions.yml"                                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta}"
    def paired = ""
    def classified_option = save_output_fastqs ? "--classified-out ${classified}" : ""
    def unclassified_option = save_output_fastqs ? "--unclassified-out ${unclassified}" : ""
    def readclassification_option = save_reads_assignment ? "--output ${prefix}.kraken2.classifiedreads.txt" : "--output /dev/null"
    def compress_reads_command = save_output_fastqs ? "pigz -p $task.cpus *.fastq" : ""

    """
    kraken2 \\
        --db $db \\
        --threads $task.cpus \\
        --report ${prefix}.kraken2.report.txt \\
        --gzip-compressed \\
        ${prefix}.unclassified.fastq \\
        ${prefix}.classified.fastq \\
        $readclassification_option \\
     //   $paired \\
        $args \\
        $reads

    $compress_reads_command

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken2: \$(echo \$(kraken2 --version 2>&1) | sed 's/^.*Kraken version //; s/ .*\$//')
        pigz: \$( pigz --version 2>&1 | sed 's/pigz //g' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta}"
    def readclassification_option = save_reads_assignment ? "--output ${prefix}.kraken2.classifiedreads.txt" : "--output /dev/null"
    def compress_reads_command = save_output_fastqs ? "pigz -p $task.cpus *.fastq" : ""

    """
    touch ${prefix}.kraken2.report.txt
    if [ "$save_output_fastqs" == "true" ]; then
        touch $classified
        touch $unclassified
    fi
    if [ "$save_reads_assignment" == "true" ]; then
        touch ${prefix}.kraken2.classifiedreads.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken2: \$(echo \$(kraken2 --version 2>&1) | sed 's/^.*Kraken version //; s/ .*\$//')
        pigz: \$( pigz --version 2>&1 | sed 's/pigz //g' )
    END_VERSIONS
    """

}


process QUAST {
    tag "$meta"
    label 'process_medium'
    cache 'lenient'

    conda "${projectDir}/conda_environments/quast.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/quast:5.2.0--py39pl5321h2add14b_1' :
        'biocontainers/quast:5.2.0--py39pl5321h2add14b_1' }"

    input:
    tuple val(meta) , path(consensus)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(gff)

    output:
    tuple val(meta), path("${prefix}")                   , emit: results
    tuple val(meta), path("${prefix}.tsv")               , emit: tsv
    tuple val(meta), path("${prefix}_transcriptome.tsv") , optional: true , emit: transcriptome
    tuple val(meta), path("${prefix}_misassemblies.tsv") , optional: true , emit: misassemblies
    tuple val(meta), path("${prefix}_unaligned.tsv")     , optional: true , emit: unaligned
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args   ?: ''
    prefix        = task.ext.prefix ?: "${meta}"
    def features  = gff             ?  "--features $gff" : ''
    def reference = fasta           ?  "-r $fasta"       : ''
    """
    quast.py \\
        --output-dir $prefix \\
        $reference \\
        $features \\
        --threads $task.cpus \\
        $args \\
        ${consensus.join(' ')}

    ln -s ${prefix}/report.tsv ${prefix}.tsv
    [ -f  ${prefix}/contigs_reports/all_alignments_transcriptome.tsv ] && ln -s ${prefix}/contigs_reports/all_alignments_transcriptome.tsv ${prefix}_transcriptome.tsv
    [ -f  ${prefix}/contigs_reports/misassemblies_report.tsv         ] && ln -s ${prefix}/contigs_reports/misassemblies_report.tsv ${prefix}_misassemblies.tsv
    [ -f  ${prefix}/contigs_reports/unaligned_report.tsv             ] && ln -s ${prefix}/contigs_reports/unaligned_report.tsv ${prefix}_unaligned.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quast: \$(quast.py --version 2>&1 | sed 's/^.*QUAST v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args      = task.ext.args   ?: ''
    prefix        = task.ext.prefix ?: "${meta}"
    def features  = gff             ? "--features $gff" : ''
    def reference = fasta           ? "-r $fasta" : ''

    """
    mkdir -p $prefix
    touch $prefix/report.tsv
    touch $prefix/report.html
    touch $prefix/report.pdf
    touch $prefix/quast.log
    touch $prefix/transposed_report.txt
    touch $prefix/transposed_report.tex
    touch $prefix/icarus.html
    touch $prefix/report.tex
    touch $prefix/report.txt

    mkdir -p $prefix/basic_stats
    touch $prefix/basic_stats/cumulative_plot.pdf
    touch $prefix/basic_stats/Nx_plot.pdf
    touch $prefix/basic_stats/genome_GC_content_plot.pdf
    touch $prefix/basic_stats/GC_content_plot.pdf

    mkdir -p $prefix/icarus_viewers
    touch $prefix/icarus_viewers/contig_size_viewer.html

    ln -s $prefix/report.tsv ${prefix}.tsv

    if [ $fasta ]; then
        touch $prefix/basic_stats/NGx_plot.pdf
        touch $prefix/basic_stats/gc.icarus.txt

        mkdir -p $prefix/aligned_stats
        touch $prefix/aligned_stats/NAx_plot.pdf
        touch $prefix/aligned_stats/NGAx_plot.pdf
        touch $prefix/aligned_stats/cumulative_plot.pdf

        mkdir -p $prefix/contigs_reports
        touch $prefix/contigs_reports/all_alignments_transcriptome.tsv
        touch $prefix/contigs_reports/contigs_report_transcriptome.mis_contigs.info
        touch $prefix/contigs_reports/contigs_report_transcriptome.stderr
        touch $prefix/contigs_reports/contigs_report_transcriptome.stdout
        touch $prefix/contigs_reports/contigs_report_transcriptome.unaligned.info
        mkdir -p $prefix/contigs_reports/minimap_output
        touch $prefix/contigs_reports/minimap_output/transcriptome.coords
        touch $prefix/contigs_reports/minimap_output/transcriptome.coords.filtered
        touch $prefix/contigs_reports/minimap_output/transcriptome.coords_tmp
        touch $prefix/contigs_reports/minimap_output/transcriptome.sf
        touch $prefix/contigs_reports/minimap_output/transcriptome.unaligned
        touch $prefix/contigs_reports/minimap_output/transcriptome.used_snps
        touch $prefix/contigs_reports/misassemblies_frcurve_plot.pdf
        touch $prefix/contigs_reports/misassemblies_plot.pdf
        touch $prefix/contigs_reports/misassemblies_report.tex
        touch $prefix/contigs_reports/misassemblies_report.tsv
        touch $prefix/contigs_reports/misassemblies_report.txt
        touch $prefix/contigs_reports/transcriptome.mis_contigs.fa
        touch $prefix/contigs_reports/transposed_report_misassemblies.tex
        touch $prefix/contigs_reports/transposed_report_misassemblies.tsv
        touch $prefix/contigs_reports/transposed_report_misassemblies.txt
        touch $prefix/contigs_reports/unaligned_report.tex
        touch $prefix/contigs_reports/unaligned_report.tsv
        touch $prefix/contigs_reports/unaligned_report.txt

        mkdir -p $prefix/genome_stats
        touch $prefix/genome_stats/genome_info.txt
        touch $prefix/genome_stats/transcriptome_gaps.txt
        touch $prefix/icarus_viewers/alignment_viewer.html

        ln -sf ${prefix}/contigs_reports/misassemblies_report.tsv ${prefix}_misassemblies.tsv
        ln -sf ${prefix}/contigs_reports/unaligned_report.tsv ${prefix}_unaligned.tsv
        ln -sf ${prefix}/contigs_reports/all_alignments_transcriptome.tsv ${prefix}_transcriptome.tsv

    fi

    if ([ $fasta ] && [ $gff ]); then
        touch $prefix/genome_stats/features_cumulative_plot.pdf
        touch $prefix/genome_stats/features_frcurve_plot.pdf
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quast: \$(quast.py --version 2>&1 | sed 's/^.*QUAST v//; s/ .*\$//')
    END_VERSIONS
    """
}


process MULTIQC {
    label 'process_single'
    cache 'lenient'

    conda "${projectDir}/conda_environments/multiqc.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.21--pyhdfd78af_0' :
        'biocontainers/multiqc:1.21--pyhdfd78af_0' }"

    input:
    path  multiqc_files, stageAs: "?/*"
    path(multiqc_config)
    path(extra_multiqc_config)
    path(multiqc_logo)

    output:
    path "*multiqc_report.html", emit: report
    path "*_data"              , emit: data
    path "*_plots"             , optional:true, emit: plots
    path "versions.yml"        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def config = multiqc_config ? "--config $multiqc_config" : ''
    def extra_config = extra_multiqc_config ? "--config $extra_multiqc_config" : ''
    def logo = multiqc_logo ? /--cl-config 'custom_logo: "${multiqc_logo}"'/ : ''
    """
    multiqc \\
        --force \\
        $args \\
        $config \\
        $extra_config \\
        $logo \\
        .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """

    stub:
    """
    mkdir multiqc_data
    touch multiqc_plots
    touch multiqc_report.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """
}


process NANOPOLISH {
    tag "$meta"
    label 'process_high'
    cache 'lenient'

    conda "${projectDir}/conda_environments/nanopolish.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nanopolish:0.14.0--h773013f_3' :
        'biocontainers/nanopolish:0.14.0--h773013f_3' }"

    input:
    tuple val(meta), val(reads), file(longreads), file(assembly), file(bam), file(bai), file(fast5)

    output:
    tuple val(meta), file('polished_genome.fa') , emit: assembly
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta}"
    """
    nanopolish index -d "${fast5}" "${longreads}"

    nanopolish variants \
        --consensus \
        -o polished.vcf \
        -r "${longreads}" \
        -b "${bam}" \
        -g "${assembly}" \
        -t "${task.cpus}" \
        --min-candidate-frequency 0.1 \
        $args

    nanopolish vcf2fasta -g "${assembly}" polished.vcf > polished_genome.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanopolish: \$( nanopolish --version | sed -e "s/nanopolish version //g" | head -n 1 )
    END_VERSIONS
    """
}


process MEDAKA {
    tag "$meta"
    label 'process_high'
    cache 'lenient'

    conda "${projectDir}/conda_environments/medaka.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/medaka:1.4.3--py38h130def0_0' :
        'biocontainers/medaka:1.4.3--py38h130def0_0' }"

    input:
    tuple val(meta), file(longreads), file(assembly)

    output:
    tuple val(meta), path('*_polished_genome.fa')   , emit: assembly
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args                    = task.ext.args ?: ''
    def prefix                  = task.ext.prefix ?: "${meta}"
    def reads_bgzip_command     = ("$longreads".endsWith('.gz')) ? "zcat $longreads | bgzip -c > ${prefix}.fastq.bgz" : ''
    def assembly_bgzip_command  = ("$assembly".endsWith('.gz'))  ? "zcat $assembly  | bgzip -c > ${prefix}.fasta.bgz" : ''
    if ("$longreads".endsWith('.gz')) { reads_bgzip_out     = "${prefix}.fastq.bgz"} else { reads_bgzip_out    = null }
    if ("$assembly".endsWith('.gz'))  { assembly_bgzip_out  = "${prefix}.fasta.bgz"} else { assembly_bgzip_out = null }

    """
    # Recompress with bgzip
    $reads_bgzip_command
    $assembly_bgzip_command

    medaka_consensus $args \
        -i ${ reads_bgzip_out ?: longreads } \
        -d ${ assembly_bgzip_out ?: assembly } \
        -o "${prefix}_polished_genome.fa" \
        -t $task.cpus

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        medaka: \$( medaka --version 2>&1 | sed 's/medaka //g' )
    END_VERSIONS
    """
}


process CUSTOM_DUMPSOFTWAREVERSIONS {
    
    publishDir "${params.outdir}", mode:'copy'
    
    input:
    path versions

    output:
    path "software_versions.yml"    , emit: yml_ch
    path "software_versions_mqc.yml", emit: mqc_yml_ch
    path "versions.yml"             , emit: versions_ch

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    template 'dumpsoftwareversions.py'
}
