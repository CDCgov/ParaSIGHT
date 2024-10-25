process PARSE_BLAST_REPORT {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::pandas=2.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.2.1':
        'biocontainers/pandas:2.2.1' }"

    input:
    tuple val(meta), path(report)

    output:
    tuple val(meta), path("*.distribution.tab"), emit: gene_content_report

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // If this is giving errors, check that the bin/ dir is in PATH
    """
    python ${projectDir}/bin/blastn2coregenome.py \\
        --input $report \\
        --output ${prefix}.distribution.tab \\
        $args
    """
}
