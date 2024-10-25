process BLAST2FASTA {
    tag "gene_content"
    label 'process_low'

    conda "conda-forge::pandas=2.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.2.1':
        'biocontainers/pandas:2.2.1' }"

    input:
    tuple val(meta), path(report)

    output:
    path("*.fasta")                   , emit: aln

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // If this is giving errors, check that the bin/ dir is in PATH
    """
    python ${projectDir}/bin/blastn2aln.py $report ./
    """
}
