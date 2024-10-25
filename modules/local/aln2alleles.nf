process ALN2ALLELES {
    tag "${meta.id}"
    label 'process_low'

    conda "bioconda::biopython=1.81"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.81':
        'biocontainers/biopython:1.81' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.alleles.fasta")  , emit: alleles_fasta
    path("*alleles_per_genome.tab")           , emit: alleles_report

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python ${projectDir}/bin/alignment2alleles.py \\
        $fasta \\
        ${prefix} \\
        ${prefix}.alleles.fasta
    """
}