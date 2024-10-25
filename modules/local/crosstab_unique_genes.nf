process CROSSTAB_UNIQUE_GENES {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::pandas=2.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.2.1':
        'biocontainers/pandas:2.2.1' }"

    input:
    tuple val(meta), path(genes)
    path(eggnog)
    path(keggs)

    output:
    tuple val(meta), path("*gene_annotations.tab"), emit: group_differences

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python ${projectDir}/bin/crosstab_unique_genes.py \\
        --genes $genes \\
        --annots $eggnog \\
        --keggs $keggs \\
        --output ${prefix}
    """
}