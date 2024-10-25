process GROUP_GENE_CONTENT_DIFFERENCES {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::pandas=2.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.2.1':
        'biocontainers/pandas:2.2.1' }"

    input:
    tuple val(meta), path(phyletic_matrix)
    path(samplesheet)
    path(samplesheet2)

    output:
    tuple val(meta), path("*.unique_group_differences.tab")  , emit: unique_group_differences
    tuple val(meta), path("*.pairwise_group_differences.tab"), emit: pairwise_group_differences
    tuple val(meta), path("*.group_unique_genes.tab")        , emit: unique_group_genes
    tuple val(meta), path("*.group_union_genes.tab")         , emit: union_group_genes

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def meta2 = samplesheet2 ? "--meta2 $samplesheet2" : ""
    """

    python ${projectDir}/bin/group_gene_content_differences.py \\
        --meta $samplesheet \\
        $meta2 \\
        --matrix $phyletic_matrix \\
        --output ${prefix}
    """
}