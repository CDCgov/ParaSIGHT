process SUMMARIZE_ALLELES {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::pandas=2.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.2.1':
        'biocontainers/pandas:2.2.1' }"

    input:
    tuple val(meta),  path(gene_distribution)
    tuple val(meta2), path(alleles_per_genome)
    tuple val(meta3), path(annotations)
    path(samplesheet)
    path(samplesheet2)

    output:
    tuple val(meta), path("*.alleles_dashboard.csv")             , emit: alleles_dashboard
    tuple val(meta), path("*.alleles_per_genome_dashboard.csv")  , emit: per_genome_dashboard

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def csv2 = samplesheet2 ? "--meta2 $samplesheet2" : ""
    def annots = annotations ? "--annots $annotations" : ""
    """

    python ${projectDir}/bin/allele_summary_and_sort.py \\
        --meta $samplesheet \\
        $csv2 \\
        $annots \\
        --genes $gene_distribution \\
        --alleles $alleles_per_genome \\
        --output ${prefix}
    """
}