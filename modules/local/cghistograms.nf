process CGHISTOGRAMS {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::seaborn=0.13.2"
    container "${projectDir}/assets/seaborn-0.11.0.sif"

    input:
    tuple val(meta), path(distribution)
    tuple val(meta), path(phyletic_matrix)
    path(samplesheet)
    path(samplesheet2)

    output:
    tuple val(meta), path("*.histogram.csv"), emit: cgHistograms

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def meta2 = samplesheet2? "--meta2 $samplesheet2" : ""
    """
    python ${projectDir}/bin/cgHistograms.py \\
        --genes $distribution \\
        --matrix $phyletic_matrix \\
        --meta $samplesheet \\
        $meta2 \\
        $args \\
        --output ${prefix}.histogram.csv
    """
}