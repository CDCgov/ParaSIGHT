process CAT_BLAST_REPORTS {
    label 'process_low'

    conda "conda-forge::perl=5.32.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/perl:5.26.2' :
        'quay.io/biocontainers/perl:5.26.2' }"

    input:
    tuple val(meta), path(txt)

    output:
	path "concat_blast.txt"               , emit: report
	
    when:
    task.ext.when == null || task.ext.when

    script:  
    def args = task.ext.args ?: ''
    """
	cat $txt | sed "s/^/${meta.id}\\t/g" >> concat_blast.txt
    """
}