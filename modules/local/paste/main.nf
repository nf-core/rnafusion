process PASTE {
    tag "$meta.id"
    label 'process_medium'

    input:
    path input

    output:
    path output

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    paste ${input} ${input} | cut -f12,16-25 > ${refflat}
    """
}
