process GET_PATH {
    tag 'getpath'
    label 'process_low'

    input:
    tuple val(meta), path(file)

    output:
    path file                , emit: file

    script:
    """
    """
}
