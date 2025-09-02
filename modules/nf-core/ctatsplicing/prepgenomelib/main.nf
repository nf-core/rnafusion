process CTATSPLICING_PREPGENOMELIB {
    tag "$meta.id"
    label 'process_single'

    container "nf-core/ctatsplicing:0.0.3"

    input:
    tuple val(meta), path(genome_lib)
    path(cancer_intron_tsv)

    output:
    tuple val(meta), path(genome_lib, includeInputs:true), emit: reference
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "CTATSPLICING_PREPGENOMELIB module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '0.0.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    /usr/local/src/CTAT-SPLICING/prep_genome_lib/ctat-splicing-lib-integration.py \\
        --cancer_introns_tsv $cancer_intron_tsv \\
        --genome_lib_dir $genome_lib

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ctatsplicing: $VERSION
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '0.0.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch $genome_lib/refGene.bed
    echo | gzip > $genome_lib/refGene.sort.bed.gz
    touch $genome_lib/refGene.sort.bed.gz.tbi
    mkdir $genome_lib/cancer_splicing_lib
    touch $genome_lib/cancer_splicing_lib/cancer_splicing.idx

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ctatsplicing: $VERSION
    END_VERSIONS
    """
}
