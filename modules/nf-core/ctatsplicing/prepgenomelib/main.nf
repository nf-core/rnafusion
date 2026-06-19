process CTATSPLICING_PREPGENOMELIB {
    tag "$meta.id"
    label 'process_single'
    stageInMode 'copy'

    container "quay.io/nf-core/ctatsplicing:0.0.3"

    input:
    tuple val(meta), path(genome_lib)
    path(cancer_intron_tsv)

    output:
    tuple val(meta), path(genome_lib, includeInputs:true), emit: reference
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple val("${task.process}"), val('ctatsplicing'), val("0.0.3"), emit: versions_ctatsplicing, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "CTATSPLICING_PREPGENOMELIB module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    """
    /usr/local/src/CTAT-SPLICING/prep_genome_lib/ctat-splicing-lib-integration.py \\
        --cancer_introns_tsv $cancer_intron_tsv \\
        --genome_lib_dir $genome_lib
    """

    stub:
    """
    touch $genome_lib/refGene.bed
    touch $genome_lib/refGene.sort.bed.gz.tbi
    """
}
