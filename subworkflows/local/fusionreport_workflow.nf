//
// Check input samplesheet and get read channels
//

include { FUSIONREPORT }                           from '../../modules/local/fusionreport/detect/main'


workflow FUSIONREPORT_WORKFLOW {
    take:
        reads
        fusionreport_ref
        fusions

    main:
        ch_versions = Channel.empty()


    tools = {fusions == [] ? "--arriba ${fusions} " : ''}
    // tools += !ericscript.empty() ? "--ericscript ${ericscript} " : ''
    // tools += !fusioncatcher.empty() ? "--fusioncatcher ${fusioncatcher} " : ''
    // tools += !pizzly.empty() ? "--pizzly ${pizzly} " : ''
    // tools += !squid.empty() ? "--squid ${squid} " : ''
    // tools += !starfusion.empty() ? "--starfusion ${starfusion} " : ''

        // def tools = !arriba ? "--arriba ${arriba_fusions} " : ''

        FUSIONREPORT( reads, fusionreport_ref, tools )
        // ch_versions = ch_versions.mix(STAR_FOR_STARFUSION.out.versions)
        // reads_junction = reads.join(STAR_FOR_STARFUSION.out.junction )

        // STARFUSION( reads_junction, starfusion_ref)
        // ch_versions = ch_versions.mix(STARFUSION.out.versions, starfusion_ref )

    emit:
        versions        = ch_versions.ifEmpty(null) // channel: [ versions.yml ]

    // versions = ch_versions.ifEmpty(null)
    // STARFUSION.out.fusions
    // STARFUSION.out.abridged
    // STARFUSION.out.coding_effect

}

