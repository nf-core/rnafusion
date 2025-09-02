//
// Extract descriptive values from BAMs
//

include { PICARD_COLLECTRNASEQMETRICS       } from '../../../modules/nf-core/picard/collectrnaseqmetrics'
include { GATK4_MARKDUPLICATES              } from '../../../modules/nf-core/gatk4/markduplicates'
include { PICARD_COLLECTINSERTSIZEMETRICS   } from '../../../modules/nf-core/picard/collectinsertsizemetrics'

workflow QC_WORKFLOW {
    take:
        ch_bam_sorted           // channel [ meta, bam        ]
        ch_refflat              // channel [ meta, refflat    ]
        ch_fasta                // channel [ meta, fasta      ]
        ch_fai                  // channel [ meta, fai        ]
        ch_rrna_interval        // channel [ meta, interval   ]

    main:
        ch_versions = Channel.empty()

        PICARD_COLLECTRNASEQMETRICS(
            ch_bam_sorted,
            ch_refflat.map{ _meta, refflat -> [ refflat ] },
            ch_fasta.map{ _meta, fasta -> [ fasta ] },
            ch_rrna_interval.map{ _meta, intervals -> [ intervals ] }.ifEmpty([])
        ) // Some chromosome or annotation may not have rRNA genes
        ch_versions = ch_versions.mix(PICARD_COLLECTRNASEQMETRICS.out.versions)
        ch_rnaseq_metrics = PICARD_COLLECTRNASEQMETRICS.out.metrics

        GATK4_MARKDUPLICATES(
            ch_bam_sorted,
            ch_fasta.map { _meta, fasta -> [ fasta ]},
            ch_fai.map { _meta, fasta_fai -> [ fasta_fai ]}
        )
        ch_versions = ch_versions.mix(GATK4_MARKDUPLICATES.out.versions)
        ch_duplicate_metrics = GATK4_MARKDUPLICATES.out.metrics

        PICARD_COLLECTINSERTSIZEMETRICS(
            ch_bam_sorted
        )
        ch_versions = ch_versions.mix(PICARD_COLLECTINSERTSIZEMETRICS.out.versions)
        ch_insertsize_metrics = PICARD_COLLECTINSERTSIZEMETRICS.out.metrics

    emit:
        versions            = ch_versions            // channel [ path       ]
        rnaseq_metrics      = ch_rnaseq_metrics      // channel [ meta, path ]
        duplicate_metrics   = ch_duplicate_metrics   // channel [ meta, path ]
        insertsize_metrics  = ch_insertsize_metrics  // channel [ meta, path ]

}
