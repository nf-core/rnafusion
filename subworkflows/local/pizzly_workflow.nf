include { KALLISTO_QUANT    }     from '../../modules/local/kallisto/quant/main'
include { PIZZLY            }     from '../../modules/local/pizzly/detect/main'
include { GET_PATH          }     from '../../modules/local/getpath/main'


workflow PIZZLY_WORKFLOW {
    take:
        reads

    main:
        ch_versions = Channel.empty()
        ch_dummy_file = file("$baseDir/assets/dummy_file_pizzly.txt", checkIfExists: true)

        if (params.pizzly || params.all) {
            if (params.pizzly_fusions) {
                ch_pizzly_fusions = params.pizzly_fusions
            } else {
                index ="${params.pizzly_ref}/kallisto"
                gtf ="${params.ensembl_ref}/Homo_sapiens.GRCh38.${params.ensembl_version}.gtf"
                transcript ="${params.ensembl_ref}/Homo_sapiens.GRCh38.${params.ensembl_version}.cdna.all.fa.gz"

                KALLISTO_QUANT( reads, index )
                ch_versions = ch_versions.mix(KALLISTO_QUANT.out.versions)

                PIZZLY( KALLISTO_QUANT.out.txt, transcript, gtf )
                ch_versions = ch_versions.mix(PIZZLY.out.versions)

                GET_PATH(PIZZLY.out.fusions)
                ch_pizzly_fusions = GET_PATH.out.file
            }
        }
        else  {
            ch_pizzly_fusions = ch_dummy_file

        }

    emit:
        fusions             = ch_pizzly_fusions
        versions            = ch_versions.ifEmpty(null)
    }

