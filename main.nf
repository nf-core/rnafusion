#!/usr/bin/env nextflow
/*
========================================================================================
    nf-core/rnafusion
========================================================================================
    Github : https://github.com/nf-core/rnafusion
    Website: https://nf-co.re/rnafusion
    Slack  : https://nfcore.slack.com/channels/rnafusion
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
========================================================================================
*/

WorkflowMain.initialise(workflow, params, log)

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/


//
// WORKFLOW: Run main nf-core/rnafusion analysis pipeline
//
workflow NFCORE_RNAFUSION {

    if (params.build_references) {

        include { BUILD_REFERENCES } from './workflows/build_references'
        BUILD_REFERENCES ()
    
    } else {

        include { RNAFUSION } from './workflows/rnafusion'
        // RNAFUSION ()

    }

}

/*
========================================================================================
    RUN ALL WORKFLOWS
========================================================================================
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    NFCORE_RNAFUSION ()
}

/*
========================================================================================
    THE END
========================================================================================
*/
