pipeline {
    agent any

    environment {
        JENKINS_API = credentials('api')
        NXF_VER = 0.32.0
    }

    stages {
        stage('Setup environment') {
            steps {
                sh "pip install nf-core"
                sh "docker pull nfcore/rnafusion:1.0.2"
            }
        }
        stage('Lint markdown') {
            steps {
                sh "markdownlint $WORKSPACE -c $WORKSPACE/.github/markdownlint.yml"
            }
        }
        stage('Nextflow legacy build') {
            steps {
                // sh "nextflow run kraken,jenkins nf-core/rnafusion"
                sh "nextflow run nf-core/rnafusion -r 1.0.2 --help"
                sh "nextflow run nf-core/rnafusion/download-references.nf -r 1.0.2 --help"
                sh "nextflow run nf-core/rnafusion/download-singularity-img.nf -r 1.0.2 --help"
            }
        }
        stage('Nextflow latest build') {
            steps {
                // sh "nextflow run kraken,jenkins nf-core/rnafusion"
                sh "NXF_VER='' nextflow run nf-core/rnafusion -r 1.0.2 --help"
                sh "NXF_VER='' nextflow run nf-core/rnafusion/download-references.nf -r 1.0.2 --help"
                sh "NXF_VER='' nextflow run nf-core/rnafusion/download-singularity-img.nf -r 1.0.2 --help"
            }
        }
    }
    
    post {
        failure {
            script {
                def response = sh(script: "curl -u ${JENKINS_API_USR}:${JENKINS_API_PSW} ${BUILD_URL}/consoleText", returnStdout: true).trim().replace('\n', '<br>')
                def comment = pullRequest.comment("## :rotating_light: Build log output:<br><summary><details>${response}</details></summary>")
            }
        }
    }
}