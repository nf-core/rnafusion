pipeline {
    agent any

    environment {
        NXF_VER = '0.32.0'
        JENKINS_API = credentials('api')
    }

    stages {
        stage('Setup environment') {
            steps {
                sh "pip install nf-core"
                sh "docker pull nfcore/rnafusion:dev"
                sh "docker tag nfcore/rnafusion:dev nfcore/rnafusion:1.0.1"
            }
        }
        stage('Lint markdown') {
            steps {
                sh "markdownlint $WORKSPACE -c $WORKSPACE/.github/markdownlint.yml"
            }
        }
        stage('Build') {
            steps {
                // sh "nextflow run kraken,jenkins nf-core/rnafusion"
                nextflow run nf-core/rnafusion --help
                nextflow run nf-core/rnafusion/download-references.nf -profile jenkins --help
                nextflow run nf-core/rnafusion/download-singularity-img.nf -profile jenkins --help
            }
        }
    }
    post {
        failure {
            script {
                def response = sh(script: "curl -u ${JENKINS_API_USR}:${JENKINS_API_PSW} ${BUILD_URL}/consoleText", returnStdout: true).trim().replace('\n', '<br>')
                def comment = pullRequest.comment("##:rotating_light: Buil log output:<br><summary><details>${response}</details></summary>")
            }
        }
    }
}