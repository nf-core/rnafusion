pipeline {
    agent any

    environment {
        NXF_VER = '0.32.0'
    }

    stages {
        stage('Setup environment') {
            steps {
                sh "pip install nf-core"
                sh "pip install pylint"
                sh "docker pull nfcore/rnafusion:dev"
                sh "docker tag nfcore/rnafusion:dev nfcore/rnafusion:1.0"
            }
        }

        stage('Lint markdown') {
            steps {
                sh "markdownlint $WORKSPACE/$JOB_NAME -c $WORKSPACE/$JOB_NAME/.github/markdownlint.yml"
            }
        }
        stage('Lint python code') {
            steps {
                sh "pylint --rcfile=$WORKSPACE/$JOB_NAME/.github/pylintrc $WORKSPACE/$JOB_NAME/bin/*/*.py --ignore=scrape_software_versions.py"
            }
        }
        /*
        stage('Build') {
            steps {
                nextflow run . -profile test,docker
            }
        }*/
    }
}