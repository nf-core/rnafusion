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
                sh "markdownlint $WORKSPACE -c $WORKSPACE/.github/markdownlint.yml"
            }
        }
        stage('Lint python code') {
            steps {
                sh "pylint --rcfile=$WORKSPACE/.github/pylintrc $WORKSPACE/bin/*/*.py --ignore=scrape_software_versions.py"
            }
        }
    }
    post {
        always {
            sh "pwd && cat ${JENKINS_HOME}/jobs/${JOB_NAME}/builds/${BUILD_NUMBER}/log >> log.txt"
            deleteDir()
        }
    }
}