pipeline {
    agent any

    environment {
        NXF_VER = '0.32.0'
        REPO_NAME = 'rnafusion'
        REPO_URL = 'https://github.com/nf-core/rnafusion.git'
    }

    stages {
        stage('Checkout') {
            steps {
                echo `pwd`
                echo "BUILD_NUMBER"
                echo "WORKSPACE"
            }
        }
    }
}