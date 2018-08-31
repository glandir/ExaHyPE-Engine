#!/usr/bin/env groovy
util = load "Miscellaneous/Jenkins/Util.groovy"

def build(config, workspace) {
    node ('Linux-Cluster') {
	ws("${env.JOB_NAME}-${config.name}-build") {
	    deleteDir()
	    def directory=sh(returnStdout: true,
			     script: "dirname ${config.exahypeFile}"
	    ).trim()
	    echo "Directory is $directory"
util.slurmBatch '''
set -x
IFS=$'\n\t'
pwd
''' + util.getModuleCode() + """
set -euo pipefail
cp -rf ${workspace}/. .
path=${config.exahypeFile}
""" + '''
ls

# Fix linker flags
export COMPILER_LFLAGS='-pthread'
export PROJECT_LFLAGS="-lrt" 

# Build toolkit
Toolkit/build.sh

# Build example
echo "Building $path"
dir="$(readlink -f $(dirname ${path}))"

java -jar Toolkit/dist/ExaHyPE.jar --not-interactive ${path}
cd $dir
make -j 32
'''
	    // Stash files for later reuse
	    // This is needed because the run step could run on another node!
	    stash includes: "${directory}/**", name: "exahype-${config.name}"
	    deleteDir()}
    }

}

def assembleBuildMatrix(combinations, workspace) {
    def buildMatrix = combinations.collectEntries {
	["${it.name}" : {
	    -> build(it, "${workspace}")
	}]
    }
    return buildMatrix
}


return this