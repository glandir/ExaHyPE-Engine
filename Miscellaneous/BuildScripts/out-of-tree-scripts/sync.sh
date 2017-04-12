#!/bin/bash
# Mini build script wrapper for ExaHyPE out of tree build system.
# Autogenerated by setup-out-of-tree.sh
# @AUTOGENINFO@
set -e
cd "$(dirname "$0")"
source "oot.env"

# -L: Copy symlinks as real files
rsync="rsync -r -L --relative --delete --exclude=.git/ --exclude=.svn/"

# assume that $oot_appdir can be quite populated with a mess
# caveats, if the variables are empty, there is --exclude=*
rsync_app_excludelist="--exclude=*.log --exclude=*.vtk --exclude=*.csv --exclude=*.mk --exclude=${oot_makesystem_prefix}* --exclude=${oot_buildprefix}*/"

echo "Build $oot_buildname: Copying files to $oot_builddir ..."

oot_codedir="$(readlink -f $oot_codedir)" # absolutize with readlink
#oot_codedir="$(pwd)/$oot_codedir/" # absolutize without readlink
cd $oot_abs_exahype # in order rsync --relative works

$rsync $oot_toolkit $oot_codedir/
$rsync $rsync_app_excludelist $oot_appdir $oot_codedir/

# before:
#cp $oot_abs_exahype/$oot_specfile @oot_codedir@/$oot_appdir/.. # the path is hacky

# after:
$rsync $oot_specfile $oot_codedir/
$rsync $oot_dependencies $oot_codedir/

# Collect Repository information for in-build --version information
buildinfo_script="./ExaHyPE/generate-buildinfo.sh"

# extract the filename line aka static_repo_info="static-repo-info.txt"
eval $(grep -E "^static_repo_info=" $buildinfo_script)
[[ x${static_repo_info} == x ]] && { echo "Failure: Somebody removed the static_repo_info definition."; }

# execute the buildinfo extraction
EXAHYPE_PATH="./ExaHyPE"
PEANO_KERNEL_PEANO_PATH="./Peano/peano"
$buildinfo_script "EXAHYPE_PATH=$EXAHYPE_PATH" "PEANO_KERNEL_PEANO_PATH=$PEANO_KERNEL_PEANO_PATH" "SYNC_OOT_FAKECALL=YES" \
	| tee $oot_codedir/$EXAHYPE_PATH/$static_repo_info \
	| tee $oot_codedir/$PEANO_KERNEL_PEANO_PATH/$static_repo_info \
	> /dev/null

echo "Finished, oot_builddir size is $(du -hs $oot_builddir | head -n1 | awk '{print $1}') in $(find $oot_builddir | wc -l) files"

