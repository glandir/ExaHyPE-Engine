#!/bin/bash
# Mini build script wrapper for ExaHyPE out of tree build system.
# Autogenerated by setup-out-of-tree.sh
# @AUTOGENINFO@
set -e
cd "$(dirname "$0")"
source "oot.env"

codedir=$(readlink -f $oot_codedir)
rm -r $codedir
rm -rf -- "$(pwd -P)"