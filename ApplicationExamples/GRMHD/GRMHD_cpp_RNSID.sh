#!/bin/bash
#
# Generate the combinations from a mexafile

MEXA="exa mexa"
infile="GRMHD_cpp_RNSID.mexa"

mkdir -p Specfiles

verbose() { echo $@; $@; } 

# todo: Obtain the combinations from the mexa file itself
#       and just iterate instead of doing this here.

for solver in DG FV Limiting; do
	solverkey="${solver}Solver"
	dim=2
		
	outfile="Specfiles/${infile/.*}-${solver}-${dim}D.exahype"
	tmpfile=$(mktemp)
	
	echo "Making $outfile ..."

	verbose $MEXA specfile  --parameter-style=adapted \
		--infile $infile --outfile $tmpfile \
		--add Solver=$solverkey
		
	if [ -e $outfile ] && ! diff $outfile $tmpfile; then
		echo "$outfile has local changes, is not overwritten. Delete to recreate."
	else
		mv $tmpfile $outfile
		#cat <(echo "/* Generated by $0 from $infile at $(date) */") $tmpfile > $outfile
	fi
done