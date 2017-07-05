exe=ExaHyPE-Euler
headerADERDG=AbstractEulerSolver_ADERDG.h
headerFV=AbstractEulerSolver_FV.h

for m in 1 2
do
  if (( m == 1 )); then
    make clean
    export SHAREDMEM=TBB
  else
    make clean
    export SHAREDMEM=None
  fi

  echo "SHAREDMEM=$SHAREDMEM"
  #read -p "press any key..."

  for p in 1 3 5 7 9
  do
    let N=2*p+1
    
    rm *.o
    sed -i -r "s,Order(\s+)= ([0-9]),Order\1= ${p}," $headerADERDG 
    sed -i -r "s,PatchSize(\s+)= ([0-9]+),PatchSize\1= ${N}," $headerFV
    
    make -j28 && \
    mv $exe $exe-p$p-$SHAREDMEM-$COMPILER
  done
done

