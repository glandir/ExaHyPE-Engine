exe=ExaHyPE-Euler
spec=multicore/Euler-output.exahype

cp $spec ${spec}_tmp

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

  for p in 3 5 7 9
  do 
    rm *.o
    sed -i -r 's,order(\s+)const(\s+)=(\s+)([0-9]+),order\1const\2=\3'$p',' $spec
    cat $spec
    multicore/configure-output.sh
    make -j28 && \
    mv $exe $exe-p$p-$SHAREDMEM-$COMPILER
  done
done

mv ${spec}_tmp $spec
