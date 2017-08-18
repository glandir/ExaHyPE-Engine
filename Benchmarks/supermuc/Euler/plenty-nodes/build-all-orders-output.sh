directory=plenty-nodes

exe=ExaHyPE-Euler
spec=$directory/Euler-output.exahype

# save originial file
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

  for p in 9 7 5 3
  do 
    rm *.o
    sed -i -r 's,order(\s+)const(\s+)=(\s+)([0-9]+),order\1const\2=\3'$p',' $spec
    cat $spec
    $directory/configure-output.sh
    make -j28 && \
    mv $exe $exe-p$p-$SHAREDMEM-$COMPILER
  done
done

# restore original file
mv ${spec}_tmp $spec
