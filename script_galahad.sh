#!/bin/bash
filename='benchmark/list_problems_1000000.dat'
n=1
while read line;
do
  # for read each line
  echo "No. $n : $line"
  sifdecoder ROSENBR
  runcutest -p gen90 ROSENBR
  n=$((n+1))
done < $filename
