#!/bin/bash
for dir in $(pwd)/*
do
  if [ -d ${dir} ]; then
     cd ${dir}
     mpirun -np 64 cp2k.popt -i input_nvt.inp -o input_nvt.log
  fi
done
