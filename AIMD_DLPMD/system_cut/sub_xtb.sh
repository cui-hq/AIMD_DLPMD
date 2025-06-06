#!/bin/bash
for dir in $(pwd)/*
do
  if [ -d ${dir} ]; then
     cd ${dir}
     sh xtb.sh > xtb.out
  fi
done
