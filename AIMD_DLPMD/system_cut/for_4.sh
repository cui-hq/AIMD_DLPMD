#!/bin/bash
source ~/.bashrc
for dir in $(pwd)/*
do
  if [ -d ${dir} ]; then
     cd ${dir}
     if [ -f xtbopt.log ]; then
	     tail -n $(head -n 1 coord_num.xyz) xtbopt.log > coord.xyz
     else
	     tail -n $(head -n 1 coord_num.xyz) coord_num.xyz > coord.xyz
     fi
     rm scoord*
     rm nohup.out
     rm charges
     rm wbo
     rm *log
     rm xtb*
     rm mdrestart
  fi
done
