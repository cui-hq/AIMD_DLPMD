#!/bin/bash
for dir in $(pwd)/*
do
  if [ -d ${dir} ]; then
     cd ${dir}
     cp ../../trans.py .
     python trans.py
  fi
done
