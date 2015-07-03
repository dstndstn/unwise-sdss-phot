#! /bin/bash

unset PBS_ARRAYID
unset PBS_O_WORKDIR

cd ${HOME}/wphot
while [ $# -gt 0 ]; do
  n="$1"
  echo "Running $n"
  python -u forcedphot.py -B 8 --dataset sdssv4 --split -v -b 1234 $n > sdssv4-logs/$n.log 2> sdssv4-logs/$n.err
  shift
done
