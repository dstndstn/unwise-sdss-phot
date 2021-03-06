#! /bin/bash

unset PBS_ARRAYID
unset PBS_O_WORKDIR

# respect memory limits... 4GB
#ulimit -v 4194304
# 16 GB
ulimit -v 16777216

cd ${HOME}/wphot
#while [ $# -gt 0 ]; do
  n="$1"
  echo "Running $n"
  python -u forcedphot.py -B 8 --dataset sdss-dr10 --split -v -b 1234 \
      --photoobjsdir /clusterfs/riemann/raid006/dr10/boss/photoObj \
      --resolvedir resolve-2010-05-23-cut \
      $n > sdss-dr10-logs/$n.log 2> sdss-dr10-logs/$n.err
#  shift
#done
