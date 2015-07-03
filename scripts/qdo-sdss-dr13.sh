#! /bin/bash

unset PBS_ARRAYID
unset PBS_O_WORKDIR

# respect memory limits... 4GB
#ulimit -v 4194304
# 16 GB
ulimit -v 16777216

cd ${HOME}/wphot

n="$1"
echo "Running $n"
python -u forcedphot.py -B 8 --dataset sdss-dr13 --split -v -b 1234 \
    --photoobjsdir /clusterfs/riemann/raid007/ebosswork/eboss/photoObj.v5b \
    --resolvedir resolve-2013-07-29 \
    $n > sdss-dr13-logs/$n.log 2> sdss-dr13-logs/$n.err
