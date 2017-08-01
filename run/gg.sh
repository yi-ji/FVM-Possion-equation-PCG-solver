#!/bin/sh

#PBS -q u-lecture
#PBS -N omp
#PBS -l select=1:ncpus=16
#PBS -Wgroup_list=gt00
#PBS -l walltime=00:10:00
#PBS -e test16.err
#PBS -o test16.lst

cd $PBS_O_WORKDIR
. /etc/profile.d/modules.sh

export OMP_NUM_THREADS=16
export KMP_AFFINITY=granularity=fine,compact,1,0
numactl ./sol20
numactl ./sol20
numactl ./sol20
numactl ./sol20
numactl ./sol20
