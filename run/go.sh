#!/bin/sh

#PBS -q u-lecture
#PBS -N omp
#PBS -l select=1:ncpus=18
#PBS -Wgroup_list=gt21
#PBS -l walltime=00:10:00
#PBS -e err
#PBS -o test18.lst

#cd $PBS_O_WORKDIR
. /etc/profile.d/modules.sh

export OMP_NUM_THREADS=18
export KMP_AFFINITY=granularity=fine,compact,1,0
numactl ./sol20
