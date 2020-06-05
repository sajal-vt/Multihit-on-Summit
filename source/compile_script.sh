#!/bin/bash
module load cuda/10.1.243
export MPI_PATH=/autofs/nccs-svm1_sw/summit/.swci/1-compute/opt/spack/20180914/linux-rhel7-ppc64le/xl-16.1.1-3/spectrum-mpi-10.3.0.1-20190611-aqjt3jo53mogrrhcrd2iufr435azcaha/
module load spectrum-mpi/10.3.0.1-20190611
#nvcc -I/home/sajal/.openmpi/include -L/home/sajal/.openmpi/lib -lmpi -code=compute_53 -arch=compute_53 3hit_mpi.cu -o solve
nvcc -I$MPI_PATH/include -L$MPI_PATH/lib -use_fast_math -code=compute_53 -arch=compute_53 -lmpi_ibm $1 -o $2
