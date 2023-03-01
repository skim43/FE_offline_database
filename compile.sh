#!/bin/bash
  
#old setup
#module load cmake
module load mvapich2/2.3.1/gcc/8.3.0
source /opt/crc/i/intel/18.0/mkl/bin/mklvars.sh intel64 lp64

export branch_name=`git branch | grep '*' | awk '{ print $2; }'` #To get the name of the current branch
export PGFEM3D_INSTALL=$PWD/build_isir_$branch_name                   #Give the name of the build directory 

make distclean
./bootstrap

./configure --prefix=$PGFEM3D_INSTALL CC=mpicc CXX=mpicxx\
 CXXFLAGS="-O1 -g" \
 MSNET_NET=isir \
--with-cnstvm=/afs/crc.nd.edu/user/s/skim43/repos_test/Generalizsed_constitutive_model \
--with-hypre=/afs/crc.nd.edu/group/cswarm/hypre/2.15.1/gcc/8.3.0/mvapich2/2.3.1 \
--with-suitesparse=/afs/crc.nd.edu/group/cswarm/SuiteSparse/4.5.5/gcc/8.3.0 \
--with-ttl=/afs/crc.nd.edu/group/cswarm/ttl/install_ttl_05182019 \
--enable-vtk=yes\
--with-vtk=/opt/crc/v/vtk/8.2.0/gcc/4.8.5\
--with-vtk-version=-8.2\
--with-tests-nprocs=16

make -j 4
make install