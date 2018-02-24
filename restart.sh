#!/bin/bash
module load cmake/3.8.0
#module load gcc/4.8.2
module load gcc/5.4.0
module load vtk
#module load openmpi/1.8.4__gcc-4.9.2
#module load gpu/cuda/8.0.44
#module load gpu/cuda/7.5.18
export BOOST_ROOT=/system/software/linux-x86_64/lib/boost/1_56_0
#export BOOST_ROOT=/data/coml-aboria/robinsonm/boost_1_58
#export EIGEN_ROOT=/system/software/linux-x86_64/lib/eigen/3.2.8/include/eigen3
export EIGEN_ROOT=/data/coml-aboria/robinsonm/eigen-eigen-5a0156e40feb/include
#export CUDA_INCLUDE_DIRS=/system/software/arcus-b/gpu/cuda/8.0.44/include
#export CUDA_INCLUDE_DIRS=/system/software/arcus-b/gpu/cuda/7.5.18/include
#export CGAL_ROOT=/data/coml-aboria/robinsonm/CGAL-4.10/install

cmake \
    #-DCGAL_ROOT=$CGAL_ROOT \
    -DCMAKE_BUILD_TYPE=RELEASE \
    -DAboria_USE_OPENMP=ON \
    -DBoost_NO_SYSTEM_PATHS=BOOL:ON \
    -DBoost_NO_BOOST_CMAKE=BOOL:ON \
    -DBOOST_LIBRARYDIR=$BOOST_ROOT/lib \
    -DBOOST_INCLUDEDIR=$BOOST_ROOT/include \
    -DEIGEN3_INCLUDE_DIR=$EIGEN_ROOT \
    ..


