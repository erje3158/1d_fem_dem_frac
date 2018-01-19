#!/bin/bash

NAME = debug

# Should be where ./Submit_Script is run
LOCAL_DIR = ${pwd}

LD_LIBRARY_PATH="${LIB}:$LD_LIBRARY_PATH"
echo $LD_LIBRARY_PATH

EXE = ${LOCAL_DIR}/source
BIN = ${LOCAL_DIR}/bin
LIB = ${LOCAL_DIR}/lib
INP = ${LOCAL_DIR}/inputs
OUT = ${LOCAL_DIR}/outputs

RUN_DIR = ${OUT}/${NAME}
mkdir -p ${RUN_DIR}
cp ${EXE}/hu_code ${RUN_DIR}
cd ${RUN_DIR}

pwd
#./hu_code ${INP}/input_boundary_file ${INP}/input_particle_file ${BIN}/qdelaunay . ${INP}/test_fem_input_2el ${INP}/test_dem_input