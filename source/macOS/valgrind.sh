#!/bin/bash

NAME="valgrind"

# Should be where ./Submit_Script is run
LOCAL_DIR=$(pwd)

EXE=${LOCAL_DIR}/source
BIN=${LOCAL_DIR}/bin
LIB=${LOCAL_DIR}/lib
INP=${LOCAL_DIR}/inputs
OUT=${LOCAL_DIR}/outputs
SRC=${LOCAL_DIR}/source

LD_LIBRARY_PATH="${LIB}:$LD_LIBRARY_PATH"

RUN_DIR=${OUT}/${NAME}

if [ -d "$RUN_DIR" ]
	then
		echo "Removing previous run directory: "$RUN_DIR
		rm -rf $RUN_DIR
fi

mkdir -p ${RUN_DIR}
cp ${EXE}/hu_code ${RUN_DIR}
cp ${SRC}/macOS/suppressions ${RUN_DIR}
cd ${RUN_DIR}

echo
echo "Output Directory : "
pwd

echo
echo
echo "Running Hierarchical Upscaling Code"
echo "-----------------------------------"
valgrind --num-callers=100 --leak-check=full --show-leak-kinds=all --gen-suppressions=all \
         --suppressions=suppressions --track-origins=yes -v --log-file=val_out.txt \
         ./hu_code ${INP}/input_boundary_file ${INP}/input_particle_file ${BIN}/qdelaunay . \
         ${INP}/valgrind_fem ${INP}/valgrind_dem
echo "-----------------------------------"
echo "End of Simulation"