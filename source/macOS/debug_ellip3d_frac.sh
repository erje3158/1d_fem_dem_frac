#!/bin/bash

NAME="debug_ellip3d_frac"

# Should be where ./Submit_Script is run
LOCAL_DIR=$(pwd)

EXE=${LOCAL_DIR}/source
BIN=${LOCAL_DIR}/bin
LIB=${LOCAL_DIR}/lib
INP=${LOCAL_DIR}/inputs/
DEM=${LOCAL_DIR}/inputs/dem_assemblies/fracture
OUT=${LOCAL_DIR}/outputs

LD_LIBRARY_PATH="${LIB}:$LD_LIBRARY_PATH"

RUN_DIR=${OUT}/${NAME}

if [ -d "$RUN_DIR" ]
	then
		echo "Removing previous run directory: "$RUN_DIR
		rm -rf $RUN_DIR
		pwd
		ls $OUT/
fi

mkdir -p ${RUN_DIR}
cp ${EXE}/hu_code ${RUN_DIR}
cd ${RUN_DIR}

echo
echo "Output Directory : "
pwd

echo
echo
echo "Running Hierarchical Upscaling Code"
echo "-----------------------------------"
./hu_code ${DEM}/input_boundary_file ${DEM}/input_particle_file ${BIN}/qdelaunay . \
          ${INP}/debug_ellip3d_frac_fem ${INP}/debug_ellip3d_frac_dem
echo "-----------------------------------"
echo "End of Simulation"