#!/bin/bash

# Currently can only be run from ./source

if [ "$#" = 1 ]
	then
	if [ "$1" = "excalibur" ]
		then
			echo "Compiling Hierarchical Upscaling Code on excalibur.arl.hpc.mil"

			module swap PrgEnv-intel PrgEnv-gnu
			module list

			rm *.o *~ ellip3d
			make
			echo "Done!"

	else
		echo "Don't recognize that platform..."
		echo "Need to specify recognized platform: excalibur or topaz"
		echo "./Make_Script.sh <platform>"
	fi
else
	echo "Incorrect number of arguments - Need to specify the machine where the code is being compiled."
	echo "./Make_Script.sh <platform>"
	echo "Current platforms: excalibur, topaz"
fi
