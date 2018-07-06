#!/bin/bash

# Currently can only be run from ./source

if [ "$#" = 1 ]
	then
	if [ "$1" = "excalibur" ]
		then
			echo "Compiling Hierarchical Upscaling Code on excalibur.arl.hpc.mil"
			cp excalibur/Makefile ./

			module swap PrgEnv-intel PrgEnv-gnu
			module list

			cd ../lib
			export LD_LIBRARY_PATH="$(pwd):$LD_LIBRARY_PATH"
			echo $LD_LIBRARY_PATH
			cd ../source

			rm *.o *~ hu_code
			make
			echo "Done!"

			rm Makefile 

	elif [ "$1" = "topaz" ]
		then
			echo "Compiling Hierarchical Upscaling Code on topaz.erdc.hpc.mil"
			cp topaz/Makefile ./

			module unload compiler/intel/16.0.0
			module load compiler/gcc/6.1.0
			module list

			cd ../lib
			MKLA=/p/home/apps/intel/parallel_studio_2017_u1/compilers_and_libraries_2017.1.132/linux/compiler/lib/intel64_lin/
			MKLB=/p/home/apps/intel/parallel_studio_2016/compilers_and_libraries_2016.0.109/linux/mkl/lib/intel64_lin/
			export LD_LIBRARY_PATH="$(pwd):$MKLA:$MKLB:$LD_LIBRARY_PATH"
			echo $LD_LIBRARY_PATH
			cd ../source

			rm *.o *~ hu_code
			make
			echo "Done!"

			rm Makefile

	elif [ "$1" = "copper" ]
        then
            echo "Compiling Hierarchical Upscaling Code on copper.ors.hpc.mil"
            cp copper/Makefile ./

            module swap PrgEnv-pgi PrgEnv-gnu
            module list

            cd ../lib
            MKLA=/opt/intel/lib/intel64/
            MKLB=/opt/intel/mkl/lib/intel64/
            export LD_LIBRARY_PATH="$(pwd):$MKLA:$MKLB:$LD_LIBRARY_PATH"
            echo $LD_LIBRARY_PATH
            cd ../source

            rm *.o *~ hu_code
            make
            echo "Done!"

            rm Makefile

	elif [ "$1" = "macOS" ]
		then
			echo "Compiling Hierarchical Upscaling Code on macOS"
			cp macOS/Makefile ./

			cd ../lib
			export LD_LIBRARY_PATH="$(pwd):$LD_LIBRARY_PATH"
			echo $LD_LIBRARY_PATH
			cd ../source

			rm *.o *~ hu_code
			make
			echo "Done!"

			rm Makefile

	elif [ "$1" = "soilblast" ]
		then
			echo "Compiling Hierarchical Upscaling Code on soilblast2.colorado.edu"
			cp soilblast/Makefile ./

			module load gcc-6.4.0
			module list

			cd ../lib
			export LD_LIBRARY_PATH="$(pwd):$LD_LIBRARY_PATH"
			echo $LD_LIBRARY_PATH
			cd ../source

			rm *.o *~ hu_code
			make
			echo "Done!"

			rm Makefile 

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
