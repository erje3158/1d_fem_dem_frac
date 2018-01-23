#!/bin/bash

if [ "$#" = 2 ]
	then
	if [ "$1" = "excalibur" ]
		then
			echo "Running Hierarchical Upscaling Code on excalibur.arl.hpc.mil"
			cd ./source/excalibur/

			if [ -f "$2" ]
			then
				cp ./"$2" ../../
				cd ../../

				qsub "$2"
				rm "$2"

			else
				echo "--> $2 is not the name of a currently existing PBS script"
				echo "--> List of available scripts: ./source/<platform>"
				ls -l *.pbs
			fi

	elif [ "$1" = "topaz" ]
		then
			echo "Running Hierarchical Upscaling Code on topaz.erdc.hpc.mil"
			cd ./source/topaz/

			if [ -f "$2" ]
			then
				cp ./"$2" ../../
				cd ../../
				
				qsub "$2"
				rm "$2"
				
			else
				echo "--> $2 is not the name of a currently existing PBS script"
				echo "--> List of available scripts: ./source/<platform>"
				ls -l *.pbs
			fi

	elif [ "$1" = "macOS" ]
		then
			echo "Running Hierarchical Upscaling Code on macOS"
			cd ./source/macOS/

			if [ -f "$2" ]
			then
				cp ./"$2" ../../
				cd ../../

				./"$2" 
				rm "$2"

			else
				echo "--> $2 is not the name of a currently existing submit script"
				echo "--> List of available scripts: ./source/<platform>"
				echo ls -l *.sh
			fi

	if [ "$1" = "soilblast" ]
		then
			echo "Running Hierarchical Upscaling Code on soilblast2.colorado.edu"
			cd ./source/soilblast/

			if [ -f "$2" ]
			then
				cp ./"$2" ../../
				cd ../../

				qsub "$2"
				rm "$2"

			else
				echo "--> $2 is not the name of a currently existing PBS script"
				echo "--> List of available scripts: ./source/<platform>"
				ls -l *.pbs
			fi
	else
		echo "--> Don't recognize that platform..."
		echo "--> Need to specify recognized platform: excalibur or topaz"
		echo "--> ./Submit_Script.sh <platform> <PBS script>"
	fi
else
	echo "--> Incorrect number of arguments - Need to specify the machine where the code is being compiled and name of an exisiting PBS script."
	echo "--> ./Submit_Script.sh <platform> <PBS script>"
	echo "--> Current platforms: excalibur, topaz"
	echo "--> PBS scripts are stored in ./source/<platform> - you will find some examples already in that directory"
	echo "--> Using the examples as a template, generate a PBS script for your job, store it in the correct directory, and run the code by specifying it's name when calling Submit_Script.sh"
fi
