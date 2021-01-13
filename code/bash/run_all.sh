#!/usr/bin/env bash

bin_folder="/Users/Manu/Documents/Projects/TheoryPaper/rigid_crosslinkers/code/cytosim/bin"
python_folder="/Users/Manu/anaconda/lib/python2.7/site-packages/leraramirez2019"
python_bin="/Users/Manu/anaconda/bin/python"
# What is called with help
usage="
=============== Help for 'run_all.sh ==============='

Usage:

    bash run_all.sh -n X -j Y run_directories

-n: specify the number of config files to generate with 'preconfig.py' (default 1)
-j: specify the number of jobs when running/analyzing simulations in parallel (default 1)
run_directories: directories in which to run the pipeline

For this script to work:

    0.1 Set the variables 'bin_folder', 'python_folder' and 'python_bin' in this script (lines 3-5) to the ABSOLUTE paths of the
    folders containing the simulation binaries, the python scripts called in this script, and your Python 2 binary.
    0.2 Make sure that the 'leraramirez2019' python module is in your python path, for the 'get_*.py' folders called by
     'master.py' to work (see the 'README.md' for more details on this).
    binaries.
    1. In each 'run_directory', there must be a template configuration file called 'config.cym.tpl'
    2. In the folder 'run_directory/..' there must be:
        a) A file called 'report_instructions.txt', which contains the report commands
        b) Files with names 'get_*.py', that are used to extract information from the reported files (see 'master.py').

Tasks performed:

    1. Delete previous simulation data that exist in the folders 'run_directory/scan/run????', and those folders.
    2. Create X new simulations folders in the form 'run_directory/scan/run????', where each of them corresponds to a
    simulation, and contains a 'config.cym' file created by 'preconfig.py' from 'run_directory/config.cym.tpl'
    3. Run the simulations, using Y jobs.
    4. Extract information from the simulation folders using 'master.py' (see 'master.py')
    5. Measure in the simulations based on the data using 'measure.py' (see 'measure.py' for each case)

"

# Stop on error
set -e

# Store the current directory, in case the paths of the targets are given as relative paths
current_dir=$(pwd)

# Default value
num_sims=1
nb_process=1

# Get keyword arguments

while getopts ":hj:n:f:" opt; do
  case $opt in
    h)  printf "$usage"

        exit
        ;;

    j)  nb_process="$OPTARG"
        ;;

    n)  num_sims="$OPTARG"
        ;;

    \?) echo "Invalid option -$OPTARG" >&2
        ;;
  esac
done

# Shift to non-keyword arguments (target folders)
shift $((OPTIND-1))

# Store current directory path

for run_directory in "$@"
do

    ## Running simulations ---------------------------------------------------------------------------------------------

    if false
    then
        echo "Running ---------> "${run_directory} "creating ${num_sims} simulations per condition and using ${nb_process} jobs for analysis and running"

        # Delete any previous data (if there was any)
        rm -f ${run_directory}/scan/run*/.DS_Store
        rm -f ${run_directory}/scan/run*/*.cym
        rm -f ${run_directory}/scan/run*/*.cmo
        rm -f ${run_directory}/scan/run*/*.txt
        rm -f ${run_directory}/scan/.DS_Store

        for j in ${run_directory}/scan/run*
        do
            if [ -d $j ]; then rmdir $j; fi
        done
        if [ -d ${run_directory}/scan ]; then rmdir ${run_directory}/scan; fi

        # Make the folder where all the simulation folders will be
        mkdir ${run_directory}/scan

        # Use preconfig.py to make a all the config files, you can change the number to the desired number of simulations

        ${python_bin} ${python_folder}/preconfig.py ${num_sims} ${run_directory}/config.cym.tpl scan

        # Use collect to put all the config files in subfolders with formatted names

        ${python_bin} ${python_folder}/collect.py ${run_directory}/scan/run%04i/config.cym scan/config????.cym

        # Run simulations in parallel

        ${python_bin} ${python_folder}/scan.py ${bin_folder}/sim nproc=${nb_process} ${run_directory}/scan/run????
    fi

	## Extracting information from simulations ------------------------------------------------------------------------
    if true
    then
	# Extract information (see master.py)

	${python_bin} ${python_folder}/master.py --nb_jobs ${nb_process} --bin_folder ${bin_folder} ${run_directory}

	# Move to data_summary directory (created by master.py)

    cd ${run_directory}/../../data_summary

    echo "-----> Measuring"
	${python_bin} measurements.py $(basename ${run_directory})

    # Move back to the original directory

    cd ${current_dir}
    fi
done

