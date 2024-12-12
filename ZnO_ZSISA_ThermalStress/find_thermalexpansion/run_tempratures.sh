#!/bin/bash

# List of specific temperatures to process
temperatures=("0000" "0025" "0050" "0075" "0100" "0125" "0150" "0200" "0250" "0300" "0350" "0400" "0450" "0500" "0550" "0600" "0700" "0800" "0900" "1000")
pressure="00" 

# Initialize previous directory GSR file
prev_dir=""
prev_gsr_file="Relax2o_GSR.nc"

# Loop through the specified temperatures
for temp in "${temperatures[@]}"; do
    dir_name="Temp_${temp}_${pressure}"
    
    # Create the directory and navigate into it
    mkdir $dir_name
    cd $dir_name

    # Copy the previous Relax2o_GSR.nc file if it's not the first iteration
    if [ -n "$prev_dir" ]; then
        if [ ! -f "$prev_gsr_file" ]; then
            cp "../$prev_dir/$prev_gsr_file" .
        fi
    fi

    # Run the Python script with temperature and pressure arguments
    python ../workflow.py $temp $pressure

    # Set the current directory as the previous directory for the next iteration
    prev_dir=$dir_name

    # Navigate back to the parent directory
    cd ..
done

