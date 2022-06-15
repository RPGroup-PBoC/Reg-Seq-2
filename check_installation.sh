#!/bin/bash

# Check if fastp is installed
conda list -n wgregseq fastp

# Check if bbmap is in the correct location
FILE=./bbmap/bbmap.sh
if test -f "$FILE"; then
    echo "$FILE exists."
else
    echo "$FILE does not exist! Please make sure it is at the correct location."
fi


