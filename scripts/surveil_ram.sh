#!/bin/bash

# Variables
wait_seconds=1
runid=$1
outfile="memory_usage_mibi_${runid}.txt"

# Presentation
echo "# Writing memory usage every ${wait_seconds} seconds to file ${outfile}"







while true; do
    
    free --mega -w \
    | awk -v date="$(echo -n "$runid ";  date -Iseconds)" '{ print date "           " $0 }' \
    >> $outfile ;
    echo -n .
    sleep $wait_seconds;

done



