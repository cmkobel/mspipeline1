#!/bin/bash

# Variables
wait_seconds=2
runid=$1
outfile="memory_usage_mebi_${runid}.txt"

# Presentation
echo "# Writing memory usage every ${wait_seconds} seconds to file ${outfile}"





# TODO: change mega into mibi. Also in the filename

while true; do
    
    free --mebi -w \
    | awk -v date="$(echo -n "$runid ";  date -Iseconds)" '{ print date "           " $0 }' \
    | grep "Mem:" \
    >> $outfile ;
    echo -n .
    sleep $wait_seconds;

done



