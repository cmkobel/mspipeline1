    
# alias snakeslurm='mkdir -p logs/old; mv logs/*.{err,out} logs/old 2> /dev/null; snakemake --profile configs/slurm'




import os
from datetime import datetime
import time
import re
#from shutil import copyfile
import glob

import pandas as pd
import re


print("/*")
print("                                                                   ")
print("                                                                   ")
print("                             ___......__             _             ")
print("                         _.-'           ~-_       _.=a~~-_         ")
print(" --=====-.-.-_----------~   .--.       _   -.__.-~ ( ___===>       ")
print("               '''--...__  (    \\ \\\\\\ { )       _.-~           ")
print("                         =_ ~_  \\\\-~~~//~~~~-=-~                 ")
print("                          |-=-~_ \\\\   \\\\                       ")
print("                          |_/   =. )   ~}                          ")
print("                          |}      ||                               ")
print("                         //       ||                               ")
print("                       _//        {{                               ")
print("                    '='~'          \\\\_    =                      ")
print("                                    ~~'                            ")
print("                                                                   ")



# Import configuration
configfile: "config.yaml"
config_batch = config["batch"]
config_d_files = config["batch_parameters"][config_batch]["d_files"]
config_database_glob = config["batch_parameters"][config_batch]["database_glob"]
config_database_glob_read = glob.glob(config_database_glob)
config_samples = config["batch_parameters"][config_batch]["samples"]



# Present configuration
print(f"config_batch:         '{config_batch}'")
print(f"config_d_files:       '{config_d_files}'")
print(f"config_database_glob: '{config_database_glob}'")
for i, j in enumerate(config_database_glob_read):
    print(f"  {i+1}) {j}")
print(f"config_samples:")
for i, j in config_samples.items():
    print(f"  {i}: {j}")
print()


# Define default workflow
rule all:
    input: expand(["output/{config_batch}/database/philosopher_database.fas", \
                   "output/{config_batch}/msfragger/output.what"], \
                   config_batch = config_batch)

rule database:
    input: glob.glob(config_database_glob)
    output: "output/{config_batch}/database/philosopher_database.fas"
    #benchmark: "output/{config_batch}/benchmarks/database.tab"
    threads: 8
    params:
        philosopher = config["philosopher_executable"],
        msfragger = config["msfragger_jar"]
    shell: """



        # Cat all database source files into one.
        cat {input} > output/{config_batch}/database/cat_database_sources.faa



        #touch {output}

        # As philosopher can't specify output files, we need to change dir.
        cd output/{config_batch}/database


        {params.philosopher} workspace \
            --nocheck \
            --clean 

        {params.philosopher} workspace \
            --nocheck \
            --init 

        rm *.fas || echo "nothing to delete" # Remove all previous databases if any.
        {params.philosopher} database \
            --custom cat_database_sources.faa \
            --nodecoys #temp speedup
            #--contam  #temp speedup

        # Manually rename the philosopher output
        mv *.fas philosopher_database.fas

        """




rule msfragger:
    input: "output/{config_batch}/database/philosopher_database.fas"
    output: "output/{config_batch}/msfragger/output.what"
    threads: 8
    params:
        config_d_files = config_d_files
    conda: "envs/openjdk.yaml"
    shell: """
        java \
            -Xmx65G \
            -jar /cluster/projects/nn9864k/shared/bin/MSFrager/bin/MSFragger-3.2/MSFragger-3.2.jar \
            --num_threads {threads} \
            --database_name {input}  {params.config_d_files}/*d


        touch {output}
        """







print("*/")