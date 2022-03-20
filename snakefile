    
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
config_samples = config["batch_parameters"][config_batch]["samples"]


config_database_glob_read = glob.glob(config_database_glob)


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
    input: expand("output/{config_batch}/database.what", config_batch = config_batch)

rule database:
    input: glob.glob(config_database_glob)
    output: "output/{config_batch}/database.what"
    conda: "envs/java.yaml"
    threads: 8
    shell: """

        touch {output}

        # Næste trin: få det til at køre på slurm. 
        """







print("*/")