    
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
config_d_base = config["batch_parameters"][config_batch]["d_base"]
config_database_glob = config["batch_parameters"][config_batch]["database_glob"]
config_database_glob_read = glob.glob(config_database_glob)
config_samples = config["batch_parameters"][config_batch]["samples"]

# Present configuration
print(f"config_batch:         '{config_batch}'")
print(f"config_d_base:        '{config_d_base}'")
print(f"config_database_glob: '{config_database_glob}:'")
for i, j in enumerate(config_database_glob_read):
    print(f"  {i+1}) {j}")
print()


# Populate dataframe
df = pd.DataFrame(data = {'sample':  config_samples.keys(),
                          'barcode': config_samples.values()})

df["path"] = config_d_base + "/" + df["barcode"]
#pd.DataFrame.from_dict({'sample': config_samples.keys()})

print(df)
print("//")

#print(df["path"].tolist()) # Debug



# Define default workflow
rule all:
    input: expand(["output/{config_batch}/database/philosopher_database.fas", \
                   "output/{config_batch}/msfragger/output.what"], \
                   config_batch = config_batch, \
                   sample = df["sample"])

rule database:
    input: glob.glob(config_database_glob)
    output: "output/{config_batch}/database/philosopher_database.fas"
    #benchmark: "output/{config_batch}/benchmarks/database.tab"
    threads: 8
    params:
        philosopher = config["philosopher_executable"]
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
    input:
        database = ["output/{config_batch}/database/philosopher_database.fas"],
        d_files = df["path"].tolist()
    output: 
        msfragger_version = "output/{config_batch}/msfragger/msfragger_version.txt",
        untouchable = "output/{config_batch}/msfragger/output.what"
    #shadow: "shallow" # Can't use shadow as it isn't a subdir.
    threads: 8
    params:
        config_d_base = config_d_base,
        msfragger_jar = config["msfragger_jar"]
    conda: "envs/openjdk.yaml"
    shell: """

        # We can't manage the output path of MSFragger, so we need to symlink it to the directory where we want the output to reside.
        # And no, we can't use snakemake-shadow, as it only works on relative subdirs.


        java \
            -Xmx64G \
            -jar {params.msfragger_jar} \
            --version > {output.msfragger_version}

        java \
            -Xmx64G \
            -jar {params.msfragger_jar} \
            --num_threads {threads} \
            --database_name {input.database} \
            --output_location "output/{wildcards.config_batch}/msfragger/" \
            {input.d_files}


        touch {output.untouchable}
        """







print("*/")