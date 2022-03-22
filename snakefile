    
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
print("                                             ______________     ")
print("                                            < MS_pipeline1 >    ")
print("                                             --------------     ")
print("                                                          \\      ")
print("                             ___......__             _     \\                     ")
print("                         _.-'           ~-_       _.=a~~-_                       ")
print(" --=====-.-.-_----------~   .--.       _   -.__.-~ ( ___===>                     ")
print("               '''--...__  (    \\ \\\\\\ { )       _.-~                         ")
print("                         =_ ~_  \\\\-~~~//~~~~-=-~                               ")
print("                          |-=-~_ \\\\   \\\\                                     ")
print("                          |_/   =. )   ~}                                        ")
print("                          |}      ||                                             ")
print("                         //       ||                                             ")
print("                       _//        {{                                             ")
print("                    '='~'          \\\\_    =                                    ")
print("                                    ~~'                                          ")
print("                                                                                 ")





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

df["basename"] = [re.sub(".d$", "", barcode) for barcode in df["barcode"]]
df["path"] = config_d_base + "/" + df["barcode"]


#df["pepXML"] = "output/" + config_batch + "/incremental_results/" + df["basename"]


print(df)
print("//")

#print(df["path"].tolist()) # Debug



# TODO: icremental_results should be renamed to "batch_results" or "prepare" or "prepare_batch"
# TODO: fix the problem that means that I can't have snakemake running batch first and then per-sample. Maybe I need to use lambda?
#       Well, asscom2 really is the other way around. It first takes samples, and then goes into samples. 

# Define default workflow
rule all:
    input: expand(["output/{config_batch}/incremental_results/philosopher_database.fas", \
                   "output/{config_batch}/incremental_results/{basename}.pepXML", \
                   "output/{config_batch}/samples/{sample}/annotate.done", \
                   "output/{config_batch}/samples/{sample}/peptideprophet-{sample}.pep.xml"], \
                   config_batch = config_batch, \
                   sample = df["sample"], \
                   basename = df["basename"])

rule philosopher_database:
    input: glob.glob(config_database_glob)
    output: 
        database = "output/{config_batch}/incremental_results/philosopher_database.fas",
        flag = touch("output/{config_batch}/incremental_results/philosopher_database.done")
    #benchmark: "output/{config_batch}/benchmarks/database.tab"
    threads: 8
    params:
        philosopher = config["philosopher_executable"]
    shell: """



        # Cat all database source files into one.
        cat {input} > output/{config_batch}/incremental_results/cat_database_sources.faa


        # As philosopher can't specify output files, we need to change dir.
        mkdir -p output/{config_batch}/incremental_results
        cd output/{config_batch}/incremental_results


        {params.philosopher} workspace \
            --nocheck \
            --clean 

        {params.philosopher} workspace \
            --nocheck \
            --init 

        rm *.fas || echo "nothing to delete" # Remove all previous databases if any.
        {params.philosopher} database \
            --custom cat_database_sources.faa \
            --contam 

        # Manually rename the philosopher output so we can grab it later
        mv *.fas philosopher_database.fas


        ls * > philosopher_database.done


        """



rule msfragger:
    input:
        database = "output/{config_batch}/incremental_results/philosopher_database.fas",  
        d_files = df["path"].tolist()
    output: 
        #msfragger_version = "output/{config_batch}/incremental_results/msfragger_version.txt",
        #untouchable = touch("output/{config_batch}/incremental_results/msfragger.done"),
        #pepXMLs = expand("output/{config_batch}/incremental_results/{pepXMLs}", config_batch = config_batch, pepXMLs = df["pepXML"])
        pepXMLs = "output/{config_batch}/incremental_results/{basename}.pepXML"
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
            --num_threads {threads} \
            --database_name {input.database} \
            --output_location "output/{wildcards.config_batch}/incremental_results/" \
            {input.d_files}


        ls output/{wildcards.config_batch}/incremental_results > output/{wildcards.config_batch}/incremental_results/msfragger.done


        # makes a .pepindex and a pepXML for each sample.

        """



rule annotate:
    input: 
        database = "output/{config_batch}/incremental_results/philosopher_database.fas"
    output: 
        flag = touch("output/{config_batch}/samples/{sample}/annotate.done")
    params:
        philosopher = config["philosopher_executable"]
    shell: """
        mkdir -p output/{config_batch}/samples/{wildcards.sample}
        cd output/{config_batch}/samples/{wildcards.sample}

        {params.philosopher} workspace \
            --nocheck \
            --clean

        {params.philosopher} workspace \
            --nocheck \
            --init

        echo "Annotating database ..."
        {params.philosopher} database \
            --annotate ../../../../{input.database}

        
        

    """



# For each sample
rule peptideprophet:
    input:
        database = "output/{config_batch}/incremental_results/philosopher_database.fas",
        flag = "output/{config_batch}/samples/{sample}/annotate.done",
        #pepXML = "output/{config_batch}/incremental_results/{basename}.pepXML"
        pepXML = lambda wildcards: "output/{config_batch}/incremental_results/" + df[df["sample"]==wildcards.sample]["basename"].values[0] + ".pepXML"
    output:
        peptide = "output/{config_batch}/samples/{sample}/peptideprophet-{sample}.pep.xml"
        #protein = "output/{config_batch}/samples/{sample}/proteinprophet-{sample}.prot.xml"
    params:
        philosopher = config["philosopher_executable"]
    shell: """

        # Since the output location of philosopher is controlled by the input location, we should copy the input file.
        cp {input.pepXML} output/{wildcards.config_batch}/samples/{wildcards.sample}/{wildcards.sample}.pepXML || echo "file exists already"

        # Because of the workspace, we're forced to change dir
        cd output/{wildcards.config_batch}/samples/{wildcards.sample}


        
        echo "Peptideprophet ..."
        {params.philosopher} peptideprophet \
            --nonparam \
            --expectscore \
            --decoyprobs \
            --ppm \
            --accmass \
            --output peptideprophet \
            --database ../../../../{input.database} \
            {wildcards.sample}.pepXML
            #../../../../{input.pepXML}



        echo "Proteinprophet ..."
        {params.philosopher} proteinprophet \
            --output proteinprophet-{wildcards.sample} \
            peptideprophet-{wildcards.sample}.pep.xml


        echo "Filter ..." 
        {params.philosopher} filter \
            --sequential \
            --razor \
            --mapmods \
            --pepxml peptideprophet-{wildcards.sample}.pep.xml \
            --protxml proteinprophet-{wildcards.sample}.prot.xml

        # Assuming that philosopher filter works in place

        {params.philosopher} report

    
    """












print("*/")