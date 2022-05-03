    
# alias smk='mv logs/*.txt logs/old 2> /dev/null; snakemake --profile profiles/slurm'


__author__: "Carl Mathias Kobel & Arturo Vera De Ponce Leon"



# who wrote this


import os
from datetime import datetime
import time
import re
#from shutil import copyfile
import glob

import pandas as pd
import re


print("/*                                                                               ")
print("                                             ______________                      ")
print("                                            < MS_pipeline1 >                     ")
print("                                             --------------                      ")
print("                                                          \\                     ")
print("                             ___......__             _     \\                    ")
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






print(df)
print("//")

#print(df["path"].tolist()) # Debug





# TODO: Ask Arturo about the weird bug that is also mentioned here https://github.com/Nesvilab/FragPipe/issues/4#issuecomment-309045306

# Define default workflow
rule all:
    input: expand(["output/{config_batch}/msfragger/philosopher_database.fas", \
                   "output/{config_batch}/msfragger/{basename}.pepXML", \
                   "output/{config_batch}/samples/{basename}/annotate.done", \
                   "output/{config_batch}/samples/{basename}/peptideprophet-{basename}.pep.xml", \
                   "output/{config_batch}/samples/{basename}/{basename}_quant.csv"], \
                   config_batch = config_batch, \
                   sample = df["sample"], \
                   basename = df["basename"])



# Build a database of the known amino acid sequences.
rule philosopher_database:
    input: glob.glob(config_database_glob)
    output: 
        database = "output/{config_batch}/msfragger/philosopher_database.fas",
        flag = touch("output/{config_batch}/msfragger/philosopher_database.done") # This flag is redundant but nice to have.
    #benchmark: "output/{config_batch}/benchmarks/database.tab"
    threads: 8
    params:
        philosopher = config["philosopher_executable"]
    shell: """

        # Show help for debugging
        {params.philosopher} database \
            --help 

        # Cat all database source files into one.
        cat {input} > output/{config_batch}/msfragger/cat_database_sources.faa


        # As philosopher can't specify output files, we need to change dir.
        mkdir -p output/{config_batch}/msfragger
        cd output/{config_batch}/msfragger


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
        mv *-decoys-contam-cat_database_sources.faa.fas philosopher_database.fas

        # clean up 
        rm cat_database_sources.faa


        ls -lt * > philosopher_database.done


        """





# I'm considering if msfragger should run individually for each sample: It seems to be fundamentally impossible. 
# TODO: Run a test without shadow: minimal.
rule msfragger:
    input:
        database = "output/{config_batch}/msfragger/philosopher_database.fas",  
        d_files = df["path"].tolist()
    output:
        pepXMLs = "output/{config_batch}/msfragger/{basename}.pepXML"
    shadow: "minimal" # The setting shadow: "minimal" only symlinks the inputs to the rule. Once the rule successfully executes, the output file will be moved if necessary to the real path as indicated by output.
    threads: 8
    params:
        config_d_base = config_d_base,
        msfragger_jar = config["msfragger_jar"]
    conda: "envs/openjdk.yaml"
    shell: """

        # We can't manage the output path of MSFragger, so we need to symlink it to the directory where we want the output to reside.
        # And no, we can't use snakemake-shadow, as it only works on relative subdirs.

        >&2 echo "MSFragger ..."
        java \
            -Xmx64G \
            -jar {params.msfragger_jar} \
            --num_threads {threads} \
            --database_name {input.database} \
            --output_location "output/{wildcards.config_batch}/msfragger/" \
            {input.d_files}



        ls output/{wildcards.config_batch}/msfragger > output/{wildcards.config_batch}/msfragger/msfragger.done


        # makes a .pepindex and a pepXML for each sample.
        # I feel like it also creates a .mgf and .mzBIN in the source directory where the .d-dirs reside
        # Should I not move the .pepXML files? No, because I'm using the --output_location argument.
        # The tutorial mentions something about moving some .tsv files after running msfragger, but I haven't seen any.
        # TODO: run crystal-c


        """


# crystalc has been disabled because the documentation is incomplete.
#rule crystalc:
#    input: "output/{config_batch}/msfragger/{basename}.pepXML"
#    output: "output/{config_batch}/crystalc/{basename}.pepXML"
#    params: 
#        crystalc_jar = config["crystalc_jar"]
#    conda: "envs/openjdk.yaml"
#    shell: """
#
#        >&2 echo "Shifting crystalc ..."
#    
#        java -Xmx64G \
#            -jar {params.crystalc_jar} \
#            --output_location output/{config_batch}/crystalc \
#            assets/crystalc.params \
#            {input}
#
#        #Since the documentiation is so rubbish for crystalc, I will not use it for now and 
#
#        # I would like to see if I can abstain from using the crystalcParameterPath at all
#
#    """




# The reason for all the subdirectories beneath is that I want to make a dir for each sample, and they all need a "workspace"
# Wouldn't it be better to annotate the database, and then copy it to each dir? No, because the binary files in each subdir are different, so that wouldn't make sense. There is going to be _a lot_ of output files for each sample further down the line, so we want to segregate early on.
# Running this job only takes a split second any way, keeping it in its own rule makes things easy to debug. So that is it.
rule annotate:
    input: 
        database = "output/{config_batch}/msfragger/philosopher_database.fas"
    output: 
        flag = touch("output/{config_batch}/samples/{basename}/annotate.done")
    params:
        philosopher = config["philosopher_executable"]
    shell: """
        mkdir -p output/{config_batch}/samples/{wildcards.basename}
        cd output/{config_batch}/samples/{wildcards.basename}

        {params.philosopher} workspace \
            --nocheck \
            --clean

        {params.philosopher} workspace \
            --nocheck \
            --init

        >&2 echo "Annotating database ..."
        {params.philosopher} database \
            --annotate ../../../../{input.database}

        
        

    """



# For each sample
rule peptideprophet:
    input:
        database = "output/{config_batch}/msfragger/philosopher_database.fas",
        flag = "output/{config_batch}/samples/{basename}/annotate.done",
        #pepXML = lambda wildcards: "output/{config_batch}/msfragger/" + df[df["sample"]==wildcards.sample]["basename"].values[0] + ".pepXML"
        pepXML = "output/{config_batch}/msfragger/{basename}.pepXML"
    output: ["output/{config_batch}/samples/{basename}/ion.tsv", \
        "output/{config_batch}/samples/{basename}/peptideprophet-{basename}.pep.xml", \
        "output/{config_batch}/samples/{basename}/peptide.tsv", \
        "output/{config_batch}/samples/{basename}/protein.fas", \
        "output/{config_batch}/samples/{basename}/proteinprophet-{basename}.prot.xml", \
        "output/{config_batch}/samples/{basename}/protein.tsv", \
        "output/{config_batch}/samples/{basename}/psm.tsv"]
        #protein = "output/{config_batch}/samples/{basename}/proteinprophet-{basename}.prot.xml"
    params:
        philosopher = config["philosopher_executable"]
    shell: """

        # Since the output location of philosopher is controlled by the input location, we should copy the input file.
        cp {input.pepXML} output/{wildcards.config_batch}/samples/{wildcards.basename}/{wildcards.basename}.pepXML || echo "file exists already"

        # Because of the workspace, we're forced to change dir
        cd output/{wildcards.config_batch}/samples/{wildcards.basename}


        
        >&2 echo "Peptideprophet ..."
        {params.philosopher} peptideprophet \
            --nonparam \
            --expectscore \
            --decoyprobs \
            --ppm \
            --accmass \
            --output peptideprophet \
            --database ../../../../{input.database} \
            {wildcards.basename}.pepXML




        >&2 echo "Proteinprophet ..."
        {params.philosopher} proteinprophet \
            --output proteinprophet-{wildcards.basename} \
            peptideprophet-{wildcards.basename}.pep.xml


        >&2 echo "Filter ..." 
        {params.philosopher} filter \
            --sequential \
            --razor \
            --mapmods \
            --pepxml peptideprophet-{wildcards.basename}.pep.xml \
            --protxml proteinprophet-{wildcards.basename}.prot.xml

        # Assuming that philosopher filter works in place
        # TODO: Ask Arturo if that is true.
        >&2 echo "Report ..."
        {params.philosopher} report



        
    """

# TODO: This rule ought to output the abundances named in the samples name and not the basename?
rule ionquant:
    input:
        irrelevant = ["output/{config_batch}/samples/{basename}/ion.tsv", \
            "output/{config_batch}/samples/{basename}/peptide.tsv", \
            "output/{config_batch}/samples/{basename}/protein.fas", \
            "output/{config_batch}/samples/{basename}/proteinprophet-{basename}.prot.xml", \
            "output/{config_batch}/samples/{basename}/protein.tsv"], 
        psm = "output/{config_batch}/samples/{basename}/psm.tsv",
        #pepxml = "output/{config_batch}/samples/{basename}/peptideprophet-{basename}.pep.xml",
        pepxml = "output/{config_batch}/msfragger/{basename}.pepXML"
    output: #touch("output/{config_batch}/samples/{basename}/ionquant.done")
        csv = "output/{config_batch}/samples/{basename}/{basename}_quant.csv"
    threads: 8
    conda: "envs/openjdk.yaml"
    params:
        ionquant_jar = config["ionquant_jar"],
        #spectral = lambda wildcards: df[df["sample"]==wildcards.sample]["path"].values[0]
        #spectral = lambda wildcards: df[df["basename"]==wildcards.basename]["path"].values[0]
        config_d_base = config_d_base


    shell: """


        >&2 echo "Ionquant ..."
        java \
            -Xmx32G \
            -jar {params.ionquant_jar} \
            --threads {threads} \
            --psm {input.psm} \
            --specdir {params.config_d_base} \
            {input.pepxml} 
            # address to msfragger pepXML file


            # TODO: Ask Arturo if it makes any sense that I'm not using the pepXML from peptideprophet, but the one directly from msfragger

            # Apparently, --specdir should point to the msfragger pepxmls. Maybe, I just need to point to the msfragger dir.
            # Or maybe I need to point directly to the file.
            # --specdir output/220315_test/msfragger/20220302_A1_Slot1-01_1_1592.pepXML 
            # Maybe the other pepxml is the culprit

            mv output/{config_batch}/msfragger/{wildcards.basename}_quant.csv output/{config_batch}/samples/{wildcards.basename}/{wildcards.basename}_quant.csv




    """
        





print(df["path"].tolist())


print("*/")


# TODO: 
#   a) get ionquant working
#   b) put msfragger into segregated jobs for each sample. (Uses a bit more resources, but is much faster!)

