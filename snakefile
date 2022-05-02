    
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


#df["pepXML"] = "output/" + config_batch + "/philosopher_workspace/" + df["basename"]


print(df)
print("//")

#print(df["path"].tolist()) # Debug



# TODO: icremental_results should be renamed to "batch_results" or "prepare" or "prepare_batch"
# TODO: fix the problem that means that I can't have snakemake running batch first and then per-sample. Maybe I need to use lambda?
#       Well, asscom2 really is the other way around. It first takes samples, and then goes into samples. 
# TODO: Ask Arturo about the weird bug that is also mentioned here https://github.com/Nesvilab/FragPipe/issues/4#issuecomment-309045306

# Define default workflow
rule all:
    input: expand(["output/{config_batch}/philosopher_workspace/philosopher_database.fas", \
                   "output/{config_batch}/msfragger/{basename}.pepXML", \
                   "output/{config_batch}/crystalc/{basename}.pepXML", \
                   "output/{config_batch}/samples/{sample}/annotate.done", \
                   "output/{config_batch}/samples/{sample}/peptideprophet-{sample}.pep.xml", \
                   "output/{config_batch}/samples/{sample}/ionquant.done"], \
                   config_batch = config_batch, \
                   sample = df["sample"], \
                   basename = df["basename"])

# TODO: Rename incremental results to philosopher_database
# TODO: Move this rule down under msfragger and crystalc so it starts looking more like the tutorial I bookmarked.
rule philosopher_database:
    input: glob.glob(config_database_glob)
    output: 
        database = "output/{config_batch}/philosopher_workspace/philosopher_database.fas",
        flag = touch("output/{config_batch}/philosopher_workspace/philosopher_database.done") # This flag is redundant but nice to have.
    #benchmark: "output/{config_batch}/benchmarks/database.tab"
    threads: 8
    params:
        philosopher = config["philosopher_executable"]
    shell: """



        # Cat all database source files into one.
        cat {input} > output/{config_batch}/philosopher_workspace/cat_database_sources.faa


        # As philosopher can't specify output files, we need to change dir.
        mkdir -p output/{config_batch}/philosopher_workspace
        cd output/{config_batch}/philosopher_workspace


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
        database = "output/{config_batch}/philosopher_workspace/philosopher_database.fas",  
        d_files = df["path"].tolist()
    output: 
        #msfragger_version = "output/{config_batch}/philosopher_workspace/msfragger_version.txt",
        #untouchable = touch("output/{config_batch}/philosopher_workspace/msfragger.done"),
        #pepXMLs = expand("output/{config_batch}/philosopher_workspace/{pepXMLs}", config_batch = config_batch, pepXMLs = df["pepXML"])
        pepXMLs = "output/{config_batch}/msfragger/{basename}.pepXML"
    #shadow: "shallow" # Can't use shadow as it isn't a subdir.
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


#print(df["basename"].tolist())
print(df["basename"].values)
# Nesvilabs software er en dårlig kombination med snakemake fordi der ikke kommer nogen output filerw
# Hvordan får jeg snakemake til at expandere ud over alle samples. Jeg vil gerne køre den for hver. Måske står det længere nede i min kode.
rule crystalc:
    input: 
        lambda wildcards: "output/{config_batch}/msfragger/" + df[df["sample"]==wildcards.sample]["basename"].values + ".pepXML"
    output: 
        #lambda wildcards: "output/{config_batch}/crystalc/" + df[df["sample"]==wildcards.sample]["basename"].values[0] + ".pepXML"
        "output/{config_batch}/crystalc/{sample}.pepXML"

    # Det her skal selvfølgelig testes og alt det der. Men til det behøver jeg nok ikke at genkøre hele msfragger. Jeg kan bare gemme output-mappen som noget bestemt, fx. msoutput, og så toggle frem og tilbage indtil koden virker pålidelig. 
    params: 
        crystalc_jar = config["crystalc_jar"]

    shell: """
    
        java -Xmx64G \
            -jar {params.crystalc_jar} \
            --output_location "output/{wildcards.config_batch}/crystalc" \
            {input}

        # I would like to see if I can abstain from using the crystalcParameterPath at all

    """




rule annotate:
    input: 
        database = "output/{config_batch}/philosopher_workspace/philosopher_database.fas"
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

        >&2 echo "Annotating database ..."
        {params.philosopher} database \
            --annotate ../../../../{input.database}

        
        

    """



# For each sample
rule peptideprophet:
    input:
        database = "output/{config_batch}/philosopher_workspace/philosopher_database.fas",
        flag = "output/{config_batch}/samples/{sample}/annotate.done",
        #pepXML = "output/{config_batch}/philosopher_workspace/{basename}.pepXML"
        pepXML = lambda wildcards: "output/{config_batch}/crystalc/" + df[df["sample"]==wildcards.sample]["basename"].values[0] + ".pepXML"
    output: ["output/{config_batch}/samples/{sample}/ion.tsv", \
        "output/{config_batch}/samples/{sample}/peptideprophet-{sample}.pep.xml", \
        "output/{config_batch}/samples/{sample}/peptide.tsv", \
        "output/{config_batch}/samples/{sample}/protein.fas", \
        "output/{config_batch}/samples/{sample}/proteinprophet-{sample}.prot.xml", \
        "output/{config_batch}/samples/{sample}/protein.tsv", \
        "output/{config_batch}/samples/{sample}/psm.tsv"]
        #protein = "output/{config_batch}/samples/{sample}/proteinprophet-{sample}.prot.xml"
    params:
        philosopher = config["philosopher_executable"]
    shell: """

        # Since the output location of philosopher is controlled by the input location, we should copy the input file.
        cp {input.pepXML} output/{wildcards.config_batch}/samples/{wildcards.sample}/{wildcards.sample}.pepXML || echo "file exists already"

        # Because of the workspace, we're forced to change dir
        cd output/{wildcards.config_batch}/samples/{wildcards.sample}


        
        >&2 echo "Peptideprophet ..."
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



        >&2 echo "Proteinprophet ..."
        {params.philosopher} proteinprophet \
            --output proteinprophet-{wildcards.sample} \
            peptideprophet-{wildcards.sample}.pep.xml


        >&2 echo "Filter ..." 
        {params.philosopher} filter \
            --sequential \
            --razor \
            --mapmods \
            --pepxml peptideprophet-{wildcards.sample}.pep.xml \
            --protxml proteinprophet-{wildcards.sample}.prot.xml

        # Assuming that philosopher filter works in place
        # TODO: Ask Arturo if that is true.
        >&2 echo "Report ..."
        {params.philosopher} report

        
    """


rule ionquant:
    input:
        irrelevant = ["output/{config_batch}/samples/{sample}/ion.tsv", \
            "output/{config_batch}/samples/{sample}/peptide.tsv", \
            "output/{config_batch}/samples/{sample}/protein.fas", \
            "output/{config_batch}/samples/{sample}/proteinprophet-{sample}.prot.xml", \
            "output/{config_batch}/samples/{sample}/protein.tsv"], 
        psm = "output/{config_batch}/samples/{sample}/psm.tsv",
        pepxml = "output/{config_batch}/samples/{sample}/peptideprophet-{sample}.pep.xml"
    output: touch("output/{config_batch}/samples/{sample}/ionquant.done")
    threads: 8
    conda: "envs/openjdk.yaml"
    params:
        ionquant_jar = config["ionquant_jar"],
        spectral = lambda wildcards: df[df["sample"]==wildcards.sample]["path"].values[0]

    shell: """


        >&2 echo "Ionquant ..."
        java \
            -Xmx32G \
            -jar {params.ionquant_jar} \
            --threads {threads} \
            --specdir {params.spectral} \
            --psm {input.psm} \
            {input.pepxml}

    """
        


#print(df["path"].tolist())


print("*/")
