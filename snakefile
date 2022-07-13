    
# alias smk='mv logs/*.txt logs/old 2> /dev/null; snakemake --profile profiles/slurm'


__author__: "Carl Mathias Kobel & Arturo Vera De Ponce Leon"


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
print("                                            < MS-pipeline1 >                     ")
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


# TODO: It looks like there is a problem with annotate. Should it really be running for each sample individually? I would think that once per batch should be plentiful.


# I don't know what happened to directory in the new snakemake version
# I could simply remove the directory() calls, but I want to keep them for now - doesn't hurt.


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
#df["path"] = config_d_base + "/" + df["barcode"]






print(df)
print("//")
print()




# Define default workflow
rule all:
    input: expand(["output/{config_batch}/metadata.tsv", \
                   "output/{config_batch}/msfragger/{basename}.pepXML", \
                   "output/{config_batch}/samples/{sample}/annotate.done", \
                   "output/{config_batch}/samples/{sample}/peptideprophet-{sample}.pep.xml", \
                   "output/{config_batch}/samples/{sample}/protein.tsv", \
                   "output/{config_batch}/samples/{sample}/{sample}_quant.csv"], \
                   config_batch = config_batch, \
                   sample = df["sample"], \
                   basename = df["basename"])

#"output/{config_batch}/samples/{sample}/{sample}_quant.csv"], \ # disabled until msfragger basename-sample conversion works

#"output/{config_batch}/msfragger/philosopher_database.fas", 











# Build a database of the known amino acid sequences.
rule philosopher_database:
    input: glob.glob(config_database_glob)
    output: 
        database = "output/{config_batch}/msfragger/philosopher_database.fas",
    benchmark: "output/{config_batch}/benchmarks/philosopher_database.tsv"
    threads: 8
    params:
        philosopher = config["philosopher_executable"]
    shell: """


        >&2 echo "Catting database files ..."
        # Cat all database source files into one.
        cat {input} > output/{config_batch}/msfragger/cat_database_sources.faa


        >&2 echo "Change dir ..."
        # As philosopher can't specify output files, we need to change dir.
        mkdir -p output/{config_batch}/msfragger
        cd output/{config_batch}/msfragger

        >&2 echo "Philosopher workspace clean ..."
        {params.philosopher} workspace \
            --nocheck \
            --clean 

        >&2 echo "Philosopher workspace init ..."
        {params.philosopher} workspace \
            --nocheck \
            --init 

        >&2 echo "Removing previous .fas ..."
        rm *.fas || echo "nothing to delete" # Remove all previous databases if any.

        >&2 echo "Philosopher database ..."
        {params.philosopher} database \
            --custom cat_database_sources.faa \
            --contam 

        >&2 echo "Move output ..."
        # Manually rename the philosopher output so we can grab it later
        mv *-decoys-contam-cat_database_sources.faa.fas philosopher_database.fas

        >&2 echo "Clean up ..."
        # clean up 
        rm cat_database_sources.faa


        """


# I don't like how msfragger writes in the input file directory. So I made this rule to (soft) link the files to where I actually want the output to be.
# I find that msfragger writes some files (...calibrated.mgf and .mzBIN). I would like to keep these files together with the rest of the pipeline outputs.
# Each sample is handled by a single job.

# I have a warning that I need to take care of. I think it explains some of the problems I've had with the pipeline wanting to run over and over even no inputs have changed.
#> The flag 'directory' used in rule link_input is only valid for outputs, not inputs.
#> The flag 'directory' used in rule link_input is only valid for outputs, not inputs.

rule link_input:
    # input:
    #     d_files = directory((config_d_base + "/" + df["barcode"]).tolist()) # Instead I should probably use some kind of flag. This definition could be a param
    output:
        dir = directory("output/{config_batch}/msfragger"), 
        linked_flag = touch("output/{config_batch}/msfragger/link_input.done"), # Might technically be unnecessary, but POIRAE
        d_files = directory("output/{config_batch}/msfragger/" + df["barcode"]) # This line is necessary for msfragger to pick up on its input.
    params:
        d_files = (config_d_base + "/" + df["barcode"]).tolist() # Instead I should probably use some kind of flag. This definition could be a param
    shell:"""

        ln -s {params.d_files} {output.dir}

        """



# I don't think rule metadata is pointing out anywhere right now.
rule metadata:
    input: "output/{config_batch}/msfragger/link_input.done"
    output: "output/{config_batch}/metadata.tsv"
    params: dataframe = df.to_csv(None, index_label = "index", sep = "\t")
    shell: """

        echo '''{params.dataframe}''' > {output}
    
    """



# The reason for all the subdirectories beneath is that I want to make a dir for each sample, and they all need a "workspace"
# Wouldn't it be better to annotate the database, and then copy it to each dir? No, because the binary files in each subdir are different, so that wouldn't make sense. There is going to be _a lot_ of output files for each sample further down the line, so we want to segregate early on.
# Running this job only takes a split second any way, keeping it in its own rule makes things easy to debug. So that is it. :)
rule annotate:
    input: 
        database = "output/{config_batch}/msfragger/philosopher_database.fas"
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








# TODO: Test the implications of using shadow: "minimal"
rule msfragger:
    input:
        linked_flag = "output/{config_batch}/msfragger/link_input.done",

        database = "output/{config_batch}/msfragger/philosopher_database.fas",  
        #d_files = (config_d_base + "/" + df["barcode"]).tolist() # deprecated
        d_files = ("output/{config_batch}/msfragger/" + df["barcode"]).tolist()

    output:
        #pepXML = "output/{config_batch}/msfragger/{sample}.pepXML", 
        #pepXMLs = lambda wildcards: "output/{config_batch}/msfragger/" + df[df["sample" == wildcards.sample]]["basename"] + ".pepXML"
        #pepXMLs = expand("output/{config_batch}/msfragger/{sample}.pepXML", sample = df["sample"], config_batch = config_batch)
        pepXMLs = "output/{config_batch}/msfragger/" + df["basename"] + ".pepXML",
        stdout = "output/{config_batch}/msfragger/msfragger.out.txt"
        #touch = touch("output/{config_batch}/msfra")

    # Use shadow to get rid of the pepindex files
    #shadow: "minimal" # The setting shadow: "minimal" only symlinks the inputs to the rule. Once the rule successfully executes, the output file will be moved if necessary to the real path as indicated by output.
    # Shadow doesn't work well with tee, as tee needs access to the log directory
    threads: 8
    params:
        config_d_base = config_d_base,
        msfragger_jar = config["msfragger_jar"],
        n_samples = len(df.index)
        #basename_out =  lambda wildcards: "output/{config_batch}/msfragger/" + df[df["sample" == wildcards.sample]]["basename"] + ".pepXML"

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
            {input.d_files} \
            | tee output/{wildcards.config_batch}/msfragger/msfragger.out.txt 

   

        # The last line extracts the number of lines that corresponds to the number of samples, to extract the total number of scans.
        # Without tee, the stdout would be lost and not written in the snakemake logs.



        ls output/{wildcards.config_batch}/msfragger > output/{wildcards.config_batch}/msfragger/msfragger.done


        # makes a .pepindex and a pepXML for each sample.
        # I feel like it also creates a .mgf and .mzBIN in the source directory where the .d-dirs reside
        # Should I not move the .pepXML files? No, because I'm using the --output_location argument.
        # The tutorial mentions something about moving some .tsv files after running msfragger, but I haven't seen any.
        # These output files should in theory be mitigated by using the shadow rule?



        """




# For each sample
rule prophet_filter:
    input:
        database = "output/{config_batch}/msfragger/philosopher_database.fas",
        flag = "output/{config_batch}/samples/{sample}/annotate.done",
        pepXML = lambda wildcards: "output/" + config_batch + "/msfragger/" + df[df["sample"] == wildcards.sample]["basename"] + ".pepXML",
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

        # Be noted that the warning about a missing file is not critical: https://github.com/cmkobel/MS_pipeline1/issues/2


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



# TODO: This rule ought to output the abundances named in the samples name and not the basename? I don't really see any neat way to do that 
rule ionquant:
    input:
        irrelevant = ["output/{config_batch}/samples/{sample}/ion.tsv", \
            "output/{config_batch}/samples/{sample}/peptide.tsv", \
            "output/{config_batch}/samples/{sample}/protein.fas", \
            "output/{config_batch}/samples/{sample}/proteinprophet-{sample}.prot.xml", \
            "output/{config_batch}/samples/{sample}/protein.tsv"], 
        psm = "output/{config_batch}/samples/{sample}/psm.tsv",
        pepXML = lambda wildcards: "output/" + config_batch + "/msfragger/" + df[df["sample"] == wildcards.sample]["basename"] + ".pepXML",

    output: #touch("output/{config_batch}/samples/{sample}/ionquant.done")
        csv = "output/{config_batch}/samples/{sample}/{sample}_quant.csv"
    threads: 8
    conda: "envs/openjdk.yaml"
    params:
        ionquant_jar = config["ionquant_jar"],
        config_d_base = config_d_base, # I think this one is global, thus does not need to be params-linked.
        basename = lambda wildcards: df[df["sample"] == wildcards.sample]["basename"].values[0]

    resources:
        mem_mb = 65536



    shell: """


        >&2 echo "Ionquant ..."
        java \
            -Xmx32G \
            -jar {params.ionquant_jar} \
            --threads {threads} \
            --psm {input.psm} \
            --specdir {params.config_d_base} \
            {input.pepXML} 
            # address to msfragger pepXML file


            # TODO: Ask Arturo if it makes any sense that I'm not using the pepXML from peptideprophet, but the one directly from msfragger

            # Apparently, --specdir should point to the msfragger pepxmls. Maybe, I just need to point to the msfragger dir.
            # Or maybe I need to point directly to the file.
            # --specdir output/220315_test/msfragger/20220302_A1_Slot1-01_1_1592.pepXML 
            # Maybe the other pepxml is the culprit

            #mv output/{config_batch}/msfragger/{wildcards.sample}_quant.csv output/{config_batch}/samples/{wildcards.sample}/{wildcards.sample}_quant.csv
            mv output/{config_batch}/msfragger/{params.basename}_quant.csv output/{config_batch}/samples/{wildcards.sample}/{wildcards.sample}_quant.csv


    """

rule rmarkdown:
    input:
        metadata = "output/{config_batch}/metadata.tsv",
        psms = "output/{config_batch}/samples/{sample}/psm.tsv",
        quants = "output/{config_batch}/samples/{sample}/{sample}_quant.csv", # This simply makes it only run if rule ionquant was successful.
    output:
        "output/{config_batch}/QC.html"
    conda: "envs/r-markdown.yaml"
    shell: """
        #Rscript --what scripts/QC.Rmd {input.metadata}

        #cp scripts/QC.Rmd QC.Rmd

        Rscript -e 'library(rmarkdown); rmarkdown::render("scripts/QC.rmd", "html_document")'

        #rm rmarkdown_template.rmd
        #mv rmarkdown_template.html ../{output}

"""

        


print("*/") # This is a language specific comment close tag that helps when you export the workflow as a graph



# TODO: Go through the whole pipeline one job at a time, and make sure that all outputs are managed in the rules.
# TODO: Export the dag and put it into the readme with a bit of documentation.