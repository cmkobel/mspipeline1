#profile: "profile/slurm-sigma2-saga/"

# alias smk='mv logs/*.txt logs/old 2> /dev/null; snakemake --profile profiles/slurm'


__author__ =  "Carl Mathias Kobel & Arturo Vera Ponce De Leon"
__version__ = "v1.0.1"


# TODO: prune these imports
from datetime import datetime
import glob
import os
import pandas as pd
import re
import time 


print("/*                                                                               ") # Helps with outputting to dot.
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


# Import configuration
configfile: "config.yaml"
config_batch = config["batch"]
config_d_base = config["batch_parameters"][config_batch]["d_base"]
config_database_glob = config["batch_parameters"][config_batch]["database_glob"]
config_database_glob_read = glob.glob(config_database_glob)
config_samples = config["batch_parameters"][config_batch]["samples"]
config_temp_dir = config["temp_dir"]


# Present configuration
print(f"config_batch:         '{config_batch}'")
print(f"config_d_base:        '{config_d_base}'")
print(f"config_database_glob: '{config_database_glob}:'")
if len(config_database_glob) < 1:
    raise Exception("Raised exception: no glob targets in config_database_glob") # Not tested yet.
k = 0
for i, j in enumerate(config_database_glob_read):
    print(f"  {i}) {j}")
    if i==29:
        print(f"and {len(config_database_glob_read)-29} more..")
        break
print()


# Create a dataframe with all inputs
df = pd.DataFrame(data = {'sample':  config_samples.keys(),
                          'barcode': config_samples.values()})

df["basename"] = [re.sub(".d$", "", barcode) for barcode in df["barcode"]]
print(df)
print("//")
print()


# Define workflow targets
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



# Save some metadata about in puts for good measure.
rule metadata:
    #input: "output/{config_batch}/msfragger/link_input.done"
    output: "output/{config_batch}/metadata.tsv"
    params: dataframe = df.to_csv(None, index_label = "index", sep = "\t")
    shell: """

        echo '''{params.dataframe}''' > {output}
    
    """




# Create a symbolic link for the input files. Msfragger writes adjacent to the input directories, so linking keeps these outputs somewhat isolated. 
# I find that msfragger writes some files (...calibrated.mgf and .mzBIN). I would like to keep these files together with the rest of the pipeline outputs.
# Update: When using slow long-term backup storage it might be better to copy the files to the wd drive instead of linking.
rule link_input:
    output:
        dir = directory("output/{config_batch}/msfragger"), 
        d_files = directory("output/{config_batch}/msfragger/" + df["barcode"] + "/"), # Bound for msfragger.
        linked_flag = touch("output/{config_batch}/msfragger/link_input.done") # Used by rule philosopher_database to wait for creation of the msfragger directory.
        # Make sure you've set write access to the directory where these files reside.
    params:
        d_files = (config_d_base + "/" + df["barcode"]).tolist() # Instead I should probably use some kind of flag. This definition could be a param.
    shell: """
        
        ln -s {params.d_files} {output.dir}
        #cp -r {params.d_files} {output.dir}
        # I'd rather manually copy the files and the link them. Otherwise snakemake will make new copies all the freakin' time.


    """



# Build a database of the known amino acid sequences.
rule philosopher_database:
    input: 
        glob = glob.glob(config_database_glob),
        linked_flag = "output/{config_batch}/msfragger/link_input.done" 
    output: 
        database = "output/{config_batch}/msfragger/philosopher_database.fas",
    benchmark: "output/{config_batch}/benchmarks/philosopher_database.tsv"
    threads: 8
    #retries: 3 # Disabled: Increasing memory doesn't burn through to slurm anyway.
    resources:
        mem_mb = lambda wildcards, attempt: 32768 * (2**attempt//2), # multiply by 1, 2, 4, 8 # This is not yet tested.  # Shows up in the log snakemake stdout but doesn't burn through to slurm.      
    params:
        philosopher = config["philosopher_executable"]
    shell: """


        #TMPDIR={config_temp_dir} # tmp_dir should be set by snakemake


        >&2 echo "Catting database files ..."
        # Cat all database source files into one.
        cat {input.glob} > output/{config_batch}/msfragger/cat_database_sources.faa


        >&2 echo "Change dir ..."
        # As philosopher can't be specified output files, we need to change dir.
        cd output/{config_batch}/msfragger

        >&2 echo "Philosopher workspace clean ..."
        {params.philosopher} workspace --nocheck --clean 

        >&2 echo "Philosopher workspace init ..."
        {params.philosopher} workspace --nocheck --init --temp $TMPDIR
        


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




# Create a philosopher "workspare" for each sample in dedicated directories.
# Creating a single workspace and copying it out is _not_ a solution as one of the binary files (in the hidden workspace sub-directory) defines the location. So we need to recalculate the same thing many times over.
rule annotate:
    input: 
        database = "output/{config_batch}/msfragger/philosopher_database.fas"
    output: 
        flag = touch("output/{config_batch}/samples/{sample}/annotate.done")
    params:
        philosopher = config["philosopher_executable"]
    resources:
        mem_mb = 65536
    shell: """

        #TMPDIR={config_temp_dir} # tmp_dir should be set by snakemake


        mkdir -p output/{config_batch}/samples/{wildcards.sample}
        cd output/{config_batch}/samples/{wildcards.sample}

        {params.philosopher} workspace --nocheck --clean

        {params.philosopher} workspace --nocheck --init --temp $TMPDIR



        >&2 echo "Annotating database ..."
        {params.philosopher} database \
            --annotate ../../../../{input.database}        


    """



# Match the PSMs to the database
# I've considered using shadow to prune some of the unneeded outputs, but since I don't know exactly which outputs I'm going to need later on, I think it is too much of a rabbit hole to dive into right now.
rule msfragger:
    input:
        # linked_flag = "output/{config_batch}/msfragger/link_input.done", # Not necessary since d_files are well defined.
        database = "output/{config_batch}/msfragger/philosopher_database.fas",  
        d_files = ("output/{config_batch}/msfragger/" + df["barcode"]).tolist()
    output:
        pepXMLs = "output/{config_batch}/msfragger/" + df["basename"] + ".pepXML",
        stdout = "output/{config_batch}/msfragger/msfragger.out.txt"
    # Use shadow to get rid of the pepindex files
    # shadow: "minimal" # The setting shadow: "minimal" only symlinks the inputs to the rule. Once the rule successfully executes, the output file will be moved if necessary to the real path as indicated by output.
    # Shadow doesn't work well with tee, as tee needs access to the log directory. Too much complexity.
    threads: 8
    params:
        config_d_base = config_d_base,
        msfragger_jar = config["msfragger_jar"],
        n_samples = len(df.index), 
    resources:
        mem_mb = 515538, # will be overwritten by set-resources in the profile, so remove that before managing it here.
        #mem_mb = 120000,
        partition = 'bigmem',
        runtime = '5-00:00:00'
        
    conda: "envs/openjdk.yaml"
    shell: """

        >&2 echo "MSFragger ..."
        java \
            -Xmx500G \
            -jar {params.msfragger_jar} \
            --num_threads {threads} \
            --database_name {input.database} \
            --output_location "output/{wildcards.config_batch}/msfragger/" \
            {input.d_files} \
            | tee output/{wildcards.config_batch}/msfragger/msfragger.out.txt 


        # We need to collect stdout so we can later compare the number of PSMs to the number of scans.

        # Get an overview of which files were created by msfragger. 
        ls -l output/{wildcards.config_batch}/msfragger > output/{wildcards.config_batch}/msfragger/msfragger.done


        # makes a .pepindex and a pepXML for each sample.
        # I feel like it also creates a .mgf and .mzBIN in the source directory where the .d-dirs reside
        # Should I not move the .pepXML files? No, because I'm using the --output_location argument.
        # The tutorial mentions something about moving some .tsv files after running msfragger, but I haven't seen any.
        # These output files should in theory be mitigated by using the shadow rule? Update: No.


        """





# Filter the raw msfragger output.
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
    resources:
        mem_mb = 64000
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
# I assume that ionquant changes the files in place?
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
        csv = "output/{config_batch}/samples/{sample}/{sample}_quant.csv",
        tsv = "output/{config_batch}/samples/{sample}/{sample}_protein.tsv"
        # Why not one of the many other files? quant is produced by msfragger you know.
    threads: 8
    conda: "envs/openjdk.yaml"
    params:
        ionquant_jar = config["ionquant_jar"],
        config_d_base = config_d_base, # I think this one is global, thus does not need to be params-linked.
        basename = lambda wildcards: df[df["sample"] == wildcards.sample]["basename"].values[0] # used to have a [0] in the end, which I just removed, and it started working again.
    resources:
        mem_mb = 65536
        #mem_mb = lambda wildcards, attempt: 16384 * (2**attempt//2) # multiply by 1, 2, 4, 8 # This is not yet tested.
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


        # mv output/{config_batch}/msfragger/{params.basename}_quant.csv output/{config_batch}/samples/{wildcards.sample}/{wildcards.sample}_quant.csv

        # Instead of simply moving that file, I might want to prepend it with its sample name:
        # This one I think is needed for rate calculation, but why don't I do it inside the msfragger job?
        # TODO move it to msfragger
        

        cp output/{config_batch}/msfragger/{params.basename}_quant.csv output/{config_batch}/samples/{wildcards.sample}/{wildcards.sample}_quant.csv
        

        # Turns out that it is really the proteins file that we're interested in
        cp output/{config_batch}/samples/{wildcards.sample}/protein.tsv output/{config_batch}/samples/{wildcards.sample}/{wildcards.sample}_protein.tsv



    """









# This is not yet implemented.
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

        


print("*/") # This is a dot-language specific comment close tag that helps when you export the workflow as a graph



# TODO: Go through the whole pipeline one job at a time, and make sure that all outputs are managed in the rules. Update: seems legit bruh.
