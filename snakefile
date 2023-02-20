# snakemake --profile profiles/slurm-sigma2-saga
# I'm experiencing some major problems with the temporary directories that i might as well fix. It seems to percolate through when I have a high amount of samples. Basically, all the jobs that use a program that uses the workspace, need to be in the same rule. Silly, but that is how it is.

# In this branch I'm not at all screwing around. I'm closely following this tutorial:
# https://fragpipe.nesvilab.org/docs/tutorial_linux.html


__author__ =  "Carl Mathias Kobel & Arturo Vera Ponce De Leon"
__version__ = "v1.0.1"


# TODO: prune these imports
from datetime import datetime
import glob
import os
import pandas as pd
import re
import time 
import random


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


# Read configuration
configfile: "config.yaml"
config_batch = config["batch"]
config_d_base = config["batch_parameters"][config_batch]["d_base"]
config_database_glob = config["batch_parameters"][config_batch]["database_glob"]
config_database_glob_read = glob.glob(config_database_glob)
config_samples = config["batch_parameters"][config_batch]["samples"]



# Present configuration
print(f"        config_batch: '{config_batch}'")
print(f"       config_d_base: '{config_d_base}'")
print(f"config_database_glob: '{config_database_glob}:'")
if len(config_database_glob) < 1:
    raise Exception("Raised exception: no glob targets in config_database_glob") # Not tested yet.

for i, j in enumerate(config_database_glob_read):
    print(f"  {i}) {j}")
    if i==19: # Only show up till 30 lines, otherwise the screen will become congested. 
        print(f"and {len(config_database_glob_read)-19} more.. ({len(config_database_glob)} in total)")
        break # Do not print more.
print()



# Create a dataframe with all inputs
df = pd.DataFrame(data = {'sample':  config_samples.keys(),
                          'barcode': config_samples.values()})

df["basename"] = [re.sub(".d$", "", barcode) for barcode in df["barcode"]]
print(df)
print("//")
print()


# Count input sizes for managing rule resources
n_genomes_database = len(config_database_glob_read)
n_samples = len(df["basename"])
print(f"n_genomes_database: {n_genomes_database}")
print(f"n_samples: {n_samples}")    


print("manifest:")
manifest = pd.DataFrame(data = {'path': config_samples.values()})
manifest['experiment'] = "" # Experiment (can be empty, alphanumeric, and _)
manifest['bioreplicate'] = "" # Bioreplicate (can be empty and integer)
manifest['data_type'] = "DDA" # Data type (DDA, DIA, GPF-DIA, DIA-Quant, DIA-Lib)
print(manifest)
print("//")


# Mock manifest file below, not sure about the many tabs between 1 and DDA. 
# /localscratch/7407431/Fragpipe.tryar_A.dir/20220506_A1_Slot1-01_1_1939.d   testa   1			DDA
# /localscratch/7407431/Fragpipe.tryar_A.dir/20220506_A2_Slot1-02_1_1940.d   testa   1			DDA
# /localscratch/7407431/Fragpipe.tryar_A.dir/20220506_A3_Slot1-03_1_1941.d   testa   1			DDA
# /localscratch/7407431/Fragpipe.tryar_A.dir/20220506_A4_Slot1-04_1_1942.d   testb   1			DDA
# /localscratch/7407431/Fragpipe.tryar_A.dir/20220506_A5_Slot1-05_1_1943.d   testb   1			DDA
# /localscratch/7407431/Fragpipe.tryar_A.dir/20220506_A6_Slot1-06_1_1944.d   testb   1			DDA




# Define workflow targets
rule all:
    input: expand(["output/{config_batch}/metadata.tsv", \
                   "output/{config_batch}/philosopher_database.fas", \
                   "output/{config_batch}/msfragger/link_input.done", \
                   "output/{config_batch}/final.flag"], \
                   config_batch = config_batch, \
                   sample = df["sample"], \
                   basename = df["basename"])



# Save some metadata about inputs for good measure.
rule metadata:
    output: "output/{config_batch}/metadata.tsv"
    params: dataframe = df.to_csv(None, index_label = "index", sep = "\t")
    shell: """

        echo '''{params.dataframe}''' > {output}
    
    """

# Link input links or copies the input data to a specific directory. Long term, this should be on the fastest possible disk ie. userwork.
rule link_input:
    output:
        dir = directory("output/{config_batch}/msfragger"), 
        d_files = directory("output/{config_batch}/msfragger/" + df["barcode"]), # Bound for msfragger.
        linked_flag = touch("output/{config_batch}/msfragger/link_input.done"), # Used by rule philosopher_database to wait for creation of the msfragger directory.
        # Make sure you've set write access to the directory where these files reside.
    params:
        d_files = (config_d_base + "/" + df["barcode"]).tolist(), # Instead I should probably use some kind of flag. This definition could be a param.
    benchmark: "output/{config_batch}/benchmarks/benchmark.link_input.{config_batch}.tsv"
    shell: """
        
        #ln -s {params.d_files} {output.dir}
        cp -r {params.d_files} {output.dir}
        # I'd rather manually copy the files and then link them with this rule. Otherwise snakemake will make new copies all the freakin' time. When the pipeline becomes stable I can change it to linking.

        # Arturo showed me a cool trick of using rsync to copy the files, because then it shows the write speed.

    """


# Make database cats all the amino acid fastas together and runs philosopher database on it.
rule make_database:
    input:
        glob = [glob.glob(config_database_glob)],
    output:
        database = "output/{config_batch}/philosopher_database.fas",
    params:
        philosopher = config["philosopher_executable"],
    retries: 4
    resources:
        mem_mb = lambda wildcards, attempt : [6000, 16000, 32000, 64000, 0, 0][attempt-1]
    benchmark: "output/{config_batch}/benchmarks/benchmark.make_database.{config_batch}.tsv"
    shell: """

        mkdir -p output/{config_batch}/
        >&2 echo "Concatenating database ..."
        cat {input.glob} > output/{config_batch}/cat_database_sources.faa

        mkdir -p output/{config_batch}/
        cd output/{config_batch}/

        {params.philosopher} workspace --init
        
        {params.philosopher} database --help

        # https://github.com/Nesvilab/philosopher/wiki/Database
        {params.philosopher} database \
            --custom cat_database_sources.faa \
            --contam 

        mv *-decoys-contam-cat_database_sources.faa.fas philosopher_database.fas # rename database file.
        rm cat_database_sources.faa # remove unneccessary .faa file.

        {params.philosopher} workspace --clean
        ls -la

    """


# def memory_msfragger(wildcards, attempt):
#     return attempt * 100000



# Run the headless fragpipe command

rule fragpipe:
    input: 
        database = "output/{config_batch}/philosopher_database.fas",
    output:
        flag = touch("output/{config_batch}/fragpipe_done.flag")
    params:
        fragpipe_workflow = f"output/{config_batch}/msfragger/fragpipe_modified.workflow",

        #msfragger_jar = config["msfragger_jar"],
        #fragpipe_base = config["fragpipe_base"],
        n_splits = 8,
        manifest = manifest.to_csv(path_or_buf=None, sep = "\t", index=False, header=False),
    threads: 8
    resources:
        #partition = "bigmem",
        #mem_mb = 40000, 
        runtime = "24:00:00"
    conda: "envs/openjdk_python.yaml"
    benchmark: "output/{config_batch}/benchmarks/benchmark.fragpipe.{config_batch}.tsv"
    shell: """


        >&2 echo "Create manifest ..."
        # It is a matter of outputting the metadata table in the correct manner.
        echo '''{params.manifest}''' > {output.manifest}


        >&2 echo "Create workflow ..."
        # Write the missing dynamic lines to the workflow, depending on the current setup.
        # Copy and modify parameter file.
        cp arturos_workflow/LFQ-MBR_TimsTOF_edit.workflow {params.fragpipe_workflow}
        echo "" >> {params.fragpipe_workflow}
        echo "num_threads = {threads}" >> {params.fragpipe_workflow}
        echo "database_name = {input.database}" >> {params.fragpipe_workflow}
        echo "output_location = output/{wildcards.config_batch}/msfragger/" >> {params.fragpipe_workflow}
        echo "" >> {params.fragpipe_workflow}



        # Debug: Let's stop it here for now.
        exit 0        

        >&2 echo "Fragpipe ..."
        # https://fragpipe.nesvilab.org/docs/tutorial_headless.html
        time fragpipe \
            --headless \
            --workflow LFQ-MBR.workflow \
            --manifest $tabmanifest.manifest.FragPipe.fp-manifest \
            --workdir $workdir \
            --ram {resources.mem_mb} \
            --threads {threads}

        >&2 echo "Fragpipe must have exited successfully ..."
    

    """





# rule msfragger:
#     input:
#         database = "output/{config_batch}/philosopher_database.fas",
#         d_files = ("output/{config_batch}/msfragger/" + df["barcode"]).tolist(),
#         #linked_flag = "output/{config_batch}/msfragger/link_input.done",
#     output: 
#         pepXMLs = "output/{config_batch}/msfragger/" + df["basename"] + ".pepXML",
#     threads: 16
#     params:
#         msfraggerparams = f"output/{config_batch}/msfragger/msfragger.params",
#         #msfraggerparams = "msfragger.params",
##         msfragger_jar = config["msfragger_jar"],
#         fragpipe_base = config["fragpipe_base"],
#         n_splits = 8 # splits is a bad name. Having 8 groups means having 7 splits, right?
#     resources:
#         partition = "bigmem",
#         #mem_mb = 70000, # Was 500000 before I used the split script
##         mem_mb = 40000,
#         #mem_mb = lambda wildcards, attempt : attempt * 100000
#         #runtime = "23:59:59" 26 samples done in 24 hours
#         #runtime = "6-00:00:00" #prev
#         runtime = "24:00:00"
#     conda: "envs/openjdk_python.yaml"
#     benchmark: "output/{config_batch}/benchmarks/benchmark.msfragger.{config_batch}.tsv"
#     shell: """

        


#         # Copy and modify parameter file.
#         cp msfragger_default.params {params.msfraggerparams}
#         echo "" >> {params.msfraggerparams}
#         echo "num_threads = {threads}" >> {params.msfraggerparams}
#         echo "database_name = {input.database}" >> {params.msfraggerparams}
#         echo "output_location = output/{wildcards.config_batch}/msfragger/" >> {params.msfraggerparams}
#         echo "" >> {params.msfraggerparams}


#         # This is the non-standard idiosyncratic msfragger-agnostic usage of the split script:
#         # python3 msfragger_pep_split.pyz 3 "java -Xmx64g -jar" msfragger.jar fragger.params *.mzML
#         #                                 ^ num_parts_str
#         #                                   ^ jvm_cmd_str
#         #                                                       ^ msfragger_jar_path_str
#         #                                                                      ^ param_path_str
#         #                                                                                    ^ *infiles_str


#         # Call msfragger with the split database script
#         # Ode to https://github.com/Nesvilab/MSFragger/issues/180#issuecomment-938323065
#         python {params.fragpipe_base}/tools/msfragger_pep_split.py {params.n_splits} \
#             "java -Xmx64g -jar" \
#             {params.msfragger_jar} \
#             {params.msfraggerparams} \
#             {input.d_files} 


#     """




# rule workspace:
#     input: 
#         pepXMLs = "output/{config_batch}/msfragger/" + df["basename"] + ".pepXML",
#         database = "output/{config_batch}/philosopher_database.fas",
#     output:
#         prot_xml = "output/{config_batch}/workspace/proteinprophet.prot.xml",
#         psm = "output/{config_batch}/workspace/psm.tsv",
#     conda: "envs/openjdk.yaml"
#     resources:
#         mem_mb = 32000
#     params:
#         philosopher = config["philosopher_executable"],
#         decoyprefix = "rev_",
#     benchmark: "output/{config_batch}/benchmarks/benchmark.workspace.{config_batch}.tsv"
#     shell: """

#         >&2 echo "mkcd ..."
#         # Make and enter workspace directory
#         mkdir -p output/{config_batch}/workspace
#         cd output/{config_batch}/workspace


#         >&2 echo "Workspace ..."
#         # Run PeptideProphet, ProteinProphet, and FDR filtering with Philosopher
#         {params.philosopher} workspace --clean
#         {params.philosopher} workspace --init


#         # I don't think I need to do anything with the database. It should already be annotated right?
#         # BTW, what the fuck happens when you annotate anyway?
#         >&2 echo "Database ..."
#         {params.philosopher} database --annotate ../../../{input.database} --prefix {params.decoyprefix}
#         >&2 ls -la

        


#         # Closed search
#         >&2 echo "Peptideprophet ..."
#         {params.philosopher} peptideprophet \
#             --nonparam --expectscore --decoyprobs --ppm --accmass \
#             --decoy {params.decoyprefix} \
#             --database ../../../{input.database} \
#             ../msfragger/*.pepXML # Take the pepXMLs directly from msfragger.

        




#         >&2 echo "Proteinprophet ..."
#         {params.philosopher} proteinprophet \
#             --maxppmdiff 2000000 \
#             --output proteinprophet \
#             ../msfragger/*.pep.xml

        

#         # closed or non-specific closed search
#         >&2 echo "Filter ..."
#         {params.philosopher} filter \
#             --sequential --razor --mapmods \
#             --tag {params.decoyprefix} \
#             --pepxml ../msfragger/*.pepXML \
#             --protxml ./proteinprophet.prot.xml 

        
        
#         >&2 echo "Reports ..."
#         # Generate reports.
#         {params.philosopher} report
#         {params.philosopher} workspace --clean

#     """




# rule ionquant:
#     input:
#         prot_xml = "output/{config_batch}/workspace/proteinprophet.prot.xml",
#         psm = "output/{config_batch}/workspace/psm.tsv",
#         pepXMLs = "output/{config_batch}/msfragger/" + df["basename"] + ".pepXML",
#     output: 
#         final_flag = touch("output/{config_batch}/final.flag"),
#         peptide = "output/{config_batch}/workspace/peptide.tsv", # These files will be missing if there is no callable output.
#         protein = "output/{config_batch}/workspace/protein.tsv",
#         ion = "output/{config_batch}/workspace/ion.tsv",
#     threads: 8
#     resources:
#         mem_mb = 128000
#     params:
#         ionquant_jar = config["ionquant_jar"],
#     conda: "envs/openjdk.yaml"
#     benchmark: "output/{config_batch}/benchmarks/benchmark.ionquant.{config_batch}.tsv"
#     shell: """


#         >&2 echo "Ionquant ..."
#         java \
#             -Xmx128G \
#             -jar {params.ionquant_jar} \
#             --threads {threads} \
#             --psm {input.psm} \
#             --specdir output/{config_batch}/msfragger \
#             {input.pepXMLs}
            

#         # protein output might be empty if there is no matches
#         touch {output.protein}

    
#     """


# #rule report:




# onsuccess:
#     shell("echo -n \"All good :)\ndepth 2 tree below:\n\"; tree -L 2 output/{config_batch}/")

# onerror:
#     shell("echo -n \"ERROR :(\ndepth 2 tree below:\n\"; tree -L 2 output/{config_batch}/")


# print("*/") # This is a dot-language specific comment close tag that helps when you export the workflow as a graph



