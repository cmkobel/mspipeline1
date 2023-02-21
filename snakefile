# snakemake --profile profiles/slurm-sigma2-saga
# I'm experiencing some major problems with the temporary directories that i might as well fix. It seems to percolate through when I have a high amount of samples. Basically, all the jobs that use a program that uses the workspace, need to be in the same rule. Silly, but that is how it is.

# In this branch I'm not at all screwing around. I'm closely following this tutorial:
# https://fragpipe.nesvilab.org/docs/tutorial_linux.html


__author__ =  "Carl Mathias Kobel & Arturo Vera Ponce De Leon"
__version__ = "v2.0.0"

# changelog:

# v2.0.0: Using Arturos pipeline


# TODO: refactor some variables and names. msfragger is not a good name for the final output dir, as it is rather fragpipe that is being called. Also, the .d files could potentially be linked to a different directory, for instance a temporary userwork dir or just a different directory, maybe called linked_inputs/.


import glob
import pandas as pd
import re
import pathlib

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


absolute_output_dir = str(pathlib.Path("output/").absolute())


# Present configuration
print(f"      abs out dir is: '{absolute_output_dir}'")
print(f"        config_batch: '{config_batch}'")
print(f"       config_d_base: '{config_d_base}'")
print(f"config_database_glob: '{config_database_glob}:'")

if len(config_database_glob_read) < 1:
    raise Exception("No glob targets in config_database_glob") # Not tested yet.

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
manifest['path'] = absolute_output_dir + "/" + config_batch + "/samples/" + manifest['path'] # Instead of using bash realpath
##manifest["path"] = "output/" + config_batch + "/msfragger/" + manifest["path"] # But then I realized that I might not need to point absolutely anyway..
#manifest["path"] = manifest["path"]
manifest["experiment"] = "experiment" # Experiment (can be empty, alphanumeric, and _) #  IonQuant with MBR requires designating LCMS runs to experiments. If in doubt how to resolve this error, just assign all LCMS runs to the same experiment name.
manifest["bioreplicate"] = "" # Bioreplicate (can be empty and integer)
manifest["data_type"] = "DDA" # Data type (DDA, DIA, GPF-DIA, DIA-Quant, DIA-Lib)
print(manifest)
print("//")






# Define workflow targets
rule all:
    input:
        metadata = f"output/{config_batch}/metadata.tsv",
        copy_input = f"output/{config_batch}/samples/copy_samples.done",
        make_database = f"output/{config_batch}/philosopher_database.fas", 
        fragpipe = f"output/{config_batch}/fragpipe_done.flag",
        



# Save some metadata about inputs for good measure.
rule metadata:
    output: "output/{config_batch}/metadata.tsv"
    params: dataframe = df.to_csv(None, index_label = "index", sep = "\t"),
    shell: """

        echo '''{params.dataframe}''' > {output}
    
    """

# Link input links or copies the input data to a specific directory. Long term, this should be on the fastest possible disk ie. userwork.
rule copy_samples: # Or place_samples, or copy_samples
    output:
        flag = touch("output/{config_batch}/samples/copy_samples.done"), # Just used to keep track of stuff - probably a bad habit.
        dir = directory("output/{config_batch}/samples"), # Why is this necessary?
        d_files = directory("output/{config_batch}/samples/" + df["barcode"]), # Bound for fragpipe.
    params:
        d_files = (config_d_base + "/" + df["barcode"]).tolist(), # Problem is that snakemake doesn't like directories as inputs, so I think it is better to define it as a param.
    benchmark: "output/{config_batch}/benchmarks/benchmark.copy_samples.{config_batch}.tsv"
    shell: """
        
        #ln -s {params.d_files} {output.dir} # Should be deprecated?
        cp -vr {params.d_files} {output.dir}

    """


# make_database cats all the amino acid fastas together and runs philosopher database on it
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

        # https://github.com/Nesvilab/philosopher/wiki/Database
        {params.philosopher} database \
            --custom cat_database_sources.faa \
            --contam 

        mv *-decoys-contam-cat_database_sources.faa.fas philosopher_database.fas # rename database file.
        rm cat_database_sources.faa # remove unneccessary .faa file.

        {params.philosopher} workspace --clean

        >&2 echo "Statistics ..."
        n_records=$(cat philosopher_database.fas | grep -E ">" | wc -l)
        echo -e "n_records_in_db\t$n_records" > n_records.tsv

    """


# def memory_msfragger(wildcards, attempt):
#     return attempt * 100000



# Run the fragpipe in headless. Define manifest and workflow on the fly.
rule fragpipe:
    input: 
        database = "output/{config_batch}/philosopher_database.fas",
    output:
        flag = touch("output/{config_batch}/fragpipe_done.flag"),
        manifest = "output/{config_batch}/fragpipe/{config_batch}.manifest"
    params:
        manifest = manifest.to_csv(path_or_buf=None, sep = "\t", index=False, header=False, lineterminator=False),
        fragpipe_workflow = f"output/{config_batch}/fragpipe/fragpipe_modified.workflow",
        n_splits = 8,

        fragpipe_executable = config["fragpipe_executable"],
        msfragger_jar = config["msfragger_jar"],
        ionquant_jar = config["ionquant_jar"],
        philosopher_executable = config["philosopher_executable"],

        fragpipe_workdir = "output/{config_batch}/fragpipe",
    threads: 8
    resources:
        #partition = "bigmem", # When using more than 178.5 GB at sigma2/saga
        mem_mb = 32768, # Arturo uses 150GB in bigmem with 12 threads.
        runtime = "24:00:00",
    conda: "envs/openjdk_python.yaml"
    #conda: "envs/openjdk_python_extra.yaml"
    benchmark: "output/{config_batch}/benchmarks/benchmark.fragpipe.tsv"
    shell: """

        echo "Create manifest ..."
        echo '''{params.manifest}''' > {output.manifest} 
        tail {output.manifest}


        echo "Create workflow ..."
        # Copy and modify parameter file with dynamic content.
        cp assets/fragpipe_workflows/LFQ-MBR.workflow {params.fragpipe_workflow}
        echo "" >> {params.fragpipe_workflow}
        echo "num_threads = {threads}" >> {params.fragpipe_workflow}
        echo "database_name = {input.database}" >> {params.fragpipe_workflow}
        echo "database.db-path = {input.database}" >> {params.fragpipe_workflow}
        echo "msfragger.misc.slice-db = {params.n_splits}" >> {params.fragpipe_workflow}
        echo "output_location = {params.fragpipe_workdir}" >> {params.fragpipe_workflow}
        echo "" >> {params.fragpipe_workflow}
        tail {params.fragpipe_workflow}


        # Convert mem_mb into gb
        mem_gb=$(({resources.mem_mb}/1024-2)) # Because there is some overhead, we subtract a few GBs. Everytime fragpipe runs out of memory, I subtract another one: that should be more effective than doing a series of tests ahead of time.
        echo "Fragpipe will be told to not use more than $mem_gb GB. In practice it usually uses a bit more."

        echo "Fragpipe ..."
        # https://fragpipe.nesvilab.org/docs/tutorial_headless.html
        {params.fragpipe_executable} \
            --headless \
            --workflow {params.fragpipe_workflow} \
            --manifest {output.manifest} \
            --workdir {params.fragpipe_workdir} \
            --ram $mem_gb \
            --threads {threads} \
            --config-msfragger {params.msfragger_jar} \
            --config-ionquant {params.ionquant_jar} \
            --config-philosopher {params.philosopher_executable}

        # Possibly do some output validation here.

    """



onstart: 
    shell("mkdir -p logs/old/; mv logs/*.log logs/old/ 2> /dev/null || exit 0") # Put old logs aside
    shell("find output/ > .onstart.txt 2> /dev/null || exit 0")

onsuccess:
    print("onsuccess: The following files were created:")
    shell("find output/ > .onsuccess.txt && diff .onstart.txt .onsuccess.txt || exit 0")


print("*/") # This is a dot-language specific comment close tag that helps when you export the workflow as a graph



