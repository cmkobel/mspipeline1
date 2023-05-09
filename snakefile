# conda activate snakemake_7_24
# snakemake --profile profiles/local
# snakemake --profile profiles/slurm-sigma2-saga


__author__ =  "Carl Mathias Kobel & Arturo Vera Ponce De Leon"
__version__ = "v2.0.0"

# changelog:
# v2.0.0: Using Arturos pipeline, resulting in the first version that actually works.



import glob
import pandas as pd
import re
import pathlib

print("/*                                                                               ") # Helps with outputting to dot format.
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
try:
    configfile: "config.yaml"
except WorkflowError:
    print("Error: Please make sure that the config.yaml file is present and correctly formatted. You can use the config_template.yaml file as a starting point.")

config_batch = config["batch"]
config_d_base = config["batch_parameters"][config_batch]["d_base"]
config_database_glob = config["batch_parameters"][config_batch]["database_glob"]
config_database_glob_read = glob.glob(config_database_glob)
config_samples = config["batch_parameters"][config_batch]["samples"]


absolute_output_dir = str(pathlib.Path("output/").absolute())


# Present configuration
print(f"        abs out dir is: '{absolute_output_dir}'")
print(f"          config_batch: '{config_batch}'")
print(f"         config_d_base: '{config_d_base}'")
print(f"  config_database_glob: '{config_database_glob}:'")

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


#print("manifest:")
manifest = pd.DataFrame(data = {'path': config_samples.values(), 'experiment': config_samples.keys()})
manifest['path'] = absolute_output_dir + "/" + config_batch + "/samples/" + manifest['path'] # Instead of using bash realpath
##manifest["path"] = "output/" + config_batch + "/msfragger/" + manifest["path"] # But then I realized that I might not need to point absolutely anyway..
#manifest["path"] = manifest["path"]
#manifest["experiment"] = "experiment" # Experiment (can be empty, alphanumeric, and _) #  IonQuant with MBR requires designating LCMS runs to experiments. If in doubt how to resolve this error, just assign all LCMS runs to the same experiment name.
manifest["bioreplicate"] = "" # Bioreplicate (can be empty and integer)
manifest["data_type"] = "DDA" # Data type (DDA, DIA, GPF-DIA, DIA-Quant, DIA-Lib)
#print(manifest); print("//")






# Define workflow targets
rule all:
    input:
        metadata = f"output/{config_batch}/metadata.tsv",
        copy_input = f"output/{config_batch}/samples/copy_samples.done",
        make_database = f"output/{config_batch}/philosopher_database.faa", 
        fragpipe = f"output/{config_batch}/fragpipe_done.flag",
        report = f"output/{config_batch}/report_MS-pipeline1_{config_batch}.html",
        zip_ = f"output/{config_batch}/MS-pipeline1_{config_batch}.zip",

        
        



# Save some metadata about inputs for good measure.
rule metadata:
    output: "output/{config_batch}/metadata.tsv"
    params: dataframe = df.to_csv(None, index_label = "index", sep = "\t"),
    shell: """

        echo '''{params.dataframe}''' > {output}

        # TODO: Write something clever to the benchmarks/ directory, so we can infer relationship between hardware allocations, input size and running time.
    
    """

# Link input links or copies the input data to a specific directory. Long term, this should be on the fastest possible disk ie. userwork.
# TODO: Add the metadata as input to this rule.
rule copy_samples: # Or place_samples, or copy_samples
    output:
        flag = touch("output/{config_batch}/samples/copy_samples.done"), # Used to keep fragpipe waiting.
        dir = directory("output/{config_batch}/samples"), # Why is this necessary?
        d_files = directory("output/{config_batch}/samples/" + df["barcode"]), # Bound for fragpipe. Update, but fragpipe uses the flag instead?
    params:
        d_files = (config_d_base + "/" + df["barcode"]).tolist(), # Problem is that snakemake doesn't like directories as inputs, so I think it is better to define it as a param.
    benchmark: "output/{config_batch}/benchmarks/benchmark.copy_samples.{config_batch}.tsv"
    shell: """

        cp -vr {params.d_files} {output.dir}

        # Enable editing of these files
        chmod -R 775 output/

    """


# make_database cats all the amino acid fastas together and runs philosopher database on it
rule make_database:
    input:
        #glob = [glob.glob(config_database_glob)], # Not sure why this one was inside bracket/list?
        glob = config_database_glob_read
    output:
        database = "output/{config_batch}/philosopher_database.faa",
    params:
        philosopher = config["philosopher_executable"],
    retries: 4 # This is some black magic voodoo shit.
    resources:
        mem_mb = lambda wildcards, attempt : [16000, 32000, 64000, 128000, 175000, 0, 0][attempt-1], # Attempt starts from 1? Confirm:
        #mem_mb = 64000,
        runtime = "24h",
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

        echo "Existing database pattern-matched files:"
        ls *-decoys-contam-cat_database_sources.faa.fas

        mv *-decoys-contam-cat_database_sources.faa.fas philosopher_database.faa # rename database file.
        rm cat_database_sources.faa # remove unneccessary .faa file.

        {params.philosopher} workspace --clean

        


    """



rule db_stats:
    input: 
        database = "output/{config_batch}/philosopher_database.faa",
    output: 
        db_stats_seqkit = "output/{config_batch}/db_stats_seqkit.tsv",
    params: 
        config_database_glob_read = config_database_glob_read,
    resources: 
        runtime = "1h",
        mem_mb = 256,
    conda: "envs/seqkit.yaml"
    shell: """

        # TODO: Database length in basepairs?
        seqkit stats \
            --tabular \
            {params.config_database_glob_read} {input.database} \
        > {output.db_stats_seqkit}

    """





# def memory_msfragger(wildcards, attempt):
#     return attempt * 100000



# Run the fragpipe in headless. Define manifest and workflow on the fly.
# Potentially we could implement retries into this rule. But since there are soo many ways fragpipe can fail, and I don't want to waste many cpu hours, I'm abstaining from implementing.
rule fragpipe:
    input: 
        copy_samples = "output/{config_batch}/samples/copy_samples.done",
        database = "output/{config_batch}/philosopher_database.faa",
    output:
        flag = touch("output/{config_batch}/fragpipe_done.flag"),
        manifest = "output/{config_batch}/fragpipe/{config_batch}.manifest",
        fragpipe_workflow = "output/{config_batch}/fragpipe/fragpipe_modified.workflow",

        #stats = "output/{config_batch}/fragpipe/fragpipe_stats.tsv", # Moved to rule fragpipe_stats

        # final results:
        #psm = "output/{config_batch}/fragpipe/experiment/psm.tsv", # A file will be created for each experiment, and really R should read each of those directly.
        
        final_ion = "output/{config_batch}/fragpipe/combined_ion.tsv",
        final_peptide = "output/{config_batch}/fragpipe/combined_peptide.tsv",
        final_protein = "output/{config_batch}/fragpipe/combined_protein.tsv",

        fragpipe_stdout = "output/{config_batch}/fragpipe/fragpipe.out.log"
    params:
        manifest = manifest.to_csv(path_or_buf=None, sep="\t", index=False, header=False), # This is a csv formatted string 
        #original_fragpipe_workflow = "assets/fragpipe_workflows/LFQ-MBR.workflow", # The path to the workflow that specifies the type of analysis
        original_fragpipe_workflow = "assets/fragpipe_workflows/LFQ-MBR_carl_no_overwrite.workflow", # The path to the workflow that specifies the type of analysis. Honestly, I don't think it matters that you just overwrite the settings in the original.
        slice_db = 16, # The number of database splits that fragpipe (msfragger) should perform. 

        fragpipe_executable = config["fragpipe_executable"],
        msfragger_jar = config["msfragger_jar"],
        ionquant_jar = config["ionquant_jar"],
        philosopher_executable = config["philosopher_executable"],

        fragpipe_workdir = "output/{config_batch}/fragpipe", # Bound for fragpipe --workdir
    threads: 16 # 8 for testing
    resources:
        #partition = "bigmem", # When using more than 178.5 GB at sigma2/saga
        #mem_mb = 32000, # for testing
        #mem_mb = 190000, # Some people like to use 150GB in bigmem with 12 threads.
        mem_mb = 453632, # Giant swap
        runtime = "36h",
    #conda: "envs/openjdk_python.yaml"
    conda: "envs/openjdk_python_extra.yaml" # TODO: Use this file, I checked it already, and you just have to install pyopenms manually. Don't want to use a previous version of python (e.g. 3.9) just to have easypqp installed, as it seems like some people do not have it too.
    benchmark: "output/{config_batch}/benchmarks/benchmark.fragpipe.{config_batch}.tsv"
    shell: """

        echo "Create manifest ..."
        echo '''{params.manifest}''' > {output.manifest} 
        tail {output.manifest}


        echo "Modifying workflow with runtime parameters ..." # TODO: Check if it matters to overwrite or not.
        # Copy and modify parameter file with dynamic content.
        cp {params.original_fragpipe_workflow} {output.fragpipe_workflow}
        echo -e "\n# Added by mspipeline1 in rule fragpipe in snakefile below ...\n" >> {output.fragpipe_workflow}

        echo "num_threads={threads}" >> {output.fragpipe_workflow}
        echo "database_name={input.database}" >> {output.fragpipe_workflow}
        echo "database.db-path={input.database}" >> {output.fragpipe_workflow}
        
        echo "output_location={params.fragpipe_workdir}" >> {output.fragpipe_workflow}
        
        # These settings minimize memory usage. 
        echo "msfragger.misc.slice-db={params.slice_db}" >> {output.fragpipe_workflow} # Default 1
        echo "msfragger.calibrate_mass=0" >> {output.fragpipe_workflow} # Default 2
        echo "msfragger.digest_max_length=35" >> {output.fragpipe_workflow} # Default 50
        # echo "msfragger.allowed_missed_cleavage_1=1" >> {output.fragpipe_workflow} # Default 2
        # echo "msfragger.allowed_missed_cleavage_2=1" >> {output.fragpipe_workflow} # Default 2
        
        echo "" >> {output.fragpipe_workflow}
        tail {output.fragpipe_workflow}


        # Convert mem_mb into gb
        mem_gb=$(({resources.mem_mb}/1024-2)) # Because there is some overhead, we subtract a few GBs. Everytime fragpipe runs out of memory, I subtract another one: that should be more effective than doing a series of tests ahead of time.
        echo "Fragpipe will be told not to use more than $mem_gb GB. In practice it usually uses a bit more."

        echo "Fragpipe ..."
        # https://fragpipe.nesvilab.org/docs/tutorial_headless.html
        {params.fragpipe_executable} \
            --headless \
            --workflow {output.fragpipe_workflow} \
            --manifest {output.manifest} \
            --workdir {params.fragpipe_workdir} \
            --ram $mem_gb \
            --threads {threads} \
            --config-msfragger {params.msfragger_jar} \
            --config-ionquant {params.ionquant_jar} \
            --config-philosopher {params.philosopher_executable} \
        | tee {output.fragpipe_stdout} # Write the log, so we can later extract the number of "scans"

       

    """

# I moved some of these stats out just to make the debugging easier. 
# This could have been tailing the fragpipe, but I just think it is easier to develop it like this. 
rule fragpipe_stats:
    input: 
        fragpipe_stdout = "output/{config_batch}/fragpipe/fragpipe.out.log",
    output:
        scans = "output/{config_batch}/fragpipe/stats_fragpipe_scans.tsv",
        #psms = "output/{config_batch}/psms.tsv", # Better to let R read these files from the zip directory.
    params:
        fragpipe_workdir = "output/{config_batch}/fragpipe", # Bound for fragpipe --workdir
    shell: """

        # Extract scans from the fragpipe stdout log. Will later be compared to the individual psm files.
        grep -E ": Scans = [0-9]+" {input.fragpipe_stdout} \
        > {output.scans}


    """
        



# Rename this to report: Do the report first, then zip the report with its outputs.
rule report:
    input: 
        "output/{config_batch}/metadata.tsv",
        "output/{config_batch}/db_stats_seqkit.tsv",
        "output/{config_batch}/fragpipe/{config_batch}.manifest",
        "output/{config_batch}/fragpipe/fragpipe_modified.workflow",
        "output/{config_batch}/fragpipe/stats_fragpipe_scans.tsv",
        "output/{config_batch}/fragpipe/combined_ion.tsv",
        "output/{config_batch}/fragpipe/combined_peptide.tsv",
        "output/{config_batch}/fragpipe/combined_protein.tsv",
        "output/{config_batch}/fragpipe/fragpipe.out.log",
    output:
        report = "output/{config_batch}/report_MS-pipeline1_{config_batch}.html",
        idrate = "output/{config_batch}/identification_rate.tsv",
    benchmark: "output/{config_batch}/benchmarks/benchmark.report.{config_batch}.tsv"
    conda: "envs/r-markdown.yaml"
    resources:
        runtime = "4h",
        mem_mb = 4096,
    shell: """

        cp scripts/QC.Rmd rmarkdown_template.rmd
        Rscript -e 'rmarkdown::render("rmarkdown_template.rmd", "html_document", output_file = "{output.report}", knit_root_dir = "output/{config_batch}/", quiet = T)' 
        rm rmarkdown_template.rmd

    """


rule zip:
    input:
        "output/{config_batch}/metadata.tsv",
        "output/{config_batch}/db_stats_seqkit.tsv",
        "output/{config_batch}/fragpipe/{config_batch}.manifest",
        "output/{config_batch}/fragpipe/fragpipe_modified.workflow",
        "output/{config_batch}/fragpipe/stats_fragpipe_scans.tsv",
        "output/{config_batch}/fragpipe/combined_ion.tsv",
        "output/{config_batch}/fragpipe/combined_peptide.tsv",
        "output/{config_batch}/fragpipe/combined_protein.tsv",
        "output/{config_batch}/report_MS-pipeline1_{config_batch}.html",
        "output/{config_batch}/identification_rate.tsv",
        "scripts/QC.Rmd"
    output:
        "output/{config_batch}/MS-pipeline1_{config_batch}.zip",
    resources:
        runtime = "1h",
    shell: """
    
        zip {output} {input} 

    """


onstart: 
    shell("mkdir -p logs/old/; (mv logs/*.log logs/old/ 2> /dev/null || exit 0) &") # Put old logs aside
    shell("find output/ > .onstart.txt 2> /dev/null || exit 0")

onsuccess:
    print("onsuccess: The following (first 100) files were created:")
    #shell("find output/ > .onsuccess.txt && diff .onstart.txt .onsuccess.txt | head -n 100 || exit 0")
    shell(""" diff .onstart.txt .onsuccess.txt > .diff.txt || head -n 10 .diff.txt; echo "$(cat .diff.txt | wc -l) files total" """)


print("*/") # This is a dot-language specific comment close tag that helps when you export the workflow as a graph



# TODO: in the benchmark directory, there should also be information about number of proteins processed, and number of scans in the tims-input. This way we might better be able to figure the resource requirements for bigger future projects.