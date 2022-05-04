# MS_pipeline1



This is a very simple pipeline to transform .d files with an appropriate database into sane protein calls with abundances.


The workflow is based on msfragger and aims to keep a tidy tree of output files.


## Installation

Prerequisites: 
  -  A conda/miniconda3 environment


Install snakemake by following the [instructions in the official documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

Set up the profiles/slurm/ configuration so that it matches your execution environment. There is also a profile for local execution without a job management system (profiles/local/)



## Usage

#### 1) Update config.yaml
Because nesvilab doesn't make their executables easily publicly available, you need to tell the pipeline where to find them on your system. Update addresses to [philosopher_executable](https://github.com/nesvilab/philosopher/releases/latest), [msfragger_jar](https://github.com/Nesvilab/MSFragger/wiki/Preparing-MSFragger#Downloading-MSFragger) and [ionquant_jar](https://github.com/Nesvilab/IonQuant#download). 


Currently the pipeline only supports the input of .d-files ([agilent/bruker](https://en.wikipedia.org/wiki/Mass_spectrometry_data_format#Proprietary_formats)). Create an item in batch_parameters where you define `d_base` which is the base directory where all .d-files reside. Define `database_glob` which is a path (or glob) to the fasta-amino acid files that you want to include in the target protein database.

Define items under the `samples` key which link sample names to the .d-files.

Lastly, set the root `batch` key to point at the batch that you want to run.

#### 2) Run

Finally, run the pipeline in your command line with:
```
$ snakemake --profile profiles/slurm/ 
```

Below is an example graph for a run involving two samples.
<img width="1465" alt="Screenshot 2022-05-04 at 15 41 50" src="https://user-images.githubusercontent.com/5913696/166693733-063ff92c-cd72-40fd-b49a-48d0fb1f897a.png">





