# MS_pipeline1



This is a very simple pipeline to transform .d files with an appropriate database into sane protein calls with abundances.


The workflow is based on msfragger and aims to keep a tidy tree of output files.


## Installation

1) Prerequisites:

A conda/miniconda3 environment


2) Install snakemake by following the [instructions in the official documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

3) Clone this repo

4) Set up the profiles/slurm/ configuration so that it matches your execution environment. There is also a profile for local execution without a job management system (profiles/local/)



## Usage

#### 1) Update config.yaml
Because nesvilab doesn't make their executables easily publicly available, you need to tell the pipeline where to find them on your system. Update addresses for the keys `philosopher_executable`, `msfragger_jar` and `ionquant_jar` which can be downloaded [here](https://github.com/nesvilab/philosopher/releases/latest), [here](https://github.com/Nesvilab/MSFragger/wiki/Preparing-MSFragger#Downloading-MSFragger) and [here](https://github.com/Nesvilab/IonQuant#download), respectively. 


Currently the pipeline only supports the input of .d-files ([agilent/bruker](https://en.wikipedia.org/wiki/Mass_spectrometry_data_format#Proprietary_formats)). Create an item in batch_parameters where you define key `d_base` which is the base directory where all .d-files reside. Define key `database_glob` which is a path (or glob) to the fasta-amino acid files that you want to include in the target protein database.

Define items under the `samples` key which link sample names to the .d-files.

Lastly, set the root `batch` key to point at the batch that you want to run.

#### 2) Run

Finally, run the pipeline in your command line with:
```
$ snakemake --profile profiles/slurm/ 
```

Below is an example graph for a run involving six samples.
<img width="1466" alt="Screenshot 2022-05-06 at 11 18 15" src="https://user-images.githubusercontent.com/5913696/167104044-77e8d2d7-20cd-4334-91bc-a730b2867e41.png">






