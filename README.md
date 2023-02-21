# MS-pipeline1


```
                                             ______________ 
                                            < MS-pipeline1 >
                                             -------------- 
                                                          \ 
                             ___......__             _     \
                         _.-'           ~-_       _.=a~~-_  
 --=====-.-.-_----------~   .--.       _   -.__.-~ ( ___===>
               '''--...__  (    \ \\\ { )       _.-~        
                         =_ ~_  \\-~~~//~~~~-=-~            
                          |-=-~_ \\   \\                    
                          |_/   =. )   ~}                   
                          |}      ||                        
                         //       ||                        
                       _//        {{                        
                    '='~'          \\_    =                 
                                    ~~'    
```



This pipeline takes 1) a list of .d files and 2) a list of fasta-amino acid files and outputs sane protein calls with abundances. It uses philosopher database and fragpipe to do the job. The snakemake pipeline maintains a nice output file tree.

### Why you should use this pipeline

Because it makes sure that all outputs are updated when you change input-parameters. It also yells at you if something fails, and hopefully makes it a bit easier to find the error.


## Installation
1) Prerequisites:
  - Preferably a HPC system, or a beefy local workstation.
  - An anaconda or [miniconda3](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) package manager on that system.

2) Clone this repo on the HPC/workstation where you want to work.
   ```
   git clone https://github.com/cmkobel/MS-pipeline1.git
   ```

3) If you don't already have an environment with snakemake and mamba installed, use the following command to install a "snakemake" environment with the bundled environment file:
   ```
   conda env create -f environment.yaml
   ```

   This environment can then be activated by typing `conda activate snakemake`


3) If needed, tweak the profiles/slurm/<file> configuration so that it matches your execution environment. There is a profile for local execution without a job management system (`profiles/local/`) as well as a few profiles for different HPC environments like PBS and SLURM. 
  

## Usage

#### 1) Update config.yaml

The file config.yaml contains all the parameters needed for this pipeline to run. You should update the parameters so they reflect which samples you wish to process.

Because nesvilab do not make their executables immediately publicly available, you need to tell the pipeline where to find them on your system. Update addresses for the keys `philosopher_executable`, `msfragger_jar`, `ionquant_jar` and `fragpipe_jar` which can be downloaded [here](https://github.com/nesvilab/philosopher/releases/latest), [here](https://github.com/Nesvilab/MSFragger/wiki/Preparing-MSFragger#Downloading-MSFragger), [here](https://github.com/Nesvilab/IonQuant#download) and [here](https://github.com/Nesvilab/FragPipe/releases), respectively. 


Currently the pipeline only supports the input of .d-files ([agilent/bruker](https://en.wikipedia.org/wiki/Mass_spectrometry_data_format#Proprietary_formats)). Create an item in batch_parameters where you define key `d_base` which is the base directory where all .d-files reside. Define key `database_glob` which is a path (or glob) to the fasta-amino acid files that you want to include in the target protein database.

Define items under the `samples` key which link sample names to the .d-files.

Lastly, set the root `batch` key to point at the batch that you want to run.

#### 2) Run

Finally, run the pipeline in your command line with:
```
$ snakemake --profile profiles/slurm/
```

Below is an example graph for a run involving two samples.

<img width="350" src="https://user-images.githubusercontent.com/5913696/212691889-f33bb808-6675-4f66-8440-84f9c215f83b.png">


## Future

This pipeline might involve an R-markdown performing trivial QC.
Also, a test data set that accelerates the development cycle. 🚴‍♀️






