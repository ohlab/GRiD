# GRiD
Growth Rate Index (GRiD) measures bacterial growth rate from reference genomes (including draft quality genomes) and metagenomic bins at ultra-low sequencing coverage (> 0.2x). GRiD is described in 

Emiola and Oh (2018) "High throughput in situ metagenomic measurement of bacterial replication at ultra-low sequencing coverage." *Nature Communications*  (https://www.nature.com/articles/s41467-018-07240-8).

GRiD algorithm consists of two modules;

1. < single > - which is applicable for growth analysis involving a single reference genome
2. < multiplex > - for the high-throughput growth analysis of all identified bacteria in a sample. Prior knowledge of microbial composition is not required. To use this module, you must download the GRiD database.
 
The comprehensive GRiD database consists of 32,819 representative bacteria genomes and can be obtained from ftp://ftp.jax.org/ohlab/Index/. 

However, to reduce runtime, we provided environment-specific database that was created using microbes mostly found in a specific microbial niche. This can be retrieved from **ftp://ftp.jax.org/ohlab/GRiD_environ_specific_database/**. For instance, if you are analyzing stool samples, it is advisable to download and extract the stool database `ftp://ftp.jax.org/ohlab/GRiD_environ_specific_database/stool_microbes.tar.gz`.   

# INSTALLATION
The easiest way to install GRiD is through miniconda which resolves all required dependencies. 

1.    If you do not have anaconda or miniconda already installed, download miniconda https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh and run the install script. Reload .bashrc environment `source ~/.bashrc`
    
    -- Set up channels --
    
    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge
          
2.    Install GRiD

`conda install grid=1.2`

To avoid compatibility issues with dependencies, it is generally advisable to create a conda environment for any software to be installed using conda (for external reference, see https://conda.io/docs/user-guide/tasks/manage-environments.html#).

As an example, you can create an environment prior to GRiD installation as below

    conda create --name GRiD
    source activate GRiD
    conda install grid=1.2

**It is highly recommended to run the example test to ensure proper installation before running GRiD on your dataset. You do not need to have downloaded the GRiD database to run the test (see "Example test" below)**.


# USAGE

    grid single <options>       GRiD using a single genome
    grid multiplex <options>    GRiD high throughput
    grid -v                     Version
    grid -h                     Display this help message

    grid single <options>
    <options>
    -r      Reads directory (single end reads)
    -o      Output directory
    -g      Reference genome (fasta)
    -l      Path to file listing a subset of reads
            for analysis [default = analyze all samples in reads directory]
    -n INT  Number of threads [default 1]
    -h      Display this message

    grid multiplex <options>
    <options>
    -r         Reads directory (single end reads)
    -o         Output directory
    -d         GRiD database directory
    -c  FLOAT  Coverage cutoff (>= 0.2) [default 1]
    -p         Enable reassignment of ambiguous reads using Pathoscope2
    -t  INT    Theta prior for reads reassignment [default 0]. Requires the -p flag
    -l         Path to file listing a subset of reads
               for analysis [default = analyze all samples in reads directory]
    -m         merge output tables into a single matrix file
    -n  INT    Number of threads [default 1]
    -h         Display this message


**NOTE: Sample reads must be in single-end format**. If reads are only available in paired-end format, use either of the mate pairs, or concatenate both pairs into a single fastq file. In addition, reads must have the .fastq extension (and not .fq). Also, when specifying the output folder using the -o flag, either select a non-existing folder (in this case, GRiD would create the folder) or select an  existing but empty folder.  

In both 'single' and 'multiplex' modules, all samples present in the reads directory would be analyzed by default. However, analysis can be restricted to a subset of samples by using the -l flag and specifying a file that lists the subset of samples.    

For the 'multiplex' module, reads mapping to multiple genomes are reassigned using Pathoscope 2 when the -p flag is set. The degree to which reads are reassigned is set by the -t (theta prior) flag. The theta prior value represents the number of non-unique reads that are not subject to reassignment. Finally, when the coverage cutoff (-c flag) is set below 1, only genomes with fragmentation levels below 90 fragments/Mbp are analyzed (see Emiola & Oh, 2018 for more details). **Note that to use the 'multiplex' module, you must have downloaded either the environ_specific_database from `[ftp://ftp.jax.org/ohlab/GRiD_environ_specific_database/]` or comprehensive database `[ftp://ftp.jax.org/ohlab/Index/]`. However, you do not need the database to run the example test**.

# OUTPUT
`single module` - two output files are generated
- A plot (.pdf) showing coverage information across the genome 
- A table of results (.txt) displaying growth rate score (GRiD), 95% confidence interval, unrefined GRiD value, species heterogeneity, genome coverage, *dnaA/ori* ratio, and *ter/dif* ratio. Species heterogeneity is a metric estimating the degree to which closely related strains/species contribute to variance in growth predictions (range between 0 - 1 where 0 indicate no heterogeneity). In most bacteria genomes, *dnaA* is located in close proximity to the *ori* whereas replication typically terminates at/near *dif* sequence. Thus, the closer *dnaA/ori* and *ter/dif* ratios are to one, the more likely the accuracy of GRiD scores.  

`multiplex module` - two output files are generated per sample
- A table of results (.txt) displaying growth rate score (GRiD) for **ALL** genomes above the coverage cutoff, unrefined GRiD value, species heterogeneity, and genome coverage. If -m flag is set, all tables will be merged into a single matrix file called "merged_table.txt".
- A heatmap (.pdf), displaying growth rate score (GRiD) for the **70 most abundant species** above the coverage cutoff with hierachical clustering. 


**Some notes to keep in mind when using the 'multiplex' module to enhance performace**
-  If possible, use the smaller `environ_specific_database` GRiD database.
-  Subsample very large files to considerably reduce runtime. 
-  Whenever possible, run your analyses with and without the Pathoscope reassignment (-p) option and compare results. You can also fine-tune the reads reassignment parameters using -t flag.
-  Also, you may want to filter your results based on `Species heterogeneity`. Typically, growth estimates with 'Species heterogeneity' values < 0.3 are very reliable.

# Example test
The test sample contain reads from *Staphylococcus epidermids*, *Lactobacillus gasseri*, and *Campylobacter upsaliensis*, each with coverage of ~ 0.5. Download the GRiD folder and run the test as shown below. 

If you created a conda environment prior to installation using the example in the "Installation" section, then always activate the environment prior to running GRiD using the command `source activate GRiD`. Otherwise, skip this step. 

`wget https://github.com/ohlab/GRiD/archive/1.2.tar.gz`

`tar xvf 1.2.tar.gz`

`cd GRiD-1.2/test`

`grid single -r . -g S_epidermidis.LRKNS118.fna -o output_single`

`grid multiplex -r . -d . -p -c 0.2 -o output_multiplex -n 16`    (This command tells GRiD to reassign ambiguous reads, calculate genomes with coverage > 0.2, and use 16 threads) 

For each module, output files (a pdf and a text file) are generated in the test folder.


# Updating GRiD database 
The database can be updated with metagenomic bins or newly sequenced bacterial genomes by running the update_database script.
 

    update_database <options>
    <options>
    -d      GRiD database directory (required)
    -g      Bacterial genomes directory (required)
    -p      Prefix for new database (required)
    -l      Path to file listing specific genomes
            for inclusion in database [default = include all genomes in directory]
    -h      Display this message

NOTE: Genomes must be in fasta format and must have either .fasta, .fa or .fna extensions.


# DEPENDENCIES
- R (tested using v 3.4.1) 
    - Required R libraries - 
    dplyr,
    getopt,
    ggplot2,
    gsubfn,
    gplots
    
- bowtie2 (tested using v 2.3.1)
- seqtk (tested using v 1.0-r31)
- samtools (tested using v 1.5)
- bedtools (tested using v 2.26.0)
- bamtools (tested using v 1.0.2)
- blast (tested using v 2.6.0)
- Pathoscope2    
- parallel
- mosdepth 

