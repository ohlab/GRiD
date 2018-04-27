# GRiD
Growth Rate Index (GRiD) measures bacterial growth rate from reference genomes (including draft quality genomes) and metagenomic bins at ultra-low sequencing coverage (> 0.2x). 

GRiD algorithm consists of two modules;

1. < single > - which is applicable for growth analysis involving a single reference genome
2. < multiplex > - for the high-throughput growth analysis of all characterized bacteria in a sample. Prior knowledge of microbial composition is not required. To use this module, download the GRiD database, consisting of 32,819 representative bacteria genomes, from ftp://ftp.jax.org/ohlab/Index/   

# REQUIREMENTS
- R (tested using v 3.2.3) 
    - Required R libraries - 
    dplyr,
    getopt,
    ggplot2,
    gsubfn,
    glots
    
- bowtie2 (tested using v 2.3.1)
- seqtk (tested using v 1.0-r31)
- samtools (tested using v 1.5)
- gcc (tested using v 7.1.0)
- bedtools (tested using v 2.26.0)
- bamtools (tested using v 1.0.2)
- usearch (tested using v 8.0)
- python (tested using v 2.7.3)
- Pathoscope2    

Add the above dependencies to your PATH enivironment. For instance, to add Pathoscope2 to your PATH, run the following commands

`$ echo 'export PATH=/path/to/pathoscope2/folder:$PATH' >> ~/.bash_profile`

`$ source ~/.bash_profile`  

Additionally, if you do not have root access, R libraries can be installed locally. But first, create a directory to store all your R libraries and run the commands below
`$ echo 'R_LIBS_USER="/path/to/R_libraries/directory/"' >>  $HOME/.Renviron
-- Fire up an R session to install a library. e.g ggplot2
$ R
install.packages("ggplot2", lib="/path/to/R_libraries/directory/")
quit()`


# USAGE

    ./grid.sh -h                     Display this help message.
    ./grid.sh single <options>       GRiD using a single genome
    ./grid.sh multiplex <options>    GRiD high throughput

    ./grid.sh single <options>
    <options>
    -r      Reads directory (single end reads)
    -o      Output directory
    -g      Reference genome (fasta)
    -l      Path to file listing a subset of reads
            for analysis [default = analyze all samples in reads directory]
    -h      Display this message

    ./grid.sh multiplex <options>
    <options>
    -r         Reads directory (single end reads)
    -o         Output directory
    -d         GRiD database directory
    -c  FLOAT  Coverage cutoff (>= 0.2) [default 1]
    -p         Enable reassignment of ambiguous reads
    -t  INT    Theta prior for reads reassignment [default 0]. Requires the -p flag
    -l         Path to file listing a subset of reads
               for analysis [default = analyze all samples in reads directory]
    -m         merge output tables into a single matrix file
    -h         Display this message


NOTE: Samples must be in single-end format. If samples are only available in paired-end format, use either of the mate pairs, or concatenate both pairs into a single fastq file. In addition, reads must have the .fastq extension (and not .fq). Using either modules, the default is to analyze all samples present in the reads directory. However, analysis can be restricted to a subset of samples by using the -l flag and specifying a file that lists the subset of samples.    

For the 'multiplex' module, reads mapping to multiple genomes are reassigned using Pathoscope 2 when the -p flag is set. The degree to which reads are reassigned is set by the -t (theta prior) flag. The theta prior value represents the number of non-unique reads that are not subject to reassignment. Finally, when the coverage cutoff (-c flag) is set below 1, only genomes with fragmentation levels below 90 fragments/Mbp are analyzed (see xxx et al. for more details). 

# Updating GRiD database 
The database can be updated with metagenomic bins or newly sequenced bacterial genomes by running the update_database.sh script (requires bowtie2).
 

    ./update_database.sh <options>
    <options>
    -d      GRiD database directory
    -g      Bacterial genomes directory
    -n      Name for new database
    -l      Path to file listing specific genomes
            for inclusion in database [default = include all genomes in directory]
    -h      Display this message

NOTE: bacteria genomes must be in fasta format and must have either .fasta, .fa or .fna extensions.

# OUTPUT
`single module` - two output files are generated
- A plot (.pdf) showing coverage information across the genome 
- A table of results (.txt) displaying growth rate (GRiD), 95% confidence interval, unrefined GRiD value, species heterogeneity, genome coverage, dnaA/ori ratio, and ter/dif ratio. Species heterogeneity is a metric estimating the degree to which closely related strains/species contributes to variance in growth predictions (range between 0 - 1 where 0 indicate no heterogeneity). In most bacteria genomes, dnaA is located in close proximity to the ori whereas replication typically terminates at/near dif sequence. Thus, the closer dnaA/ori and ter/dif ratios are to one, the more likely the accuracy of GRiD scores.  

`multiplex module` - two output files are generated per sample
- A heatmap (.pdf), displaying growth rate (GRiD) from genomes above the coverage cutoff with hierachical clustering. 
- A table of results (.txt) displaying growth rate (GRiD) of genomes above the coverage cutoff, unrefined GRiD value, species heterogeneity, and genome coverage. If -m flag is set, all tables will be merged into a single matrix file called "merged_table.txt".
