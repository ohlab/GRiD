# GRiD
Growth Rate Index (GRiD) measures bacterial growth rate from draft genomes and metagenomic bins at ultra-low sequencing coverage (> 0.2x)

REQUIREMENTS
- R (tested using v 3.2.3) 
    - Required R libraries - 
    dplyr
    getopt
    ggplot2
    gsubfn
- bowtie2
- seqtk
- samtools
- gcc (tested using v 7.1.0)
- bedtools
    

USAGE
Download "GRiD_setup.qsub" and "GRiD.R" files 

Edit the file "GRiD_setup.qsub" and specify paths to your Reads directory, output directory, GRiD scripts directory, and finally, bowtie index file of genome. Submit the qsub file.

NOTE: The script assumes reads are paired-end and annotated as name_1.fastq and name_2.fastq, respectively. If delimiter separating reads isn't the underscore symbol, modify lines 19 and 23 accordingly. Likewise, if using 
