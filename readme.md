# Purpose of Dada Pipel

Dada Pipel is a wrapper that uses the Dada2 package to get from paired reads in fastQ files to the count table (phyloseq object/ amplicon sequence variant (ASV or SV) table). 


# How to use it

NB: currently the wrappers are implemented for cases in which the FW and RV reads of each sample are in an individual folder. These "sample" folders are all collected in one "project" folder. FW and RV reads can be distinguished by character patterns. Important for dada2: *The forward and reverse fastqs contain reads in matched order*

## Step 1: Do a quality check of the sequence reads to decide on filtering parameters with (**Dada_QualityCheck.R**)

Dada_QualityCheck.R runs the **Dada2_QualityCheck** function. The Mean and Median quality scores at the different cycles over all reads are determined and plotted. The data is stored in newly generated folders Dada_Data and Dada_Plot: 

Specifically, the function Dada2_QualityCheck is a wrapper of:

- the ShortRead:::qa function (quality assessment on short reads)
    - in addition Dada2_QualityCheck: 
        - checks if all reads were of the same length
        - calculates mean and median quality score for each cycle (nucleotide)

To run the function from the terminal, use:

-  Rscript Dada_QualityCheck.R
    - make sure to adjust the pathnames in Dada_QualityCheck.R before calling it.


## Step 2: Run the actual Dada pipeline (**Dada_Wrapper.R**)

Dada_Wrapper.R runs the Dada2_wrap function. 

Specifically, the function Dada2_wrap is a wrapper of:

- Read quality check: the ShortRead:::qa function (quality assessment on short reads) (see Step 1 for additions)
- Filtering: dada2:::fastqPairedFilter
    - the filtered reads are saved in the folder Dada_FilteredFastqs using "SampleName"_F_Filtered.fastq.gz and "SampleName"_R_Filtered.fastq.gz
    - 



# Features

# Possible Improvements:

- Dada2 now offerst filterAndTrim which is a multithreaded convenience interface for the fastqFilter and fastqPairedFilter filtering functions
    - it offers multithread, If TRUE, input files are filtered in parallel via mclapply. So I think it is just a matter of speed no outcome changes. 