# Purpose of Dada Pipel

Dada Pipel is a wrapper that uses the Dada2 package to get from paired reads in fastQ files to the count table (phyloseq object)


# How to use it

NB: currently the wrappers are implemented for cases in which the FW and RV reads of each sample are in an individual folder. These "sample" folders are all in one folder. FW and RV reads can be distinguished by character patterns. 

## Step 1: Do a quality check of the sequence reads to decide on filtering parameters with (Dada_QualityCheck.R)

Dada_QualityCheck.R runs the Dada2_QualityCheck function. The Mean and Median quality scores at the different cycles over all reads are determined and plotted. The data is stored in newly generated folders: 

Specifically, the function Dada2_QualityCheck is a wrapper of:

- the ShortRead:::qa function (quality assessment on short reads)
    - in addition: 
        - checks if all reads were of the same length
        - calculates mean and median quality score for each cycle (nucleotide)



To run the function from the terminal, use:

-  Rscript Dada_QualityCheck.R
    - make sure to adjust the pathnames in Dada_QualityCheck.R


# Features