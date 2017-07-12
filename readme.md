# Purpose of Dada Pipel

Dada Pipel is a wrapper that uses the Dada2 package to get from paired reads in fastQ files to the count table (phyloseq object/ amplicon sequence variant (ASV or SV) table). 


# How to use it

NB: currently the wrappers are implemented for cases in which the FW and RV reads of each sample are in individual folders. These "sample" folders are all collected in one "project" folder. FW and RV reads can be distinguished by character patterns. Important for dada2: *The forward and reverse fastqs contain reads in matched order*

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

# Steps of the Dada2_wrap function

- 1.) Construct character vectors to the FW and RV fastq files
- 2.) Determine the quality scores and save the stats in Data folder
- 3.) Generate and save some quality plots
- 4.) Filtering:
    - uses fastqPairedFilter and saves filtered files in created Dada_FilteredFastqs folder
- 5.) Estimate err_F matrix and then err_R matrix
    - uses learnErrorsAdj which is basically learnErrors (NB: runs derepFastq() and dada() inside) from dada2 but saves dds and drps in the outcome in case all samples are used for error estimation. This allows in that case to not repeat derepFastq() and dada() in step 6. NB: to allow this the randomize option of learnErrors was removed. 
- 6.) dereplication, denoising, merging, (bimera identification old version)
- 6a.) In case all samples had been used for error estimation, i.e. dd_F, drp_F, dd_R, drp_R were saved as drp and dds in errorsRV and errorsFW. These are lists of the drp and dd of all samples.
    - The number of filtered reads, unique reads (drp), and denoised reads are saved for each sample using sapply.
    - mergePairs is used directly on the lists to merge the reads of the different samples: See: **UnderstandMergers.R** for understanding the resulting abundances.
- 6b.) err_F and err_R have been estimated only on a subset of samples, so errorsFW$dd and errorsRV$dd are NULL
    - derepFastq() and dada() (with the estimated (or given) err_F and err_R) and then mergePairs are run individually on the filtered reads of each sample.
- NB: on the chimera removal: our old version was based on isBimeraDenovo (probably similar to removeBimeraDenovo(method = 'pooled') and was based on bimera identification separately in the FW and RV reads, and only mergers were kept were both FW and RV were not identified as bimera. Since dada2 tutorials are now all with removeBimeraDenovo(method = "conensus") after seqtab generation, we do this now, and bimera removal is in step 8 (the new version removes normally some viewer bimeras).
- 7.) Generate sequence table, just makeSequenceTable(mergers)
- 8.) remove bimeras:
    - 

# Features

# Possible Improvements:

- Dada2 now offers filterAndTrim which is a multithreaded convenience interface for the fastqFilter and fastqPairedFilter filtering functions
    - it offers multithread, If TRUE, input files are filtered in parallel via mclapply. So I think it is just a matter of speed no outcome changes. 