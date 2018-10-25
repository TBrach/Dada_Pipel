# What to update or recent updates:

- 2018-10-25: added the randomize argument in again, that allows the samples to be randomized for error_matrix estimation (in learnErrorsAdj)
- 2018-10-25: added filteredQualityStats = TRUE option that simply allows to also make quality plots of the filtered Reads.
- the "pool" option of the dada() command has been added: pool = TRUE in Dada_Wrapper 
	- <https://benjjneb.github.io/dada2/pool.html>
	- NB: pooling for e.g. 70 samples makes the code really slow
	- it might however tackle dada2 assumed issue of alpha diversity by ignoring singletons.
- learnErrors uses now nbases instead of nreads, depending on the trimming you get indeed shorter reads in some cases than in others, but I think it is fine to not directly update to nbases in learnErrorsAdj
- you should also plot the quality of the filtered reads
- add the randomize option to learnErrorsAdj
	
# OF NOte

 - I removed the equal length requirment in the QualityCheck.R. So you should use the quality stats to record whether the reads from all samples had the same length!

# Requirements

- current code requires R version 3.4 up, and dada2 1.4
- Data requirement: FW and RV reads of each sample are in individual folders. These "sample" folders are all collected in one "project" folder. FW and RV reads can be distinguished by character patterns. Important for dada2: *The forward and reverse fastqs contain reads in matched order*

# Purpose of Dada_Pipel

Dada_Pipel offers wrapper function around the Dada2 package to use with Rscript. 
- Dada_QualityCheck.R: generates quality check plots of your fastQ data.
- DadaWrapper.R: runs the pipeline from the paired reads in your fastQ files to the count table (i.e. amplicon sequence variant (ASV or SV) table).
- Dada_WrapperAssignTaxonomyAddSpecies.R: a dada2 wrapper for taxonomy assignment, i.e. assignTaxonomy() and addSpecies()

# How to use it

## Step 1: **Dada_QualityCheck.R**

- runs the **Dada2_QualityCheck** function. The Mean and Median quality scores at the different cycles over all reads are determined and plotted. The data is stored in newly generated folders Dada_Data and Dada_Plot.

Specifically, the function Dada2_QualityCheck is a wrapper of:

- the ShortRead:::qa function (quality assessment on short reads)
- in addition: 
    - mean and median quality scores for each cycle (nucleotide) are calculated and plotted

To run the function from the terminal, use:

-  Rscript Dada_QualityCheck.R
    - make sure to adjust the pathnames in Dada_QualityCheck.R before calling it.
- nohup Rscript Dada_QualityCheck.R &


## Step 2: Run the actual Dada pipeline (**DadaWrapper.R**)

- on hulk:
	- nohup /usr/local/R-3.4.1/bin/Rscript DadaWrapper.R &
- on porus (here Rscript is directly linked to the new R.version 3.4.1)
	- nohup Rscript DadaWrapper.R &

### Steps of the Dada2_wrap function

- 1.) Construct character vectors to the FW and RV fastq files
    - Requires that your samples are in individual folders within the path folder. 
    - Generates F_fastq and R_fastq.
- 2.) Determine the quality scores and save the stats in Data folder
    - uses ShortRead:qa funtion on F_fastq[i] and R_fastq[i] to get a quality overview of each cycle. It checks whether the reads within a sample are all of the same length, but not whether the reads between different samples are of the same length.
    - generates F_QualityStats and R_QualityStats
- 3.) Generate and save some quality plots
    - uses F_QualityStats and R_QualityStats to generate some overview plots of the cycle quality.
- 4.) Filtering:
    - uses fastqPairedFilter (could use filterAndTrim wrapper but makes no difference) and saves filtered files in created Dada_FilteredFastqs folder
    - uses inputs trimLeft, truncLen, maxEE, maxN, truncQ
    - first uses truncLen, i.e. trims at the end (e.g. when your read is 250 bp and you set truncLen to 230 it will cleave of the last 20 bp), then it uses trimLeft, e.g. when you set it to 10 it will cut the first 10 bp, so leaving 220 from the intitially 250, then probably checking maxEE and so on. 
- 5.) Estimate err_F matrix and then err_R matrix
    - uses learnErrorsAdj which is basically learnErrors (NB: runs derepFastq() and dada() inside) from dada2 but saves dds and drps in the outcome in case all samples are used for error estimation. This allows in that case to not repeat derepFastq() and dada() in step 6. NB: It is possible to randomize the samples for error matrix estimation.
    - NB also: this is basically calling dada() without pool. NB also it is the number of unique reads that nreadsLearn refers to. This is on a per sample basis, so reads from different samples could be the same.
- 6.) dereplication, denoising, merging
    - dereplicate all filtered reads for a sample and then make use of the error matrix to denoise. 
        - so first derepFastq is used for dereplication resulting in a derep_class object drp_F, where length(drp_F$map) tells you how many filtered reads there were, and length(drp_F$uniques) how many unique reads there are. Also the average quality of each unique read at each cycle is saved in drp_F$quals.
        - then dada() is run using err_F or err_R, resulting in a dada-class object dd_F. NB: dd_F$denoised contains the denoised sequencing, sum(dd_F$denoised) <= sum(drp_F$uniques) showing that some unique reads are considered errors but are not assigned to a denoised sequence, so maybe these are the singletons that get removed. 
        - then mergePairs is used that uses drp_F, dd_F, drp_R, and dd_R as inputs. It basically first tries to pair each denoised FW sequence with each denoised RW sequence. So not each FW could be paired with several RV sequences. How you get to the actual abundances is via the two maps, you can see it in **20170712_UnderstandmergePairs.svg**  and 20170712_UnderstandmergePairsandMergers.R
    - 6a.) In case all samples have been used for error estimation, then dd_F, drp_F, dd_R, drp_R were saved as drp and dds in errorsRV and errorsFW. These are lists of the drp and dd of all samples and mergePairs allows list inputs.
        - mergePairs is used directly on the lists to merge the reads of the different samples: 
    - 6b.) In case err_F and err_R have been estimated only on a subset of samples, so errorsFW$dd and errorsRV$dd are NULL
        - derepFastq() and dada() (with the estimated (or given) err_F and err_R) and then mergePairs are run individually on the filtered reads of each sample.
- NB: on the chimera removal: 
    - Since dada2 tutorials are now all with removeBimeraDenovo(method = "conensus") after seqtab generation, we do this now and bimera removal is in step 8 in Dada2_wrap
    -  Dada2_wrap_BimFWRV is our old version, it uses isBimeraDenovo (probably similar to removeBimeraDenovo(method = 'pooled') separately on the FW and RV reads. Then mergers.nochim is created where only amplicons are kept for which both FW and RV reads were not identified as bimera. This usually identifies more amplicons as bimeras than does the new version in Dada2_wrap (step 8)
- 7.) Generate sequence table, just makeSequenceTable(mergers)
- 8.) remove bimeras: 
    - (NB: from tutorial: Fortunately, the accuracy of the sequences after denoising makes identifying chimeras simpler than it is when dealing with fuzzy OTUs: all sequences which can be exactly reconstructed as a bimera (two-parent chimera) from more abundant sequences.)
    -  If a majority of reads failed to pass the chimera check, you may need to revisit the removal of primers, as the ambiguous nucleotides in unremoved primers interfere with chimera identification.
    - See 20170117_UnderstandremoveBimeraDenovo




## Parameters to Consider


## Step 3: Dada_WrapperAssignTaxonomyAddSpecies.R

- assigns taxonomy to the found SV (from seqtab.nochim).
- assignTaxonomy() implements the RDP naive Bayesian classifier method described in Wang et al. 2007. 
    - assigns down to genus level based on reference database, currently silva_nr_v128_train_set.fa.gz or rdp_train_set_16.fa.gz
- assignSpecies (addSpecies) uses exact string matching vs a reference database, currently silva_species_assignment_v128.fa.gz or rdp_species_assignment_16.fa.gz, to assign species level taxonomy.

### Parameters to think about:

- minBoot for assignTaxonomy: sets the minimum bootstrapping support required to return a taxonomic classification. The original paper recommended a threshold of 50 for sequences of 250nts or less (as these are) but a threshold of 80 generally. The higher this value the more NA you get. 
    - I would stick to 50, I tried 80 and it was a bit more NA, I tried 25 and it was almost as 50, no real gain.
- allowMultiple for assignSpecies: usually you get NA when there is ambiguity, when allowMultiple = T, you get all the species that fitted, you can also use an integer allowMultiple = 3, then you get NA when sequence fitted to more than 3 species.
- tryRC: is tryRC in assignTaxonomy, should the reverse complement be checked if it fits better to the reference data base. Usually would only make sense if a lot of your assignments look really weird, like Eukaryota NA NA...
- NB: **We unfortunately usually get far lower percentages assigned to species level as they claim in their tutorial**: For human microbiome data and V4 sequencing, we typically see a bit under half of the abundant sequence variants (SVs) assigned unambiguously to species-level, and a bit under 2/3rds assigned ambiguously (ie. allowing multiple species assignment). This fraction drops off in less sampled environments and for rare SVs.
    - Question is of course, what are abundant sequence variants? The ones with a prevalence above 25\%?
- I was trying a longer trimLeft (maybe some remaining primer sequence prevents good assignment) but NB: filterAndTrim has by default trimLeft = 0.

Links:

- <https://benjjneb.github.io/dada2/assign.html>

## Step 4: construct_phylogenetic_tree

- 


# Possible Improvements:

- Dada2 now offers filterAndTrim which is a multithreaded convenience interface for the fastqFilter and fastqPairedFilter filtering functions
    - it offers multithread, If TRUE, input files are filtered in parallel via mclapply. So I think it is just a matter of speed no outcome changes. 
    
# Possible issues 

## alpha diversity confounded by total amplicons in a way that cannot be corrected for by rarefaction

- dada non pooled produces often data without singletons, rarefaction curves of this data plateau very early (e.g. at 10000 reads)
- in other words the data implies almost no influence of total amplicons above 15000 reads on richness
- there was however so far always a clear trend in the data that samples with higher total amplicons had higher richness
- I therefore think it is more likely to get rare SVs for samples with higher total amplicons, but still these amplicons get counts above 3 or so.
- to illustrate the trend: split your data samples into 2 groups high and low richness and test if total amplicons (taxa_sums) differs
- Solutions: 
    - test alpha diversity on residuals to lm richness ~ total amplicons
    - rarefy on filtered reads
    - use pooled
 