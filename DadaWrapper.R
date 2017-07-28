# ---- Sourcing the plot and wrapper functions ----
datapath <- "/emc/cbmr/data/MICROBIOME/raw/mouse/stool/2017-07-13_DK_age_ManiAging/"
load(file.path(datapath, "Dada_Data/QualityStats.RData"))

namesF <- names(filtFs)
namesR <- names(filtRs)
filtFs <- file.path("/emc/cbmr/data/MICROBIOME/raw/mouse/stool/2017-07-13_DK_age_ManiAging/Dada_FilteredFastqs", basename(filtFs))
filtRs <- file.path("/emc/cbmr/data/MICROBIOME/raw/mouse/stool/2017-07-13_DK_age_ManiAging/Dada_FilteredFastqs", basename(filtRs))
names(filtFs) <- namesF
names(filtRs) <- namesR


# CHANGE pathToFunctions here:
# pathToFunctions <- "/Users/jvb740/MarieCurie_Work/BackgroundKnowledge/16S_Learning/Dada_Pipel/Functions"
pathToFunctions <- "/home/jvb740/Dada_Pipel/Functions"

source(file.path(pathToFunctions, "Dada_PlotFunctions.R"))
source(file.path(pathToFunctions, "Dada_WrapFunctions.R"))
# ----

# ---- Calling the wrap function (Adjust INPUTS) ----
Dada2_wrap(path = "/emc/cbmr/data/MICROBIOME/raw/mouse/stool/2017-07-13_DK_age_ManiAging/Clean",
           F_pattern = "1.fq.gz", 
           R_pattern = "2.fq.gz",
           path2 = "/emc/cbmr/data/MICROBIOME/raw/mouse/stool/Pooled2",
           trimLeft = c(10,10), # how many nucleotides will be clipped from the beginning of the FW and RV reads, respectively
           truncLen = c(230, 165), # how many nucleotides will be clipped from the end of the FW and RV reads, respectively
           maxEE = 0.8, # After truncation, reads with higher than maxEE "expected errors" will be discarded
           maxN = 0, # After truncation, sequences with more than maxN Ns will be discarded, NB: Dada2 requires no Ns!
           truncQ = 2,# Truncate reads at the first instance of a quality score less than or equal to truncQ 
           # NB: see: https://github.com/benjjneb/dada2/issues/140 
           # ? I did not understand this because I thought it should clash with dada2 not allowing sequences of variable length, but this i now supported:
           # https://github.com/benjjneb/dada2/issues/55 (at the end)
           nreadsLearn = 1.2e+06, # the number of reads (distributed over far less unique reads) used to learn the error matrixes, i.e. nreads in dada2:::learnErrors
           err_F = err_F, # when error matrix given, the error matrix estimation is skipped
           err_R = err_R,
           minOverlap = 20, # minOverlap from the mergePairs command
           maxMismatch = 0, # maxMismatch from the mergePairs command
           F_QualityStats = F_QualityStats, # if given e.g. from Dada_QualityCheck the quality stats collection part is jumped over
           R_QualityStats = F_QualityStats,
           filtFs = filtFs, # if given and FilteredFolder with files exist, Filtering can be jumped over
           filtRs = filtRs,
           pool = TRUE)

# Then call on terminal: Rscript DadaWrapper.R


# path alternatives
# path = "/home/tdr438/16s/Danfund/Clean"
# path2 = "/home/jvb740/DanFunD_Test"
# "/Users/jvb740/MarieCurie_Work/MouseProject/ResultsAndProtocols/ManiAging_Results/16S_Sequencing/2017-07-13_DK_age_ManiAging_Pretest"