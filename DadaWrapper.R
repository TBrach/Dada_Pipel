# # ---- Sourcing the plot and wrapper functions ----
# # datapath <- "/emc/cbmr/data/MICROBIOME/raw/mouse/stool/2017-07-13_DK_age_ManiAging/Dada2_Analysis"
# datapath <- "/Users/jvb740/MarieCurie_Work/DanFunDProject/Pretest/Pretest_Dada_Analysis"
# load(file.path(datapath, "Dada_Data/QualityStats.RData"))
# 
# namesF <- names(filtFs)
# namesR <- names(filtRs)
# 
# #filtFs <- file.path("/emc/cbmr/data/MICROBIOME/raw/mouse/stool/2017-07-13_DK_age_ManiAging/Dada2_Analysis_filt_rare_5000/Dada_FilteredReads_rarefied", basename(filtFs))
# #filtRs <- file.path("/emc/cbmr/data/MICROBIOME/raw/mouse/stool/2017-07-13_DK_age_ManiAging/Dada2_Analysis_filt_rare_5000/Dada_FilteredReads_rarefied", basename(filtRs))
# filtFs <- file.path("/Users/jvb740/MarieCurie_Work/DanFunDProject/Pretest/Pretest_Dada_Analysis_FilterRarefied/Dada_FilteredReads_rarefied", basename(filtFs))
# filtRs <- file.path("/Users/jvb740/MarieCurie_Work/DanFunDProject/Pretest/Pretest_Dada_Analysis_FilterRarefied/Dada_FilteredReads_rarefied", basename(filtRs))
# names(filtFs) <- namesF
# names(filtRs) <- namesR


# CHANGE pathToFunctions here:
pathToFunctions <- "/Users/jvb740/MarieCurie_Work/BackgroundKnowledge/16S_Learning/Dada_Pipel/Functions"
# pathToFunctions <- "/home/jvb740/Dada_Pipel/Functions"

source(file.path(pathToFunctions, "Dada_PlotFunctions.R"))
source(file.path(pathToFunctions, "Dada_WrapFunctions.R"))
# ----

# ---- Calling the wrap function (Adjust INPUTS) ----
Dada2_wrap(path = "/Users/jvb740/MarieCurie_Work/Galaxy_Mouse_Project/16S_Analysis/BGI_results/Clean_Data",
           F_pattern = "1.fq.gz", 
           R_pattern = "2.fq.gz",
           path2 = "/Users/jvb740/MarieCurie_Work/Galaxy_Mouse_Project/16S_Analysis/BGI_results/Dada_Analysis", #"/emc/cbmr/data/MICROBIOME/raw/mouse/stool/Pooled2"
           trimLeft = c(10,10), # how many nucleotides will be clipped from the beginning of the FW and RV reads, respectively
           truncLen = c(230, 190), # how many nucleotides will be clipped from the end of the FW and RV reads, respectively
           maxEE = 1, # After truncation, reads with higher than maxEE "expected errors" will be discarded
           maxN = 0, # After truncation, sequences with more than maxN Ns will be discarded, NB: Dada2 requires no Ns!
           truncQ = 2,# Truncate reads at the first instance of a quality score less than or equal to truncQ 
           # NB: see: https://github.com/benjjneb/dada2/issues/140 
           # ? I did not understand this because I thought it should clash with dada2 not allowing sequences of variable length, but this i now supported:
           # https://github.com/benjjneb/dada2/issues/55 (at the end)
           nreadsLearn = 1.2e+06, # the number of reads (distributed over far less unique reads) used to learn the error matrixes, i.e. nreads in dada2:::learnErrors
           err_F = NULL, # when error matrix given, the error matrix estimation is skipped
           err_R = NULL,
           minOverlap = 30, # minOverlap from the mergePairs command
           maxMismatch = 0, # maxMismatch from the mergePairs command
           F_QualityStats = NULL, # if given e.g. from Dada_QualityCheck the quality stats collection part is jumped over
           R_QualityStats = NULL,
           filtFs = NULL, # if given and FilteredFolder with files exist, Filtering can be jumped over
           filtRs = NULL,
           pool = FALSE)

# Then call on terminal: Rscript DadaWrapper.R

# path = "/emc/cbmr/data/MICROBIOME/raw/mouse/stool/2017-07-13_DK_age_ManiAging/Clean"
# path2 = "/emc/cbmr/data/MICROBIOME/raw/mouse/stool/2017-07-13_DK_age_ManiAging/Dada2_Analysis_filt_rare_5000",
# path = "/Users/jvb740/MarieCurie_Work/MouseProject/ResultsAndProtocols/ManiAging_Results/16S_Sequencing/2017-07-13_DK_age_ManiAging/Clean"
# path2 = "/Users/jvb740/MarieCurie_Work/MouseProject/ResultsAndProtocols/ManiAging_Results/16S_Sequencing/2017-07-13_DK_age_ManiAging/Dada2_Analysis_filt_rare_10000",



# path alternatives
# path = "/home/tdr438/16s/Danfund/Clean"
# path2 = "/home/jvb740/DanFunD_Test"
# "/Users/jvb740/MarieCurie_Work/MouseProject/ResultsAndProtocols/ManiAging_Results/16S_Sequencing/2017-07-13_DK_age_ManiAging_Pretest"