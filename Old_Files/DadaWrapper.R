# ---- Sourcing the plot and wrapper functions ----

# CHANGE pathToFunctions here:
pathToFunctions <- "/Users/jvb740/MarieCurie_Work/BackgroundKnowledge/16S_Learning/Dada_Pipel"
#pathToFunctions <- "/home/jvb740/Dada_Pipel"

source(file.path(pathToFunctions, "Dada_PlotFunctions.R"))
source(file.path(pathToFunctions, "Dada_WrapFunctions.R"))
# ----

# ---- Calling the wrap function (Adjust INPUTS) ----
Dada2_wrap(path = "/Users/jvb740/MarieCurie_Work/BackgroundKnowledge/16S_Learning/DanFunD",
           F_pattern = "1.fq.gz", 
           R_pattern = "2.fq.gz",
           path2 = NULL,
           trimLeft = c(10,10), # how many nucleotides will be clipped from the beginning of the FW and RV reads, respectively
           truncLen = c(220, 160), # how many nucleotides will be clipped from the end of the FW and RV reads, respectively
           maxEE = 1, # After truncation, reads with higher than maxEE "expected errors" will be discarded
           maxN = 0, # After truncation, sequences with more than maxN Ns will be discarded, NB: Dada2 requires no Ns!
           truncQ = 2,# Truncate reads at the first instance of a quality score less than or equal to truncQ 
           # NB: see: https://github.com/benjjneb/dada2/issues/140 
           # ? I did not understand this because I thought it should clash with dada2 not allowing sequences of variable length, but this i now supported:
           # https://github.com/benjjneb/dada2/issues/55 (at the end)
           NSAM.LEARN = 40, # the number of samples used to estimate the F and R error matrixes. When NULL all samples are used. Number of filtered reads should be 1 million
           err_F = NULL, # when error matrix given, the error matrix estimation is skipped
           err_R = NULL,
           minOverlap = 20, # minOverlap from the mergePairs command
           maxMismatch = 0, # maxMismatch from the mergePairs command
           F_QualityStats = NULL, # if given e.g. from Dada_QualityCheck the quality stats collection part is jumped over
           R_QualityStats = NULL,
           filtFs = NULL, # if given and FilteredFolder with files exist, Filtering can be jumped over
           filtRs = NULL)

# Then call on terminal: Rsript DadaWrapper.R


# path alternatives
# path = "/home/tdr438/16s/Danfund/Clean"
# path2 = "/home/jvb740/DanFunD_Test"