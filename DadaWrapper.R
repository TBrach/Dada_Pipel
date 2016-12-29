### Sourcing the plot and wrapper functions

# CHANGE pathToFunctions here:
pathToFunctions <- "/Users/jvb740/MarieCurie_Work/BackgroundKnowledge/16S_Learning/Dada_Pipel"
#pathToFunctions <- "/home/jvb740/Dada_Pipel"

source(file.path(pathToFunctions, "Dada_PlotFunctions.R"))
source(file.path(pathToFunctions, "Dada_WrapFunctions.R"))

### Calling the wrap function (Adjust INPUTS)
Dada2_wrap(path = "/Users/jvb740/MarieCurie_Work/BackgroundKnowledge/16S_Learning/DanFunD",
           F_pattern = "1.fq.gz", 
           R_pattern = "2.fq.gz",
           path2 = NULL,
           trimLeft = c(10,10),
           truncLen = c(220, 160),
           maxEE = 1,
           maxN = 0,
           truncQ = 2,
           NSAM.LEARN = 40,
           err_F = NULL,
           err_R = NULL,
           minOverlap = 20,
           maxMismatch = 0,
           F_QualityStats = NULL,
           R_QualityStats = NULL,
           filtFs = NULL,
           filtRs = NULL)

## Then call on terminal
# Rsript DadaWrapper.R


## path alternatives
# path = "/home/tdr438/16s/Danfund/Clean"
# path2 = "/home/jvb740/DanFunD_Test"