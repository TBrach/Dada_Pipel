### Sourcing the plot and the wrapper function

# ATTENTION: change pathToFunctions here if necessary#
#pathToFunctions <- "/Users/jvb740/MarieCurie_Work/BackgroundKnowledge/16S_Learning/Dada_Pipel"
pathToFunctions <- "jvb740@porus01:/home/jvb740/Dada_Pipel"

source(file.path(pathToFunctions, "Dada_PlotFunctions.R"))
source(file.path(pathToFunctions, "Dada_WrapFunctions.R"))

### Calling the wrap function (Adjust INPUTS)
Dada2_wrap(path = "jvb740@porus01:/home/tdr438/16s/Danfund/Clean",
           F_pattern = "1.fq.gz", 
           R_pattern = "2.fq.gz",
           path2 = "jvb740@porus01:/home/jvb740/DanFunD_Test",
           trimLeft = c(10,10),
           truncLen = c(220, 160),
           maxEE = 1,
           maxN = 0,
           truncQ = 2,
           NSAM.LEARN = 45,
           err_F = NULL,
           err_R = NULL,
           minOverlap = 20,
           maxMismatch = 0)

## Then call on terminal
# Rsript DadaWrapper.R