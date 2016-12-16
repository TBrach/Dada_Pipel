### Sourcing the plot and the wrapper function

# ATTENTION: change pathToFunctions here if necessary#
pathToFunctions <- "/Users/jvb740/MarieCurie_Work/BackgroundKnowledge/16S_Learning/Dada_Pipel"

source(file.path(pathToFunctions, "Dada_PlotFunctions.R"))
source(file.path(pathToFunctions, "Dada_WrapFunction.R"))

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
           NSAM.LEARN = 2,
           err_F = NULL,
           err_R = NULL,
           minOverlap = 20,
           maxMismatch = 0)

## Then call on terminal
# Rsript DadaWrapper.R