# ---- Source the functions ----

# ATTENTION: change pathToFunctions here if necessary#
pathToFunctions <- "/Users/jvb740/MarieCurie_Work/BackgroundKnowledge/16S_Learning/Dada_Pipel"
#pathToFunctions <- "/home/jvb740/Dada_Pipel"

source(file.path(pathToFunctions, "Dada_PlotFunctions.R"))
source(file.path(pathToFunctions, "Dada_WrapFunctions.R"))

# ----

# ---- Call the wrap function (Adjust INPUTS) ----
Dada2_QualityCheck(path = "/Users/jvb740/MarieCurie_Work/BackgroundKnowledge/16S_Learning/DanFunD",
           F_pattern = "1.fq.gz", 
           R_pattern = "2.fq.gz",
           path2 = NULL)
# ----

# Then call on terminal Rscript Dada_QualityCheck.R