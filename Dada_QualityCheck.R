# ---- Source the functions ----

# ATTENTION: change pathToFunctions here if necessary#
pathToFunctions <- "/Users/jvb740/MarieCurie_Work/BackgroundKnowledge/16S_Learning/Dada_Pipel/Functions"
#pathToFunctions <- "/home/jvb740/Dada_Pipel"

source(file.path(pathToFunctions, "Dada_PlotFunctions.R"))
source(file.path(pathToFunctions, "Dada_WrapFunctions.R"))

# ----

# ---- Call the wrap function (Adjust INPUTS) ----
Dada2_QualityCheck(path = "/Users/jvb740/MarieCurie_Work/Galaxy_Mouse_Project/16S_Analysis/BGI_results/Clean_Data", #"/Users/jvb740/MarieCurie_Work/BackgroundKnowledge/16S_Learning/DanFunD",
           F_pattern = "1.fq.gz", 
           R_pattern = "2.fq.gz",
           path2 = "/Users/jvb740/MarieCurie_Work/Galaxy_Mouse_Project/16S_Analysis/BGI_results/Dada_Analysis")
# ----

# Then call on terminal Rscript Dada_QualityCheck.R