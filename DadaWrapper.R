source("/Users/jvb740/MarieCurie_Work/BackgroundKnowledge/16S_Learning/Dada_Pipel/161107_Dada_PlotFunctions.R")
source("/Users/jvb740/MarieCurie_Work/BackgroundKnowledge/16S_Learning/Dada_Pipel/161209_DadaWrapFunction.R")
# 
# Dada2_wrap(path = "/Users/jvb740/MarieCurie_Work/BackgroundKnowledge/16S_Learning/DanFunD",
#            F_pattern = "*F.fastq.gz", R_pattern <- "*R.fastq.gz")

#source("/home/jvb740/DanFunD_Dada/161107_Dada_PlotFunctions.R")
#source("/home/jvb740/DanFunD_Dada/161209_DadaWrapFunction.R")

Dada2_wrap(path = "/Users/jvb740/MarieCurie_Work/BackgroundKnowledge/16S_Learning/DanFunD",
           #path = "/Users/jvb740/MarieCurie_Work/BackgroundKnowledge/16S_Learning/DanFunD",
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
