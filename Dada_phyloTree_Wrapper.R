# source seqtab data and function
datapath <- "/Users/jvb740/MarieCurie_Work/Project_GalaxyMouse/DE_LiverDisease_GalaxyMouse_20180911_Reads/BGI_results/Dada_Analysis_Pooled/Dada_Data"
load(file.path(datapath, "DenoisedData.RData"))

functpath <- "/Users/jvb740/MarieCurie_Work/BackgroundKnowledge/16S_Learning/Dada_Pipel/Functions"

source(file.path(functpath, "Dada_TaxonomyFunctions.R"))

if(!exists("seqtab.nochim")){
        stop("no seqtab.nochim has been loaded")
}

construct_phylogenetic_tree(seqtab.nochim = seqtab.nochim, 
                            savepath = "/Users/jvb740/MarieCurie_Work/Project_GalaxyMouse/DE_LiverDisease_GalaxyMouse_20180911_Reads/BGI_results/Dada_Analysis_Pooled/Dada_phylogenetic_tree")

message("DONE!")