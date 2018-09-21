# source seqtab data and function
datapath <- "/Users/jvb740/MarieCurie_Work/NormalizationProjectNastya/Results/16S_Sequencing_Pre/Dada_Analysis"
load(file.path(datapath, "Dada_Data/DenoisedData.RData"))

functpath <- "/Users/jvb740/MarieCurie_Work/BackgroundKnowledge/16S_Learning/Dada_Pipel/Functions"

source(file.path(functpath, "Dada_TaxonomyFunctions.R"))

if(!exists("seqtab.nochim")){
        stop("no seqtab.nochim has been loaded")
}

construct_phylogenetic_tree(seqtab.nochim = seqtab.nochim, 
                            savepath = "/Users/jvb740/MarieCurie_Work/NormalizationProjectNastya/Results/16S_Sequencing_Pre/Dada_Analysis/Dada_phylogenetic_tree/")

message("DONE!")