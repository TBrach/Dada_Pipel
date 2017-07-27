# source seqtab data and function
datapath <- "/Users/jvb740/MarieCurie_Work/MouseProject/ResultsAndProtocols/ManiAging_Results/16S_Sequencing/2017-07-13_DK_age_ManiAging/Dada2_Analysis_trimLeft30"
load(file.path(datapath, "Dada_Data/DenoisedData.RData"))

functpath <- "/Users/jvb740/MarieCurie_Work/BackgroundKnowledge/16S_Learning/Dada_Pipel/Functions"

source(file.path(functpath, "Dada_TaxonomyFunctions.R"))

if(!exists("seqtab.nochim")){
        stop("no seqtab has been loaded")
}

assignTaxonomyaddSpecies(seqtab = seqtab.nochim, 
                         minBoot = 50,
                         allowMultiple = TRUE,
                         PathToRefs = "/Users/jvb740/MarieCurie_Work/BackgroundKnowledge/16S_Learning/AssignTaxonomy",
                         RefDataBase = "silva_nr_v128_train_set.fa.gz",
                         SpeciesDB = "silva_species_assignment_v128.fa.gz",
                         PathToSave = "/Users/jvb740/MarieCurie_Work/MouseProject/ResultsAndProtocols/ManiAging_Results/16S_Sequencing/2017-07-13_DK_age_ManiAging/Dada2_Analysis_trimLeft30/Dada_Taxonomy",
                         tryRC = FALSE)