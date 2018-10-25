# source seqtab data and function
datapath <- "/Users/jvb740/MarieCurie_Work/Project_GalaxyMouse/DE_LiverDisease_GalaxyMouse_20180911_Reads/BGI_results/Dada_Analysis_Pooled/"
load(file.path(datapath, "Dada_Data/DenoisedData.RData"))

functpath <- "/Users/jvb740/MarieCurie_Work/BackgroundKnowledge/16S_Learning/Dada_Pipel/Functions"

source(file.path(functpath, "Dada_TaxonomyFunctions.R"))

if(!exists("seqtab.nochim")){
        stop("no seqtab.nochim has been loaded")
}

assignTaxonomyaddSpecies(seqtab = seqtab.nochim, 
                         minBoot = 50,
                         allowMultiple = TRUE,
                         PathToRefs = "/Users/jvb740/MarieCurie_Work/BackgroundKnowledge/16S_Learning/AssignTaxonomy",
                         RefDataBase = "silva_nr_v132_train_set.fa.gz",
                         SpeciesDB = "silva_species_assignment_v132.fa.gz",
                         PathToSave = "/Users/jvb740/MarieCurie_Work/Project_GalaxyMouse/DE_LiverDisease_GalaxyMouse_20180911_Reads/BGI_results/Dada_Analysis_Pooled/Dada_Taxonomy/Silva_v132_minBoot50/",
                         tryRC = FALSE)
