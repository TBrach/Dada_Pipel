# source seqtab data and function
datapath <- "/Users/jvb740/MarieCurie_Work/BackgroundKnowledge/16S_Learning/DanFunD_243_maxEE-1"
load(file.path(datapath, "Dada_Data/DenoisedData.RData"))

functpath <- ""

source(file.path(functpath, "Dada_TaxonomyFunctions.R"))

assignTaxonomyaddSpecies(seqtab, 
                         minBoot = 80,
                         allowMultiple = 3,
                         PathToRefs = NULL,
                         RefDataBase = "silva_nr_v123_train_set.fa.gz",
                         SpeciesDB = "silva_species_assignment_v123.fa.gz",
                         PathToSave = NULL)