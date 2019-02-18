# source seqtab data and function
# datapath <- "/Users/jvb740/MarieCurie_Work/Project_Normalization/DK_Healthy_Normalization_20181023_Reads/Dada_Analysis_filteredReads42500/Dada_Data/"
# load(file.path(datapath, "DenoisedData.RData"))

datapath <- "/Users/jvb740/MarieCurie_Work/Project_DanFundIgA/Data_Analysis_DanFundIgA/PhyseqsAndOtherData"

seqtab.nochim <- readRDS(file = file.path(datapath, "seqtab.nochim.merge_DF_82.rds"))

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
                         PathToSave = "/Users/jvb740/MarieCurie_Work/Project_DanFundIgA/20180928_DK_Age_DanFundIgA/Dada_Analysis_20190205/Dada_Taxonomy_Combined82/Silva_v132_minBoot50/",
                         tryRC = FALSE)
