# source seqtab data and function
datapath <- "/Users/jvb740/MarieCurie_Work/DanFunDProject/Pretest/Pretest_DadaAnalysis_FilterRarified_Pooled/"
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
                         RefDataBase = "silva_nr_v128_train_set.fa.gz",
                         SpeciesDB = "silva_species_assignment_v128.fa.gz",
                         PathToSave = "/Users/jvb740/MarieCurie_Work/DanFunDProject/Pretest/Pretest_DadaAnalysis_FilterRarified_Pooled/Dada_Taxonomy/Silva_v128",
                         tryRC = FALSE)