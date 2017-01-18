#######################################
### FUNCTION: assignTaxonomyaddSpecies
#######################################

# Function that simply calls assignTaxonomy and then addSpecies from dada2
## Input
# seqtab: your sequence table from dada2
# minBoot: minBoot for assignTaxonomy
# allowMultiple: default 3, allowMultiple for addSpecies
# PathToRefs: path to the folder with bot the RefDataBase and SpeciesDB in
# RefDataBase: the reference database for assignTaxonomy
# SpeciesDB: the reference database for addSpecies
# PathToSave: the folder into which taxa and taxa.species will be saved, if NULL working directory is used
## Output:
# taxa, taxa.species are saved as Taxonomy.RData


assignTaxonomyaddSpecies <- function(seqtab, 
                                     minBoot = 80, 
                                     allowMultiple = 3,
                                     PathToRefs = NULL,
                                     RefDataBase = "silva_nr_v123_train_set.fa.gz",
                                     SpeciesDB = "silva_species_assignment_v123.fa.gz",
                                     PathToSave = getwd()){
        
        ## Packages
        try(library(dada2), biocLite("dada2"))
        
        
        ## Reference Databases
        if(is.null(PathToRefs)){
                stop("PathToRefs must be given")
        }
        
        RefDB <- file.path(PathToRefs, RefDataBase)
        
        if(!file.exists(RefDB)){
                stop("could not find the reference database for assignTaxonomy")
        }
        
        SpecDB <- file.path(PathToRefs, SpeciesDB)
        
        if(!file.exists(SpecDB)){
                stop("could not find the reference database for addSpecies")
        }
        
        ## list the input
        InputSave <- list(seqtab = seqtab, minBoot = minBoot,
                                              allowMultiple = allowMultiple,
                                              PathToRefs = PathToRefs, 
                                              RefDataBase = RefDataBase,
                                              SpeciesDB = SpeciesDB,
                                              PathToSave = PathToSave)
        
        
        ## run assignTaxonomy and save
        ptm <- proc.time()
        taxa <- assignTaxonomy(seqtab, refFasta = RefDB, verbose = TRUE, minBoot = minBoot)
        TimeForassignTaxonomy <- (proc.time()-ptm)[3]
        
        InputSave[[length(InputSave) + 1]] <- TimeForassignTaxonomy
        
        save(taxa, InputSave, file = file.path(PathToSave, "Taxonomy.RData"))
        
        message("assignTaxonomy has been run")
        
        ## run addSpecies and save
        ptm <- proc.time()
        taxa.species <- addSpecies(taxa, refFasta = SpecDB, verbose = TRUE, allowMultiple = allowMultiple)
        TimeForaddSpecies <- (proc.time()-ptm)[3]
        
        InputSave[[length(InputSave) + 1]] <- TimeForaddSpecies
        
        save(taxa, taxa.species, InputSave, file = file.path(PathToSave, "Taxonomy.RData"))


}


#######################################
### FUNCTION: NAPerCForTaxLevel
#######################################

# Function to determine how many sequences could not be unambiguously (did not pass minBoot crtierion) assigned using a given minBoot value
## Input
# taxa: the output of the assignTaxonomy (dada2) command. A matrix with nrow = number of sequences and ncol usually 6, the taxonomic levels. When 
# the minBoot criterion was not fulfilled at a taxonomic level, assignTaxonomy assigns an NA. 
## Output
# NaPerC: data frame that states for each taxonomic level in taxa the number of NA and the percentage of NA. 

NAPerCForTaxLevel <- function(taxa){
        
        CountNA <- apply(taxa, 2, function(x){sum(is.na(x))})
        CountNA['Total'] <- nrow(taxa)
        NaPerC <- as.data.frame(CountNA)
        NaPerC$PerCNA <- 100*(NaPerC$CountNA/nrow(taxa))
        CountNonNA <- apply(taxa, 2, function(x){sum(!is.na(x))})
        CountNonNA['Total'] <- nrow(taxa)
        NaPerC$CountNonNA <- CountNonNA
        NaPerC$PerCNonNA <- 100*(NaPerC$CountNonNA/nrow(taxa))
        return(NaPerC)
        
}


#######################################
### FUNCTION: TaxLevelvsPenetrance
#######################################

# Function to determine how many sequences could be assigned to the different taxonomic levels compared to the number of smaples the amplicons are present
## Input
# taxa: the output of the assignTaxonomy (dada2) and maybe addSpecies command. nrow = number of sequences and ncol number taxonomic levels.
# seqtab: The abundance table output from Dada2_wrap
## Output
# TLvsPenetrance: list of 2 data frames that state for each penetrance level the assigned taxonomic levels, first data frame as numbers, second as percentage 
## requires: 
# dplyr!

TaxLevelvsPenetrance <- function(taxa, seqtab){
        
        TLvsPenetrance <- as.data.frame(cbind(colSums(seqtab != 0), apply(taxa, 2, function(x){!is.na(x)})))
        colnames(TLvsPenetrance)[1] <- "InNoSamples"
        TLvsPenetrance <- dplyr::group_by(TLvsPenetrance, InNoSamples)
        
        if("Species" %in% colnames(taxa)){
                
                TLvsPenetrance <- as.data.frame(dplyr::summarise(TLvsPenetrance, NoAmplicons = n(), Kingdom = sum(Kingdom), Phylum = sum(Phylum), Class = sum(Class),
                                                                 Order = sum(Order), Family = sum(Family), Genus = sum(Genus), Species = sum(Species)))
                TLPC <- dplyr::mutate(TLvsPenetrance, Kingdom = Kingdom/NoAmplicons, Phylum = 100*(Phylum/NoAmplicons), Class = 100*(Class/NoAmplicons), Order =
                                              100*(Order/NoAmplicons), Family = 100*(Family/NoAmplicons), Genus = 100*(Genus/NoAmplicons), Species = 100*(Species/NoAmplicons))
                
        } else {
                
                TLvsPenetrance <- as.data.frame(dplyr::summarise(TLvsPenetrance, NoAmplicons = n(), Kingdom = sum(Kingdom), Phylum = sum(Phylum), Class = sum(Class),
                                                                 Order = sum(Order), Family = sum(Family), Genus = sum(Genus)))
                TLPC <- dplyr::mutate(TLvsPenetrance, Kingdom = Kingdom/NoAmplicons, Phylum = 100*(Phylum/NoAmplicons), Class = 100*(Class/NoAmplicons), Order =
                                              100*(Order/NoAmplicons), Family = 100*(Family/NoAmplicons), Genus = 100*(Genus/NoAmplicons))
                                     
        }
        
        TLvsPenetrance <- list(Counts = TLvsPenetrance, PerC = TLPC)
        
        return(TLvsPenetrance)
        
}


#######################################
### FUNCTION: TaxLevelvsAbundance
#######################################

# Function to determine how many sequences could be assigned to the different taxonomic levels compared to the number of smaples the amplicons are present
## Input
# taxa: the output of the assignTaxonomy (dada2) and maybe addSpecies command. nrow = number of sequences and ncol number taxonomic levels.
# seqtab: The abundance table output from Dada2_wrap
## Output
# TLvsAb: list, first entry is the data.frame, second entry is a trellis object if Level was given

TaxLevelvsAbundance <- function(taxa, seqtab, Level = NULL){
        
        TLvsAbundance <- as.data.frame(cbind(colSums(seqtab), !is.na(taxa)))
        colnames(TLvsAbundance)[1] <- "Abundance"
        
        TLvsAb <- list(TLvsAbundance)
        
        if(!is.null(Level)){
                
                Tr <- ggplot(TLvsAbundance, aes_string(x = "Abundance", y = Level))
                Tr <- Tr + geom_point(col = "#E69F00", size = 3) +
                        scale_y_continuous(breaks = c(0,1), labels = c('No', "Yes")) +
                        theme_bw() +
                        theme(panel.grid.minor = element_blank())
                
                TLvsAb[[2]] <- Tr
        }

        return(TLvsAb)
        
}


#######################################
### FUNCTION: KeepAmplicons
#######################################

# Function to find the amplicons that are present in at least "Percentage" of samples
## Input
# taxa: the output of the assignTaxonomy (dada2) and maybe addSpecies command. nrow = number of sequences and ncol number taxonomic levels.
# seqtab: The abundance table output from Dada2_wrap
# percentage: only amplicons present in at least this percentage of samples will be kept
## Output
# taxseq: list with the shortened taxa and seqtab

KeepAmplicons <- function(taxa, seqtab, Percentage = 10){
        
        if(!identical(colnames(seqtab), rownames(taxa))){
                stop("taxa and seqtab do not fit togetehr")
        }
        
        
        PerCSampleValue <- ceiling((Percentage/100)*dim(seqtab)[1])
        KeepAmplis <- colnames(seqtab)[colSums(seqtab != 0) >= PerCSampleValue]
        seqtab.keep <- seqtab[,KeepAmplis]
        taxa.keep <- taxa[KeepAmplis,]
        taxseq <- list(taxa.keep, seqtab.keep)
        return(taxseq)
}






#######################################
### gm_own: calculate geometric mean
#######################################
# see commented below, this function comes from 
# <http://stackoverflow.com/questions/2602583/geometric-mean-is-there-a-built-in>

## Input:
# x numeric vector
# na.rm: if FALSE you get NA as soon as an NA is in your data, if TRUE the NA get basically treated as 0 (but NOTE these zeros always count also when zero.count = FALSE)
# zero.count: This is IMPORTANT, if TRUE 0s count, so if x contains a 0 the GM will be lower than when zero.count = FALSE, in 
# which case the gm is calculated for the same vector in which all 0 have been removed.
gm_own = function(x, na.rm=FALSE, zero.count = TRUE){
        if(any(x < 0, na.rm = TRUE)){
                return(NaN)
        }
        if(zero.count){
                exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
        } else {
                exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x[x > 0]))
        }
}

# Explanation of original gm function:
# Note zero.propagate just makes sure you get a 0 if there is any zero in your vector
# if there is no zero exp(mean(log(x), na.rm = TRUE)) == exp(sum(log(x[x > 0]), na.rm=TRUE) / length(x))
# I think this zero.propagate is silly, I therefore decided to make my own function, gm_own below

gm = function(x, na.rm=TRUE, zero.propagate = FALSE){
        if(any(x < 0, na.rm = TRUE)){
                return(NaN)
        }
        if(zero.propagate){
                if(any(x == 0, na.rm = TRUE)){
                        return(0)
                }
                exp(mean(log(x), na.rm = na.rm))
        } else {
                exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
        }
}






normalisealaDeSeq <- function(physeq) {
        
        
        
        
}