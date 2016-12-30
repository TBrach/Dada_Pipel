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