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
### FUNCTION: TaxLevelvsprevalence
#######################################

# Function to determine how many sequences could be assigned to the different taxonomic levels compared to the number of smaples the amplicons are present
## Input
# taxa: the output of the assignTaxonomy (dada2) and maybe addSpecies command. nrow = number of sequences and ncol number taxonomic levels.
# seqtab: The abundance table output from Dada2_wrap
## Output
# TLvsprevalence: list of 2 data frames that state for each prevalence level the assigned taxonomic levels, first data frame as numbers, second as percentage 
## requires: 
# dplyr!

TaxLevelvsprevalence <- function(taxa, seqtab){
        
        TLvsprevalence <- as.data.frame(cbind(colSums(seqtab != 0), apply(taxa, 2, function(x){!is.na(x)})))
        colnames(TLvsprevalence)[1] <- "InNoSamples"
        TLvsprevalence <- dplyr::group_by(TLvsprevalence, InNoSamples)
        
        if("Species" %in% colnames(taxa)){
                
                TLvsprevalence <- as.data.frame(dplyr::summarise(TLvsprevalence, NoAmplicons = n(), Kingdom = sum(Kingdom), Phylum = sum(Phylum), Class = sum(Class),
                                                                 Order = sum(Order), Family = sum(Family), Genus = sum(Genus), Species = sum(Species)))
                TLPC <- dplyr::mutate(TLvsprevalence, Kingdom = Kingdom/NoAmplicons, Phylum = 100*(Phylum/NoAmplicons), Class = 100*(Class/NoAmplicons), Order =
                                              100*(Order/NoAmplicons), Family = 100*(Family/NoAmplicons), Genus = 100*(Genus/NoAmplicons), Species = 100*(Species/NoAmplicons))
                
        } else {
                
                TLvsprevalence <- as.data.frame(dplyr::summarise(TLvsprevalence, NoAmplicons = n(), Kingdom = sum(Kingdom), Phylum = sum(Phylum), Class = sum(Class),
                                                                 Order = sum(Order), Family = sum(Family), Genus = sum(Genus)))
                TLPC <- dplyr::mutate(TLvsprevalence, Kingdom = Kingdom/NoAmplicons, Phylum = 100*(Phylum/NoAmplicons), Class = 100*(Class/NoAmplicons), Order =
                                              100*(Order/NoAmplicons), Family = 100*(Family/NoAmplicons), Genus = 100*(Genus/NoAmplicons))
                                     
        }
        
        TLvsprevalence <- list(Counts = TLvsprevalence, PerC = TLPC)
        
        return(TLvsprevalence)
        
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
# zeros.count: This is IMPORTANT, if TRUE 0s count, so if x contains a 0 the GM will be lower than when zero.count = FALSE, in 
# which case the gm is calculated for the same vector in which all 0 have been removed.
## Output:
# the geometric mean of x
gm_own = function(x, na.rm=FALSE, zeros.count = TRUE){
        if(any(x < 0, na.rm = TRUE)){
                return(NaN)
        }
        if(zeros.count){
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



#######################################
### filttaxa_by_prevalence
#######################################
# filters taxa by prevalence and provides plots that illustrate how the sample_sums
# were changed by the filtering

## Input:
# physeq
# prevalence: in %, only take stay that are present in more than prevalence % of the samples unless MaxCountCheck is used
# MaxCountCheck, MaxCount: if max count is used a taxa can remain even if not fulfilling the prevalence
# criterion if its max count in one of the samples is bigger than MaxCount 
# the idea is that some samples might depend very much on a rare taxa

## Output:
# list with the filtered physeq, the sparsity measures, and 3 different plots (3-5)


filttaxa_by_prevalence <- function(physeq, prevalence = 30, MaxCountCheck = FALSE, MaxCount = 0.2*(min(sample_sums(physeq)))) {
        
        if (MaxCountCheck) {
                physeq_f <- filter_taxa(physeq, function(x){(sum(x != 0) > (prevalence/100)*length(x)) || (max(x) > MaxCount)}, prune = TRUE)
        } else {
                physeq_f <- filter_taxa(physeq, function(x){(sum(x != 0) > (prevalence/100)*length(x))}, prune = TRUE)
        }
        Sparsity <- 100*(sum(otu_table(physeq) == 0)/(ntaxa(physeq)*nsamples(physeq)))
        Sparsity_f <- 100*(sum(otu_table(physeq_f) == 0)/(ntaxa(physeq_f)*nsamples(physeq_f)))
        PCCountsRemoved <- 100*(1-(sample_sums(physeq_f)/sample_sums(physeq)))
        DF <- data.frame(Sample = names(PCCountsRemoved), PCRemoved = PCCountsRemoved)
        DF$Sample <- as.factor(DF$Sample)
        DF <- dplyr::arrange(DF, PCRemoved)
        # relevel so samples are shown from min NoReads to max NoReads
        LevelsWant <- as.character(DF$Sample) 
        for (i in 1:length(LevelsWant)) {
                DF$Sample <- relevel(DF$Sample, ref = LevelsWant[i])
        }
        
        DF2 <- data.frame(before = sample_sums(physeq), after = sample_sums(physeq_f))
        DF2$Sample <- rownames(DF2)
        DF2$Sample <- as.factor(DF2$Sample)
        DF2 <- dplyr::arrange(DF2, before)
        LevelsWant <- as.character(DF2$Sample) 
        for (i in 1:length(LevelsWant)) {
                DF2$Sample <- relevel(DF2$Sample, ref = LevelsWant[i])
        }
        DF2 <- tidyr::gather(DF2, key = physeq, value = count, -Sample)
        
        Tr <- ggplot(DF, aes(x = Sample, y = PCRemoved))
        Tr <- Tr + geom_point(col = "#E69F00", size = 2) +
                ylab("% of amplicon reads removed") +
                xlab("") +
                theme_bw() +
                theme(panel.grid.minor = element_blank(),
                      panel.grid.major.y = element_blank(),
                      panel.grid.major.x = element_line(color = "#999999", size = .15),
                      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6))
        
        Tr2 <- ggplot(data = DF, aes(x = PCRemoved))
        Tr2 <- Tr2 + geom_histogram(binwidth = 2, col = "black", fill = "#E69F00") +
                geom_rug() +
                ylab("No Samples") + 
                xlab("% of amplicon reads removed") +
                theme_bw() + 
                theme(panel.grid.minor = element_blank(),
                      panel.grid.major.y = element_blank(),
                      panel.grid.major.x = element_line(color = "#999999", size = .15))
        
        cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
        
        Tr3 <- ggplot(data = DF2, aes(x = Sample, y = count, color = physeq))
        Tr3 <- Tr3 + geom_point(size = 2) +
                ylab("sample_sums") +
                xlab("") +
                scale_color_manual(values = cbPalette[c(4,2)]) +
                theme_bw() +
                theme(panel.grid.minor = element_blank(),
                      panel.grid.major.y = element_blank(),
                      panel.grid.major.x = element_line(color = "#999999", size = .15),
                      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
                      legend.title = element_blank())
        
        
        list(FilteredPhyseq = physeq_f, Sparsity = data.frame(Before = Sparsity, After = Sparsity_f), PlotReadsRemoved = Tr,
             HistReadsRemoved = Tr2, PlotSampleSumsBeforeAfter = Tr3)
        
}

#######################################
### adjust_LS
#######################################
# own implementation of DESeq2 library size adjustment with little tweaks

## Input:
# physeq
# zeros.count: zeros.count of gm_own, if TRUE 0 will be considered when calculating the geometric mean
# if FALSE not and thus the geometric means will be bigger (see gm_own)
# percentile: the percentile to select the size factor SF of a sample based on its count ratios to the reference sample. In
# DESeq percentile = 50, i.e. the median is used. NOTE when ignore.zero.ratios = FALSE, you might get SF = 0 in case that more 
# than 50% of the taxa of in a sample are 0. Function will through a warning then. s
# ignore.zero.ratios: if TRUE, the SF of each sample (i.e. the percentile) will be calculated based on the pool of non-zero ratios of the sample
# plots: if TRUE SFs and plots will be given in an output list, otherwise just the adjusted physeq is returned


adjust_LS <- function(physeq, zeros.count = FALSE, percentile = 50, ignore.zero.ratios = TRUE, plots = FALSE)  {
        
        # ---- Step 1: Use the geometric means for each taxa over all samples to calculate a "representative" Reference Sample ------
        if(taxa_are_rows(physeq)){
                GM <- apply(otu_table(physeq), 1, gm_own, zeros.count = zeros.count)   
        } else {
                GM <- apply(otu_table(physeq), 2, gm_own, zeros.count = zeros.count) 
        }
        
        # ---- Step 2: Calculate Size factors --------
        # -- 2a: Calculate Ratio Matrix by dividing counts of each sample by Reference Sample
        
        if(taxa_are_rows(physeq)){
                RatioMatrix <- sweep(otu_table(physeq), 1, GM, "/")
        } else {
                RatioMatrix <- sweep(otu_table(physeq), 2, GM, "/") 
        }
        
        # -- 2b: Insert: Calculate 0 Percentage for each sample
        
        if(taxa_are_rows(physeq)){
                zeroPC <- 100*(colSums(otu_table(physeq) == 0)/ntaxa(physeq))
        } else {
                zeroPC <- 100*(rowSums(otu_table(physeq) == 0)/ntaxa(physeq))
        }
        
        
        # -- 2c: calculate SFs for each sample depending on given percentile and whether to ignore.zero.ratios
        # if ignore.zero.ratios = TRUE, the percentile refers only to the non zero ratios in each sample
        # NB: DESEQ uses ignore.zero.ratios = TRUE
        if(ignore.zero.ratios){
                if (taxa_are_rows(physeq)) {
                        SFs <- apply(RatioMatrix, 2, function(x){x <- x[x>0]; quantile(x, probs = percentile/100, na.rm = T)})
                } else {
                        SFs <- apply(RatioMatrix, 1, function(x){x <- x[x>0]; quantile(x, probs = percentile/100, na.rm = T)})
                }  
        } else {
                if (taxa_are_rows(physeq)) {
                        SFs <- apply(RatioMatrix, 2, quantile, probs = percentile/100, na.rm = T)
                } else {
                        SFs <- apply(RatioMatrix, 1, quantile, probs = percentile/100, na.rm = T) 
                }   
        }
        
        if(min(SFs) == 0) { warning("in at least one sample the Size Factor was 0!") }
        
        
        # --- 3: calculate the new counts and put into a physeq object
        
        if(taxa_are_rows(physeq)){
                if (!identical(colnames(otu_table(physeq)), names(SFs))) {stop("names SFs do not fit to physeq")}
                Mat <- sweep(otu_table(physeq),2,SFs, "/")
                phynew <- phyloseq(otu_table(Mat, taxa_are_rows = TRUE), sample_data(physeq), tax_table(physeq))
        } else {
                if (!identical(rownames(otu_table(physeq)), names(SFs))) {stop("names SFs do not fit to physeq")}
                Mat <- sweep(otu_table(physeq),1,SFs, "/") 
                phynew <- phyloseq(otu_table(Mat, taxa_are_rows = FALSE), sample_data(physeq), tax_table(physeq))
        }
        
        if (plots){
                
                # ---- step 4: Generate Plots that illustrate the process
                # -- 4a: calculate for each sample a histogram of the ratios used to determine the SF
                
                cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
                
                # Define the histos function
                histos <- function(x, SF, zPC, SampleName) {
                        x <- data.frame(x = x)
                        Tr <- ggplot(x, aes(x = x))
                        Tr <- Tr + geom_histogram(binwidth = 0.3, col = "black", fill = "#E69F00") +
                                geom_rug() +
                                geom_vline(xintercept = SF, col = "#009E73") +
                                ylab("Frequency") + 
                                xlab("Count/(Count of Ref Sample)") +
                                ggtitle(paste("Sample: ", SampleName, "; zeroPC: ", round(zPC, 2), "; size factor: ", round(SF,2), sep = "")) +
                                theme_bw() + 
                                theme(panel.grid.minor = element_blank(),
                                      panel.grid.major.y = element_blank(),
                                      panel.grid.major.x = element_line(color = "#999999", size = .15))
                }
                
                
                if(taxa_are_rows(physeq) && (!all.equal(colnames(RatioMatrix), names(SFs), names(zeroPC)))) {warning("names messed up!")}
                if(!taxa_are_rows(physeq) && (!all.equal(rownames(RatioMatrix), names(SFs), names(zeroPC)))) {warning("names messed up!")}
                
                if(ignore.zero.ratios){
                        if(taxa_are_rows(physeq)){
                                HistoList <- lapply(1:nsamples(physeq), function(i){histos(x = RatioMatrix[,i][RatioMatrix[,i] > 0], SF = SFs[i], zPC = zeroPC[i], SampleName = names(SFs)[i])})
                        } else {
                                HistoList <- lapply(1:nsamples(physeq), function(i){histos(x = RatioMatrix[i,][RatioMatrix[i,] > 0], SF = SFs[i], zPC = zeroPC[i], SampleName = names(SFs)[i])})
                        }  
                } else {
                        if(taxa_are_rows(physeq)){
                                HistoList <- lapply(1:nsamples(physeq), function(i){histos(x = RatioMatrix[,i], SF = SFs[i], zPC = zeroPC[i], SampleName = names(SFs)[i])})
                        } else {
                                HistoList <- lapply(1:nsamples(physeq), function(i){histos(x = RatioMatrix[i,], SF = SFs[i], zPC = zeroPC[i], SampleName = names(SFs)[i])}) 
                        }   
                }
                
                # -- 3b: compare calculated SFs to library sizes of the samples
                if(!identical(names(SFs), names(sample_sums(physeq)))){warning("Probably some Mix Up CHECK!")}
                comp <- data.frame(Sample = names(SFs), TotAmps = sample_sums(physeq), SFs = SFs)
                comp$SFsNormed <- comp$SFs/median(comp$SFs)
                comp$TotAmpsNormed <- comp$TotAmps/mean(comp$TotAmps)
                comp <- dplyr::arrange(comp, desc(TotAmpsNormed))
                comp$Sample <- as.factor(comp$Sample)
                LevelsWant <- as.character(comp$Sample)
                for(i in 1:length(LevelsWant)){
                        comp$Sample <- relevel(comp$Sample, LevelsWant[i])
                }
                
                comp <- comp[c(1,4,5)]
                names(comp)[c(2,3)] <- c("SizeFactor_DESeq", "TotalAmplicons_relAb")
                comp <- tidyr::gather(comp, key = Corrector, value = NormedValue, -Sample)
                Tr <- ggplot(comp, aes(x = Sample, y = NormedValue, color = Corrector))
                Tr <- Tr + geom_point(size = 2) +
                        xlab("") +
                        ylab("correction value (normalized to median)") +
                        ggtitle(paste("Median SF: ", round(median(SFs),3), " Median TotAmp: ", round(median(sample_sums(physeq)),3), sep = "")) +
                        scale_color_manual(values = cbPalette[c(4,2)]) +
                        theme_bw() +
                        theme(panel.grid.minor = element_blank(),
                              panel.grid.major.y = element_blank(),
                              panel.grid.major.x = element_line(color = "#999999", size = .15),
                              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
                              legend.title = element_blank())
                
                # -- 3c: save Histograms of the SFs and sample_sums
                histo2 <- function(x, xtitle, gtitle) {
                        x <- data.frame(x = x)
                        Tr <- ggplot(x, aes(x = x))
                        Tr <- Tr + geom_histogram(binwidth = diff(range(x))/60, col = "black", fill = "#E69F00") +
                                geom_rug() +
                                geom_vline(xintercept = median(x$x), col = "#009E73", size = 1) +
                                ylab("No Samples") + 
                                xlab(xtitle) +
                                ggtitle(gtitle) +
                                theme_bw() + 
                                theme(panel.grid.minor = element_blank(),
                                      panel.grid.major.y = element_blank(),
                                      panel.grid.major.x = element_line(color = "#999999", size = .15))
                }
                
                Tr2 <- histo2(SFs, xtitle = "Size Factors", gtitle = "Size Factors a la DESeq")
                Tr3 <- histo2(sample_sums(physeq), xtitle = "Total Amplicons", gtitle = "Size Factors a la relative abundance")
                
                List <- list(Physeq = phynew, SFs = SFs, SizeFactorCompare = Tr, SFHistos = list(Tr2, Tr3), RefSample = GM, RatioMatrix = RatioMatrix, zeroPC = zeroPC, HistoList = HistoList)
                
                
        } else {
                
                phynew
                
        }
        
        
}