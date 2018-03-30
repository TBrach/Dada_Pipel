#######################################
### simulate_counts_withTPTaxa##
#################
## simulates a count table from a given template physeq where a defined number of true positive taxa will
# be more abundant either in grp1 or in grp2. The degree of abundance enrichment is defined by foldeffects. 
## requires: phyloseq, uses simulate_count_samples (NB for all physeqs generated taxa_are_rows = FALSE)
## Input: 
# - templatelist: a list with physeq template objects named with sampletypes (I usually use only one sampletype, 
# so only one template)
# - simparams: character vector with the simulation parameters to use named in simparamslabels, specifically
#     - n: the median number of liberary sizes (TotalAmplicons) used to sample from for the simulated samples
#     - min and maxEffectSize: the fold effects will be sampled from this range
#     - J1 and J2: the number of samples that will be generated for grp1 and grp2 respectively
# - NB: uses more inputs that should be defined before in the code
## Output: a list with the simulated phyloseq objects containing the TP taxa

simulate_counts_withTPTaxa <- function(templatelist, simparams, simparamslabels, nTP){
        simlist <- list()
        counter <- 0
        for (i in simparams) {
                counter <- counter + 1
                
                params = strsplit(i, comdelim)[[1]]
                names(params) <- simparamslabels
                
                # Initialize/reset to NULL
                n = sim = sim1 = sim2 = n1 = n2 = NULL
                
                n = as.integer(params["nreads"])
                sampletypei = params["SampleType"]
                
                # The number of samples for group 1
                Ji1 = as.integer(params["nsamplesgr1"])
                # The number of samples for group 2
                Ji2 = as.integer(params["nsamplesgr2"])
                # the actual physeq template
                templatei = templatelist[[sampletypei]]
                
                ##.## NB he has now a security loop because rarely a simulation fails
                ##.## see in his code the line starting with:
                # Rarely a simulation has a weird value and fails.
                ##.## I work without this security here
                
                # == sample LibrarySizes (TotalAmplicons) for the simulated samples ==
                scaledSums = round(n*(sample_sums(templatei)/median(sample_sums(templatei))))
                # in case of n = median(sample_sums(physeqi)) = default, scaledSums = sample_sums(templatei)
                n1   = sample(scaledSums, size = Ji1, replace = TRUE)
                n2   = sample(scaledSums, size = Ji2, replace = TRUE)
                # ====
                
                # ----  Set the TP taxa ----
                
                TPTaxa = sample(taxa_names(templatei), nTP, replace=FALSE)
                # randomly decide how many TruePositives will be more abundant in grp1:
                upGr1 <- sample(nTP, 1)
                # see which TP will be more abundant in Gr1 and which less abundant in Grp1 (more abundant in grp2)
                TPTaxa1 <- sample(TPTaxa, upGr1)
                TPTaxa2 <- setdiff(TPTaxa, TPTaxa1)
                
                wh.TP1 = taxa_names(templatei) %in% TPTaxa1
                wh.TP2 = taxa_names(templatei) %in% TPTaxa2
                
                # sample the foldeffects for the different TP so you can put the FE into the TP taxa names
                # NB: I need if statement to check whether minEffectSize<maxEffectSize because sample samples from 1:x if x 
                # is a single numeric, and since 3:3 = 3 would produce values between 1 and 3 instead of all 3!!
                if(as.numeric(params["minEffectSize"]) < as.numeric(params["maxEffectSize"])){
                        effectsizes1 <- sample(params["minEffectSize"]:params["maxEffectSize"],
                                               size = length(TPTaxa1), replace = TRUE)
                        effectsizes2 <- base::sample(params["minEffectSize"]:params["maxEffectSize"],
                                                     size = length(TPTaxa2), replace = TRUE)
                } else if(as.numeric(params["minEffectSize"]) == as.numeric(params["maxEffectSize"])) {
                        effectsizes1 <- rep(as.integer(params["minEffectSize"]), length(TPTaxa1))
                        effectsizes2 <- rep(as.integer(params["minEffectSize"]), length(TPTaxa2))
                } else {
                        stop("Something is wrong with your foldeffect choice!")
                }
                
                
                newname1 = paste0(taxa_names(templatei)[wh.TP1], "-TP-U_", effectsizes1)
                # U for upregulated (more abundant) in grp1
                newname2 = paste0(taxa_names(templatei)[wh.TP2], "-TP-D_", effectsizes2)
                # D for downregulated (less abundant) in grp1
                
                taxa_names(templatei)[wh.TP1] <- newname1
                taxa_names(templatei)[wh.TP2] <- newname2
                
                # --------
                
                # ----  Simulate the data separately for grp1 and grp2 using simulate_count_samples function ----
                
                sim1 = simulate_count_samples(postfix = paste0(sampletypei, ";grp1"), physeq = templatei, NoSamples = Ji1, SampleSizes = n1,
                                              TPNames = newname1, FEffects = effectsizes1)
                sim2 = simulate_count_samples(postfix = paste0(sampletypei, ";grp2"), physeq = templatei, NoSamples = Ji2, SampleSizes = n2,
                                              TPNames = newname2, FEffects = effectsizes2)
                
                # -------------
                
                if( is.null(sim1) || is.null(sim2) || 
                    is.null(n1) || is.null(n2)) {
                        stop("Error found during simulation. Need to investigate cause.")  
                }
                
                sim = merge_phyloseq(sim1, sim2) 
                
                # NB: the merge_phyloseq changes the typeof(otu_table(sim)) to double! Therefore
                # I change to integer here (because DESEQ uses integers), but note costs a second per 50 simulations with 20 samples
                Mat <- matrix(as.integer(round(as(otu_table(sim), "matrix"))), nrow = nsamples(sim))
                rownames(Mat) <- rownames(otu_table(sim))
                colnames(Mat) <- colnames(otu_table(sim))
                otu_table(sim) <- otu_table(Mat, taxa_are_rows = FALSE)
                
                simlist[[i]] <- sim
        }
        names(simlist) <- simparams
        simlist
}

#######################################
### simulate_count_samples
#################
## simulates count samples using sparsity_subsample given parameters defined within simulate_counts_withTPTaxa
## Input: 
# - postfix: just for labelling
# - physeq: a phyloseq object
# - NoSamples: how many simulated samples shall the output phyloseq object have
# - SampleSizes: the number of total reads of the simulated samples (n of length 1, so all
# NoSamples samples will have the same length, or n of length NoSamples so defining length of each
# sample) 
# TPNames: the names of the taxa that will be upregulated (over abundant)
# FEffects: the TPTaxa will be upregulated by the given FEffects
## Output: a phyloseq object with the simulated otu table, some very simple sample data
# to understand the simulation parameters, and the tax_table from physeq

simulate_count_samples = function(postfix="sim", physeq, NoSamples, SampleSizes, TPNames, FEffects){
        # ---- generate the urn matrix ----
        
        pi = taxa_sums(physeq)
        
        # add the effect at the TPTaxa
        pi[which(names(pi) %in% TPNames)] <- pi[which(names(pi) %in% TPNames)]*FEffects
        
        # new (22.04.2017): correct pi with NoSamples/nsamples(physeq) ---
        # idea behind: the urn should be the sum of as many real samples as samples will be simulated
        pi <- round(pi*(NoSamples/nsamples(physeq)))
        
        ##.## NB: smaller numbers in pi also SPEED UP THE SIMULATION
        ## BECAUSE THE RAREFACTION_SUBSAMPLE FUNCTION BELOW IS MUCH FASTER
        ## THE SMALLER THE NUMBERS IN PI WITHOUT CHANGING THE RESULT REALLY, SEE
        ## MAIN TEXT.
        
        if(length(NoSamples) != 1){ stop("Length of NoSamples should be 1.") }
        if(length(SampleSizes) != 1 && length(SampleSizes) != NoSamples){
                stop("SampleSizes should be length 1, or equal to NoSamples.")
        }
        
        # == change pi into a matrix to deal with sparsity problem ==
        
        # New 23.04.2017: Remember old version produced data with very low sparsity
        piMat <- matrix(rep(pi, NoSamples), ncol = NoSamples)
        # each column will be used as an urn
        # Determine the sparsity of each taxa 
        if(taxa_are_rows(physeq)) { physeq <- t(physeq) }
        TaxaSpars <- colSums(otu_table(physeq) == 0)/nrow(otu_table(physeq))
        NoZerosPerTaxa <- floor(TaxaSpars*NoSamples)
        # put zeros into the pimatrix
        ZeroIndexList <- lapply(NoZerosPerTaxa, FUN = function(x){sample(NoSamples,x)})
        piMat <- t(sapply(1:nrow(piMat), function(e){piMat[e,ZeroIndexList[[e]]] = 0;
        piMat[e,]}))
        # the next line is for the special case that NoSamples is = 1
        if(ncol(piMat) != NoSamples && nrow(piMat) == NoSamples){piMat <- t(piMat)}
        if(ncol(piMat) != NoSamples){ stop("NoSamples unequal to the number of columns in piMat") }
        # result: piMat with as many zeros (roughly) as the template physeq
        
        # ----    simulate the count table using sparsity_subsample    ----  
        
        simct = sapply(1:NoSamples, function(i){sparsity_subsample(piMat, SampleSizes, i)})
        # NB: only works when 
        
        # --------
        
        # ---- transform into physeq object ----
        
        # make samples rows and taxa columns (just preference)
        simct = t(simct)
        
        colnames(simct) <- names(pi)
        # Add new simulated sample_names to the row (sample) indices
        # old: rownames(simct) <- paste0(1:NoSamples, "::", 1:nrow(simct), postfix, sep="")
        rownames(simct) <- paste0("SimSample", "::", 1:nrow(simct), "_", postfix, sep="")
        
        # Put simulated abundances together with metadata as a phyloseq object
        OTU = otu_table(simct, taxa_are_rows=FALSE)
        
        # Define sample_data
        SDF = data.frame(sample = sample_names(OTU), TableNumber = 1:NoSamples, type = "simulated",
                         postfix = postfix)
        rownames(SDF) <- sample_names(OTU)
        SD  = sample_data(SDF)
        
        OutputPhyseq <- phyloseq(OTU, SD, tax_table(physeq))
        
}

#######################################
### sparsity_subsample##
#################
## Function very similar to phyloseq:::rarefaction_subsample, however with the difference that replace is
# always true in the sample command. 
# The new thing (23.04.2017) is that a Probability Matrix is given, each column is a probability urn,
# the index decides which column is chosen
## Input
# ProbMatrix: columns define the probabilities with which the taxa will be generated
# sample.sizes: vector with samples sizes to draw from
# index: decides which column in ProbMat and sample.size will be used

# NB: got this from see body(phyloseq:::rarefaction_subsample), the replace = TRUE version
# with replace = T you can also generate samples with higher count numbers than pi (sum(pi)).
# it is much faster also, see also replace in ?rarefy_even_depth


sparsity_subsample <- function (ProbMat, sample.sizes, index) {
        # attention, only works if nrow(ProbMat) is not NULL
        simsample <- integer(nrow(ProbMat)) 
        suppressWarnings(draws <- sample(1:nrow(ProbMat), sample.sizes[index], replace = TRUE, prob = ProbMat[,index]))
        drawtable <- table(draws)
        simsample[as(names(drawtable), "integer")] <- drawtable
        return(simsample)
}




#######################################
### overviewSimulation##
#################
## 
## Input: 
# - templatelist: a list with physeq template objects named with sampletypes (I usually use only one sampletype, 
# so only one template)
# - simparams: character vector with the simulation parameters to use named in simparamslabels, specifically
#     - n: the median number of liberary sizes (TotalAmplicons) used to sample from for the simulated samples
#     - min and maxEffectSize: the fold effects will be sampled from this range
#     - J1 and J2: the number of samples that will be generated for grp1 and grp2 respectively
# - NB: uses more inputs that should be defined before in the code
## Output: a list with the simulated phyloseq objects containing the TP taxa

## record interesting data about the different simulations, i.e. the number of TP, how many upregulated (grp1) and downregulated (upregulated in grp2), and the percentage of 0s in the table, compared to the input to estimate sparsity differences
# first keep some info of the template (NB: directly TEMPLATE here not templatelist)

overviewSimulation <- function(template, simlist, comdelim, plot = TRUE){
        
        # collect Overview of general data -----------------------
        
        if(taxa_are_rows(template)) { template <- t(template) } # from the tidy data point of view
        # I prefer taxa_are_rows = FALSE so rows (= observations = samples), and colums = variables = taxa
        # even though the tidy would be Sample column, Taxa column, Count column!
        # on top in the simulations taxa_are_rows = FALSE
        
        OTUTT <- as(otu_table(template), "matrix")
        OTUTTNA <- OTUTT
        OTUTTNA[OTUTTNA == 0] <- NA
        # also get a relative abundance Table
        OTUTT_RA <- OTUTT/rowSums(OTUTT)
        OTUTTNA_RA <- OTUTT_RA
        OTUTTNA_RA[OTUTTNA_RA == 0] <- NA
        
        
        Overview <- data.frame(NoSamples = nsamples(template),
                               NoTaxa = ntaxa(template),
                               MedianSampleSum = round(median(sample_sums(template))),
                               Sparsity = round(100*(sum(otu_table(template) == 0)/(ntaxa(template)*nsamples(template))), 3),
                               MaxCount = max(OTUTT),
                               nTP = 0, 
                               nTP_U = 0, 
                               nTP_D = 0,
                               MedianFE = 0,
                               # ---- test of taxa and sample variations between template and simulations ----
                               # == test of variation of the taxa over all samples ==
                               MedianTaxaSD = round(median(apply(OTUTT, 2, sd), na.rm = TRUE), 3),
                               MedTaxaSDNoZ = round(median(apply(OTUTTNA, 2, sd, na.rm = TRUE), na.rm = TRUE), 3),
                               # na.rm = TRUE in case in some Taxa was only present in only 1 sample
                               
                               # # NB: SD correlates with taxa_sums (= colSums(OTUTT))
                               # # plot(colSums(OTUTT), apply(OTUTT, 2, sd))
                               # # in principle you could correct by taxa_sums
                               # MedianTaxaSD_Cor = round(median(apply(OTUTT*(max(colSums(OTUTT))/colSums(OTUTT)), 2, sd), na.rm = TRUE), 3),
                               # MedianTaxaSDNoZ_Cor = round(median(apply(OTUTTNA*(max(colSums(OTUTTNA))/colSums(OTUTTNA)), 2, sd, na.rm = TRUE), na.rm = TRUE), 3),
        
                               # get same taxa variation for relative abundance table
                               MedianTaxaSD_RA = round(median(100*apply(OTUTT_RA, 2, sd), na.rm = TRUE), 3),
                               MedTaxaSDNoZ_RA = round(median(100*apply(OTUTTNA_RA, 2, sd, na.rm = TRUE), na.rm = TRUE), 3),
                               # # NB again: SD correlates with taxa_sums (= colSums(OTUTT))
                               # # plot(colSums(OTUTT_RA), apply(OTUTT_RA, 2, sd))
                               # # but I do not think correction over taxa_sums makes sense
                               # =====
                               
                               # == test of variation of the samples over all taxa ==
                               MedianSampleSD = round(median(apply(OTUTT, 1, sd), na.rm = TRUE), 3),
                               MedianSampleSDNoZ = round(median(apply(OTUTTNA, 1, sd, na.rm = TRUE), na.rm = TRUE), 3),
                               # makes in principle only sense on relative abundance
                               MedianSampleSD_RA = round(median(100*apply(OTUTT_RA, 1, sd), na.rm = TRUE), 3),
                               MedianSampleSDNoZ_RA = round(median(100*apply(OTUTTNA_RA, 1, sd, na.rm = TRUE), na.rm = TRUE), 3)
                               )
        
        
        Mat <- t(sapply(simlist, function(sim){
                NSamples <- nsamples(sim)
                NTaxa <- ntaxa(sim)
                nTPs <- sum(grepl("TP", taxa_names(sim)))
                nTPU <- sum(grepl("TP-U", taxa_names(sim)))
                nTPD <- sum(grepl("TP-D", taxa_names(sim)))
                OTUT <- as(otu_table(sim), 'matrix')
                OTUT_RA <- OTUT/rowSums(OTUT)
                Sparsity <- round(((sum(OTUT == 0)/(dim(OTUT)[1]*dim(OTUT)[2]))*100), 2)
                MedianSampleSum <- round(median(sample_sums(sim)))
                MaxCount <- max(OTUT)
                
                TPList <- strsplit(taxa_names(sim)[grepl("TP", taxa_names(sim))],"_")
                MedianFE <- median(as.numeric(sapply(TPList, function(TPName){TPName[[length(TPName)]]})))
                
                MedianTaxaSD <- round(median(apply(OTUT, 2, sd), na.rm = TRUE), 3)
                MedianTaxaSD_RA <- round(median(100*apply(OTUT_RA, 2, sd), na.rm = TRUE), 3)
                MedianSampleSD <- round(median(apply(OTUT, 1, sd), na.rm = TRUE), 3)
                MedianSampleSD_RA <- round(median(100*apply(OTUT_RA, 1, sd), na.rm = TRUE), 3)
                OTUT[OTUT == 0] <- NA
                OTUT_RA[OTUT_RA == 0] <- NA
                MedTaxaSDNoZ <- round(median(apply(OTUT, 2, sd, na.rm = TRUE), na.rm = TRUE), 3)
                MedTaxaSDNoZ_RA <- round(median(100*apply(OTUT_RA, 2, sd, na.rm = TRUE), na.rm = TRUE), 3)
                MedianSampleSDNoZ <- round(median(apply(OTUT, 1, sd, na.rm = TRUE), na.rm = TRUE), 3)
                MedianSampleSDNoZ_RA = round(median(100*apply(OTUT_RA, 1, sd, na.rm = TRUE), na.rm = TRUE), 3)
                
                c(NSamples, NTaxa, MedianSampleSum, Sparsity, MaxCount, nTPs, nTPU, nTPD, MedianFE,
                  MedianTaxaSD, MedTaxaSDNoZ, MedianTaxaSD_RA, MedTaxaSDNoZ_RA, MedianSampleSD, MedianSampleSDNoZ, 
                  MedianSampleSD_RA, MedianSampleSDNoZ_RA)
        }))
        
        colnames(Mat) <- colnames(Overview)
        Overview <- rbind(Overview, Mat)
        
        if (plot) {
                
                # ---- plots ----
                # ==  Tr: illustration that Sample SD correlates with Sample sum  ==
                # so in consequence you should use relative abundances to compare Sample SDs
                
                cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
                               "#0072B2", "#D55E00", "#CC79A7")
                # SampleSD <- 100*apply(OTUTT, 1, sd)
                # Temp <- data.frame(SSums = sample_sums(template), SSDs = SampleSD)
                # Tr <- ggplot(Temp, aes(x = SSums, y = SSDs))
                # Tr <- Tr + geom_point(size = 3, color = cbPalette[2]) +
                #         geom_smooth(method = 'lm') +
                #         theme_bw() +
                #         xlab("Sample Sums Template") +
                #         ylab("SD Samples")
        
                
                # == Tr2: Histogram of template sample SD on relative abundances == 
                Temp <- data.frame(Sample = sample_names(template), Sample_SD_RA = 100*apply(OTUTT_RA, 1, sd))
                Tr2 <- ggplot(Temp, aes(x = Sample_SD_RA))
                Tr2 <- Tr2 + geom_histogram(binwidth = diff(range(Temp$Sample_SD_RA))/30, fill = cbPalette[2]) +
                        geom_rug() +
                        geom_vline(xintercept = median(Temp$Sample_SD_RA), col = "#009E73", lty = "dashed") +
                        ggtitle(paste("Sample SD (RA), distribution of ", nsamples(template), " template samples (median green line)", sep = "")) +
                        theme_bw() +
                        xlab("Sample_SD_RA")
                
                # # == Tr3: Histogram of template sample SD on relative abundances excluding zeros == 
                # Temp <- data.frame(Sample = sample_names(template), Sample_SD_RA_NZ = 100*apply(OTUTTNA_RA, 1, sd, na.rm = TRUE))
                # Tr3 <- ggplot(Temp, aes(x = Sample_SD_RA_NZ))
                # Tr3 <- Tr3 + geom_histogram(binwidth = diff(range(Temp$Sample_SD_RA_NZ))/30, fill = cbPalette[2]) +
                #         geom_rug() +
                #         geom_vline(xintercept = median(Temp$Sample_SD_RA_NZ), col = "#009E73", lty = "dashed") +
                #         ggtitle(paste("Sample SD (RA, NZ), distribution of ", nsamples(template), " template samples (median green line)", sep = "")) +
                #         theme_bw() +
                #         xlab("Sample_SD_RA_NZ")
                
                # == get a similar histogram of the Sample SD (RA) for each similation ==
                
                TrList <- lapply(simlist, function(sim){
                        OTUT_RA <- as(otu_table(sim), 'matrix')/sample_sums(sim)
                        Simu <- data.frame(Sample = sample_names(sim), Sample_SD_RA = 100*apply(OTUT_RA, 1, sd))
                        Tr <- ggplot(Simu, aes(x = Sample_SD_RA))
                        Tr <- Tr + geom_histogram(binwidth = diff(range(Simu$Sample_SD_RA))/30, fill = cbPalette[2]) +
                                geom_rug() +
                                geom_vline(xintercept = median(Simu$Sample_SD_RA), col = "#009E73", lty = "dashed") +
                                geom_vline(xintercept = median(Temp$Sample_SD_RA), col = cbPalette[7], lty = "dashed") +
                                ggtitle(paste("Sample SD (RA), distribution of ", nsamples(sim), " template samples (median green line)", sep = "")) +
                                theme_bw() +
                                xlab("Sample_SD_RA")
                        Tr  
                })
                
                
                # == Max Values in the template Samples ==
                
                OutTemp <- data.frame(x = 1:nsamples(template), Max = sort(apply(OTUTT, 1, max)))
                ot <- ggplot(OutTemp, aes(x = x, y = Max))
                ot <- ot + 
                        geom_hline(yintercept = max(OutTemp$Max), lty = "dashed", col = cbPalette[4]) + 
                        geom_hline(yintercept = median(OutTemp$Max), lty = "dashed", col = cbPalette[6]) +
                        geom_point(size = 3, color = cbPalette[2]) +
                        xlab("Sample") +
                        ylab("Max Count") +
                        ggtitle("Max counts in template samples") +
                        theme_bw()
                
                
                TrListMax <- lapply(simlist, function(sim){
                        OTUT <- as(otu_table(sim), 'matrix')
                        OutSim <- data.frame(x = 1:nsamples(sim), Max = sort(apply(OTUT, 1, max)))
                        Tr <- ggplot(OutSim, aes(x = x, y = Max))
                        Tr <- Tr + 
                                geom_hline(yintercept = max(OutSim$Max), lty = "dashed", col = cbPalette[4]) + 
                                geom_hline(yintercept = median(OutTemp$Max), lty = "dashed", col = cbPalette[7]) +
                                geom_hline(yintercept = median(OutSim$Max), lty = "dashed", col = cbPalette[6]) +
                                geom_point(size = 3, color = cbPalette[2]) +
                                xlab("Sample") +
                                ylab("Max Count") +
                                ggtitle("Max counts in simulation") +
                                theme_bw()
                        Tr  
                })
                
                
                ## ==== Heatmap plots template ====
                
                # NB: the goal was here to have 0s as a special color (red, or white)
                # then to have the option to give Count limits, so that every count above that value gets the max colour
                # The latter point can be done with the scale package and the option: oob = swish
                # please read here: <https://github.com/tidyverse/ggplot2/issues/866>
                # see also the na.value option not used

                DF_OTUTT <- as.data.frame(t(OTUTT))
                colnames(DF_OTUTT) <- paste0("S", 1:nsamples(template))
                DF_OTUTT$Taxa <- rownames(DF_OTUTT)
                DF_OTUTT <- tidyr::gather(DF_OTUTT, key = Sample , value = Count, - Taxa)
                DF_OTUTT$Taxa <- as.factor(DF_OTUTT$Taxa)
                LW <- taxa_names(template)
                for(z in 1:ntaxa(template)){
                        DF_OTUTT$Taxa <- relevel(DF_OTUTT$Taxa, ref = LW[z])
                }
                DF_OTUTT$Sample <- as.factor(DF_OTUTT$Sample)
                LW <- rev(paste0("S", 1:nsamples(template)))
                for(z in 1:nsamples(template)){
                        DF_OTUTT$Sample <- relevel(DF_OTUTT$Sample, ref = LW[z])
                }
                
                
                color.palette=viridis
                hmt <- ggplot(DF_OTUTT, aes(x=Sample, y=Taxa, fill=Count))
                hmt <- hmt + 
                        geom_raster() + # supposedly works faster than geom_tile
                        #scale_fill_viridis(limits = c(0,300), oob=squish)+
                        scale_fill_gradientn(limits = c(0,300), colors = c("red",color.palette(5)),values = c(0, 1e-12, 0.25, 0.5, 0.75, 1), oob = squish) +
                        scale_x_discrete(position = "top") +
                        #coord_equal() +
                        labs(x=NULL, y=NULL, title="Template") +
                        theme_tufte(base_family="Helvetica") +
                        theme(plot.title=element_text(hjust=0)) +
                        theme(axis.ticks=element_blank()) +
                        theme(axis.text=element_text(size=7)) +
                        theme(legend.title=element_text(size=8)) +
                        theme(legend.text=element_text(size=6))
                
                DF_OTUTT <- dplyr::group_by(DF_OTUTT, Sample)
                DF_OTUTT <- dplyr::mutate(DF_OTUTT, Count = Count/sum(Count))
                
                hmtRA <- ggplot(DF_OTUTT, aes(x=Sample, y=Taxa, fill=Count))
                hmtRA <- hmtRA + 
                        geom_raster() + # supposedly works faster than geom_tile
                        #scale_fill_viridis(limits = c(0,300), oob=squish)+
                        scale_fill_gradientn(limits = c(0,.2), colors = c("red",color.palette(5)),values = c(0, 1e-12, 0.15, 0.3, 0.45, 1), oob = squish) +
                        scale_x_discrete(position = "top") +
                        #coord_equal() +
                        labs(x=NULL, y=NULL, title="Template") +
                        theme_tufte(base_family="Helvetica") +
                        theme(plot.title=element_text(hjust=0)) +
                        theme(axis.ticks=element_blank()) +
                        theme(axis.text=element_text(size=7)) +
                        theme(legend.title=element_text(size=8)) +
                        theme(legend.text=element_text(size=6))
                
                
                ## ==== Heatmap plots simulation ====
                
                hmList <- vector(mode = 'list', length = length(simlist))
                hmRAList <- vector(mode = 'list', length = length(simlist))
                hmListTP <- vector(mode = 'list', length = length(simlist))
                hmRAListTP <- vector(mode = 'list', length = length(simlist))
                
                
                for (i in 1:length(simlist)){
                        ct <- t(as(otu_table(simlist[[i]]), "matrix"))
                        ct <- as.data.frame(ct)
                        #rownames(ct) <- paste0("T", 1:ntaxa(simlist[[i]]))
                        colnames(ct) <- paste0("S", 1:nsamples(simlist[[i]]))
                        ct$Taxa <- rownames(ct)
                        ct <- tidyr::gather(ct, key = Sample , value = Count, - Taxa)
                        ct$Taxa <- as.factor(ct$Taxa)
                        LW <- taxa_names(simlist[[i]])
                        for(z in 1:ntaxa(simlist[[i]])){
                                ct$Taxa <- relevel(ct$Taxa, ref = LW[z])
                        }
                        ct$Sample <- as.factor(ct$Sample)
                        LW <- rev(paste0("S", 1:nsamples(simlist[[i]])))
                        for(z in 1:nsamples(simlist[[i]])){
                                ct$Sample <- relevel(ct$Sample, ref = LW[z])
                        }
                        # add color to the taxa so you see up and down
                        colyaxis <- vector(mode = "character", length = ntaxa(simlist[[i]]))
                        colyaxis[] <- "black"
                        colyaxis[grepl("TP-U", levels(ct$Taxa))] <- "#E69F00"
                        colyaxis[grepl("TP-D", levels(ct$Taxa))] <- "#009E73"
                        
                        ngrp1 <- unlist(strsplit(names(simlist)[i], comdelim))[6]
                        colxaxis <- vector(mode = "character", length = nsamples(simlist[[i]]))
                        colxaxis[] <- "#009E73"
                        colxaxis[1:ngrp1] <- "#E69F00"
                        
                        hms <- ggplot(ct, aes(x=Sample, y=Taxa, fill=Count))
                        hms <- hms + 
                                geom_raster() + # supposedly works faster than geom_tile
                                #scale_fill_viridis(limits = c(0,300), oob=squish)+
                                scale_fill_gradientn(limits = c(0,300), colors = c("red",color.palette(5)),values = c(0, 1e-12, 0.25, 0.5, 0.75, 1), oob = squish) +
                                #coord_equal() +
                                scale_x_discrete(position = "top") +
                                labs(x=NULL, y=NULL, title= paste("Simulation ", i, sep = "")) +
                                theme_tufte(base_family="Helvetica") +
                                theme(plot.title=element_text(hjust=0)) +
                                theme(axis.ticks=element_blank()) +
                                theme(axis.text=element_text(size=7)) +
                                theme(legend.title=element_text(size=8)) +
                                theme(legend.text=element_text(size=6)) +
                                theme(axis.text.y = element_text(colour = colyaxis),
                                      axis.text.x = element_text(colour = colxaxis))
                        
                        hmList[[i]] <- hms
                        
                        ctRA <- dplyr::group_by(ct, Sample)
                        ctRA <- dplyr::mutate(ctRA, Count = Count/sum(Count))
                        
                        hmsRA <- ggplot(ctRA, aes(x=Sample, y=Taxa, fill=Count))
                        hmsRA <- hmsRA + 
                                geom_raster() + # supposedly works faster than geom_tile
                                #scale_fill_viridis(limits = c(0,300), oob=squish)+
                                scale_fill_gradientn(limits = c(0,.2), colors = c("red",color.palette(5)),values = c(0, 1e-12, 0.15, 0.3, 0.45, 1), oob = squish) +
                                #coord_equal() +
                                scale_x_discrete(position = "top") +
                                labs(x=NULL, y=NULL, title= paste("Simulation ", i, sep = "")) +
                                theme_tufte(base_family="Helvetica") +
                                theme(plot.title=element_text(hjust=0)) +
                                theme(axis.ticks=element_blank()) +
                                theme(axis.text=element_text(size=7)) +
                                theme(legend.title=element_text(size=8)) +
                                theme(legend.text=element_text(size=6)) +
                                theme(axis.text.y = element_text(colour = colyaxis),
                                      axis.text.x = element_text(colour = colxaxis))
                        
                        hmRAList[[i]] <- hmsRA
                        
                        # now only plot the TP --
                        ctTP <- dplyr::filter(ct, grepl("TP", Taxa))
                        colyaxis2 <- colyaxis[colyaxis != "black"]
                        hmsTP <- ggplot(ctTP, aes(x=Sample, y=Taxa, fill=Count))
                        hmsTP <- hmsTP + 
                                geom_raster() + # supposedly works faster than geom_tile
                                #scale_fill_viridis(limits = c(0,300), oob=squish)+
                                scale_fill_gradientn(limits = c(0,300), colors = c("red",color.palette(5)),values = c(0, 1e-12, 0.25, 0.5, 0.75, 1), oob = squish) +
                                #coord_equal() +
                                scale_x_discrete(position = "top") +
                                labs(x=NULL, y=NULL, title= paste("Simulation ", i, sep = "")) +
                                theme_tufte(base_family="Helvetica") +
                                theme(plot.title=element_text(hjust=0)) +
                                theme(axis.ticks=element_blank()) +
                                theme(axis.text=element_text(size=7)) +
                                theme(legend.title=element_text(size=8)) +
                                theme(legend.text=element_text(size=6)) +
                                theme(axis.text.y = element_text(colour = colyaxis2),
                                      axis.text.x = element_text(colour = colxaxis))
                        
                        hmListTP[[i]] <- hmsTP
                        
                        ctTPRA <- dplyr::filter(ctRA, grepl("TP", Taxa))
                        hmsTPRA <- ggplot(ctTPRA, aes(x=Sample, y=Taxa, fill=Count))
                        hmsTPRA <- hmsTPRA + 
                                geom_raster() + # supposedly works faster than geom_tile
                                #scale_fill_viridis(limits = c(0,300), oob=squish)+
                                scale_fill_gradientn(limits = c(0,.2), colors = c("red",color.palette(5)),values = c(0, 1e-12, 0.15, 0.3, 0.45, 1), oob = squish) +
                                #coord_equal() +
                                scale_x_discrete(position = "top") +
                                labs(x=NULL, y=NULL, title= paste("Simulation ", i, sep = "")) +
                                theme_tufte(base_family="Helvetica") +
                                theme(plot.title=element_text(hjust=0)) +
                                theme(axis.ticks=element_blank()) +
                                theme(axis.text=element_text(size=7)) +
                                theme(legend.title=element_text(size=8)) +
                                theme(legend.text=element_text(size=6)) +
                                theme(axis.text.y = element_text(colour = colyaxis2),
                                      axis.text.x = element_text(colour = colxaxis))
                        
                        hmRAListTP[[i]] <- hmsTPRA
                        
                } 
                
                names(hmList) <- names(simlist)
                names(hmListTP) <- names(simlist)
                names(hmRAList) <- names(simlist)
                names(hmRAListTP) <- names(simlist)
                
                Overvlist <- list(Overview = Overview, Hist_Template_SampleSD_RA = Tr2, Hist_Sims_SampleSD_RA = TrList,
                                  Template_MaxCounts = ot, Sims_MaxCounts = TrListMax, HM_Template = hmt, HM_RA_Template = hmtRA, HM_Sims = hmList, HM_RA_Sims = hmRAList,
                                  HM_Sims_TP = hmListTP, HM_RA_Sims_TP = hmRAListTP)
                
                
                
        } else {
                
                Overvlist <- list(Overview = Overview)
                
        }
        
        
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



#######################################
### adj_LS
#######################################
# implementation of DESeq2 library size similar to estimateSizeFactorsForMatrix
# in difference to adjust_LS (obsolete) ignore.zero.ratios is always TRUE (so 0 ratios are always ignored)

## Input:
# physeq
# zeros.count: zeros.count of gm_own, if TRUE 0 will be considered when calculating the geometric mean
# if FALSE not and thus the geometric means will be bigger (see gm_own)
# percentile: the percentile to select the size factor SF of a sample based on its count ratios to the reference sample. In
# DESeq percentile = 50, i.e. stats::median is used. 
# plots: if TRUE SFs and plots will be given in an output list, otherwise just the adjusted physeq is returned


adj_LS <- function(physeq, zeros.count = FALSE, percentile = 50, plots = FALSE)  {
        
        # ---- Step 1: Calculate Geometric mean for each taxa over all samples ------
        # NB: these GM is basically the reference sample
        if(taxa_are_rows(physeq)){
                GM <- apply(otu_table(physeq), 1, gm_own, zeros.count = zeros.count)   
        } else {
                GM <- apply(otu_table(physeq), 2, gm_own, zeros.count = zeros.count) 
        }
        
        # ---- Step 2: Calculate Size factors --------
        
        # NB: x/y = exp(log(x) - log(y))
        
        if (taxa_are_rows(physeq)) {
                SFs <- apply(as(otu_table(physeq), "matrix"), 2, function(sample_cnts){exp(quantile((log(sample_cnts)-log(GM))[sample_cnts > 0], probs = percentile/100, na.rm = T))})
                SFs <- SFs/exp(mean(log(SFs)))
        } else {
                SFs <- apply(as(otu_table(physeq), "matrix"), 1, function(sample_cnts){exp(quantile((log(sample_cnts)-log(GM))[sample_cnts > 0], probs = percentile/100, na.rm = T))})
                SFs <- SFs/exp(mean(log(SFs)))
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
                
                cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
                
                
                # compare calculated SFs to library sizes of the samples
                if(!identical(names(SFs), names(sample_sums(physeq)))){warning("Probably some Mix Up CHECK!")}
                comp <- data.frame(Sample = names(SFs), TotAmps = sample_sums(physeq), SFs = SFs)
                comp$SFsNormed <- comp$SFs/median(comp$SFs)
                # comp$TotAmpsNormed <- comp$TotAmps/mean(comp$TotAmps)
                comp$TotAmpsNormed <- comp$TotAmps/median(comp$TotAmps)
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
                
                List <- list(Physeq = phynew, SFs = SFs, SizeFactorCompare = Tr, SFHistos = list(Tr2, Tr3), RefSample = GM)
                
                
        } else {
                
                phynew
                
        }
}



#######################################
### mtApply
#######################################
# a simple lapply on a list of physeq objects to do phyloseq::mt and add further info
# phyloseq::mt is a phyloseq wrapper for mt.minP from "multtest" package
# I added a fisher exact test on the number of present and absent samples for each taxon
## Input
# simlist (list of phyloseqs)
# classlabel (as in mt.minP): refers to a factor in sample_data(phyloseq) that defines the two groups
# method and test of mt.minP or mt
## Output
# list of dataframes with the results

mtApply <- function(simlist, classlabel = "postfix",  test = "wilcoxon", method = "fdr"){
        lapply(simlist, FUN = function(sim){
                mat <- mt(sim, classlabel = classlabel, method = method, test = test) # a phyloseq wrapper for mt.minP from "multtest" package
                mat <- mat[,c("teststat", "rawp", "adjp", "fdr")] # see in explanation in next section: rawp and adjp come from a permutation method (https://github.com/joey711/phyloseq/issues/439), fdr is simpy p.adjust used on rawp
                # - add median counts of the groups -
                CT <- as(otu_table(sim), "matrix")
                groupFactor <- sample_data(sim)[[classlabel]]
                Median_grp1 <- apply(CT[grepl(levels(groupFactor)[1], rownames(CT)),], 2, median)
                Median_grp2 <- apply(CT[grepl(levels(groupFactor)[2], rownames(CT)),], 2, median)
                n1 <- sum(grepl(levels(groupFactor)[1], rownames(CT)))
                n2 <- sum(grepl(levels(groupFactor)[2], rownames(CT)))
                Mean_grp1 <- apply(CT[grepl(levels(groupFactor)[1], rownames(CT)),], 2, mean)
                Mean_grp2 <- apply(CT[grepl(levels(groupFactor)[2], rownames(CT)),], 2, mean)
                Zeros_grp1 <- apply(CT[grepl(levels(groupFactor)[1], rownames(CT)),], 2, function(cnts){sum(cnts == 0)})
                Zeros_grp2 <- apply(CT[grepl(levels(groupFactor)[2], rownames(CT)),], 2, function(cnts){sum(cnts == 0)}) 
                Present_grp1 <- n1 - Zeros_grp1
                Present_grp2 <- n2 - Zeros_grp2
                Sparsity_grp1 <- 100*(Zeros_grp1/n1)
                Sparsity_grp2 <- 100*(Zeros_grp2/n2)
                # -- add fisher exact test of presence differences (should be none in simulation) --
                Fisher <- t(sapply(1:ntaxa(sim), FUN = function(i){
                        fisherMat <- matrix(c(Present_grp1[i], Zeros_grp1[i], Present_grp2[i],
                                              Zeros_grp2[i]), ncol = 2, dimnames = list(c("Present", "Absent"), c("grp1", "grp2")))
                        Test <- fisher.test(fisherMat)
                        cbind(Test$p.value, Test$estimate)
                }))
                rownames(Fisher) <- taxa_names(sim)
                fdrFisher <- p.adjust(Fisher[,1], method = "fdr")
                cbind(id = rownames(mat), mat, Median_grp1 = Median_grp1[rownames(mat)], Median_grp2 = Median_grp2[rownames(mat)],
                      Mean_grp1 = Mean_grp1[rownames(mat)], Mean_grp2 = Mean_grp2[rownames(mat)], Present_grp1 = Present_grp1[rownames(mat)],
                      Present_grp2 = Present_grp2[rownames(mat)], Zeros_grp1 = Zeros_grp1[rownames(mat)], Zeros_grp2 = Zeros_grp2[rownames(mat)], 
                      Sparsity_grp1 = Sparsity_grp1[rownames(mat)], Sparsity_grp2 = Sparsity_grp2[rownames(mat)], pFisher = Fisher[rownames(mat),1],
                      fdrFisher = fdrFisher[rownames(mat)], oddsRatioFisher = Fisher[rownames(mat),2])
        })
}




#######################################
### wilcoxTestApply
#######################################
# a simple lapply on a list of physeq objects to do wilcoxon test and add further info on the physeq objects
# in addition it performs a fisher exact test on the sparsity proportions
# with excludeZeros you can decide on whether 0 counts should be excluded for the wilcox.test, the fisher sparsity test is
# of course not affected. 
# NB: in case in one of the two groups all counts are 0 and excludeZeros = T, one value in that group is set to the mean of
# all counts over all samples plus 1. So in case you have an all 0 taxon, you get p-value = NA, wilcox.test(x=1, y = 1)
# The provided Median and Mean values are affected by Exclude 0
# The teststatistic is based on the standardized teststatistic, equation provided by multtest::mt.minP (compare with mtApply)
# (see equation for standStat in the code)
## Input
# simlist (list of phyloseqs)
# classlabel (as in mt.minP): refers to a factor in sample_data(phyloseq) that defines the two groups
# excludeZeros: decides on whether 0s should be considered when comparing the groups in a wilcox.test
## Output
# list of dataframes with the results

wilcoxTestApply <- function(simlist, classlabel = "postfix", excludeZeros = FALSE) {
        lapply(simlist, FUN = function(sim){
                if(taxa_are_rows(sim)){sim <- t(sim)}
                mat <- as(otu_table(sim), "matrix")
                groupFactor <- sample_data(sim)[[classlabel]]
                
                res_mat <- apply(mat, 2, function(taxon_counts){
                        x <- taxon_counts[grepl(levels(groupFactor)[1], names(taxon_counts))]
                        Zeros_grp1 <- sum(x == 0)
                        Sparsity_grp1 <- 100*(Zeros_grp1/length(x))
                        Present_grp1 <- length(x)-Zeros_grp1
                        if(excludeZeros){
                                # in case there are only 0s in x
                                if(all(x == 0)){x[1] <- ceiling(mean(taxon_counts))+1} 
                                x <- x[x != 0]
                        }
                        Median_grp1 <- median(x, na.rm = T)
                        Mean_grp1 <- mean(x, na.rm = T)
                        y <- taxon_counts[grepl(levels(groupFactor)[2], names(taxon_counts))]
                        Zeros_grp2 <- sum(y == 0)
                        Present_grp2 <- length(y)-Zeros_grp2
                        Sparsity_grp2 <- 100*(Zeros_grp2/length(y))
                        if(excludeZeros){
                                if(all(y == 0)){y[1] <- ceiling(mean(taxon_counts))+1}
                                y <- y[y != 0]
                        }
                        Median_grp2 <- median(y, na.rm = T)
                        Mean_grp2 <- mean(y, na.rm = T)
                        wilcTest <- wilcox.test(x = x, y = y, alternative = "two", paired = F, exact = F)
                        pValue <- wilcTest$p.value
                        W <- wilcTest$statistic
                        # calculate standardized rank sum Wilcoxon statistics as in multtest::mt.minP
                        Ranks <- rank(c(x, y))
                        n1 <- length(x)
                        n2 <- length(y)
                        # Wx <- sum(Ranks[1:n1])-(n1*(n1+1)/2) # would be the same as W
                        standStat <- -1*((sum(Ranks[1:n1]) - n1*(n1+n2+1)/2)/sqrt(n1*n2*(n1+n2+1)/12))
                        # check that multtest::mt.minP would give the same statistic
                        # mati <- matrix(c(x,y), nrow = 1)
                        # grFac <- as.factor(c(rep("a",8), rep("b",8)))
                        # testStat <- multtest::mt.minP(mati, grFac, test = "wilcoxon")$teststat
                        # identical(standStat, testStat) # TRUE
                        
                        # -- add fisher exact test of presence differences (should be none in simulation) --
                        fisherMat <- matrix(c(Present_grp1, Zeros_grp1, Present_grp2, Zeros_grp2), 
                                            ncol = 2, dimnames = list(c("Present", "Absent"), c("grp1", "grp2")) )
                        Test <- fisher.test(fisherMat)
                        
                        c(teststat = standStat, rawp = pValue, Median_grp1 = Median_grp1, Median_grp2 = Median_grp2, 
                          Mean_grp1 = Mean_grp1, Mean_grp2 = Mean_grp2, n1 = n1, n2 = n2, Present_grp1 = Present_grp1, 
                          Present_grp2 = Present_grp2, Zeros_grp1 = Zeros_grp1, Zeros_grp2 = Zeros_grp2, 
                          Sparsity_grp1 = Sparsity_grp1, Sparsity_grp2 = Sparsity_grp2, W, 
                          pFisher = Test$p.value, Test$estimate)
                })
                res_mat <- t(res_mat)
                fdr <- p.adjust(res_mat[,"rawp"], method = "fdr")
                fdrFisher <- p.adjust(res_mat[,"pFisher"], method = "fdr")
                DF <- data.frame(id = rownames(res_mat), res_mat, adjp = fdr, fdr = fdr, fdrFisher = fdrFisher)
                DF <- dplyr::select(DF, 1:3, 19:20, 4:7, 10:15, 8:9, 16, 17, 21, 18)
                colnames(DF)[21] <- "oddsRatioFisher"
                DF <- dplyr::arrange(DF, desc(abs(teststat)))
        })
}




#######################################
### DESeq2Apply
#######################################
# a simple lapply on a list of physeq objects to do DESeq2 pipeline and add further info
## Input
# simlist (list of phyloseqs)
# classlabel (as in mt.minP): refers to a factor in sample_data(phyloseq) that defines the two groups
# method and test of mt.minP or mt
## Output
# list of dataframes with the results

DESeq2Apply <- function(simlist, classlabel = "postfix", type = "ratio"){
        lapply(simlist, FUN = function(physeq){
                if(taxa_are_rows(physeq)){stop("Taxa should not be rows in the simulations")}
                # change classlabel to factor so phyloseq_to_deseq2 does not need to do it
                sample_data(physeq)[[classlabel]] = as.factor(sample_data(physeq)[[classlabel]])
                DES = phyloseq_to_deseq2(physeq, formula(paste("~", classlabel)))
                if(type == "ratio"){
                        GM <- apply(otu_table(physeq), 2, gm_own, zeros.count = FALSE)
                }
                dds <- estimateSizeFactors(DES, type = type, geoMeans = GM)
                # NB: geoMeans is ignored when type = "iterate"
                # NB2: "iterate" takes much longer than "ratio", and when using GM via gm_own with "ratio" the size factors 
                # correlate with above 99% with size factors from "iterate"
                dds <- estimateDispersions(dds) # NOTE you could add quiet = TRUE
                dds <- nbinomWaldTest(dds)
                res <- as.data.frame(results(dds, contrast=c(classlabel,levels(sample_data(physeq)[[classlabel]])[1],levels(sample_data(physeq)[[classlabel]])[2])))
                # contrast had to be used to have grp1/grp2 instead of the other way around
                res$fdr <- p.adjust(res$pvalue, method = "fdr")
                
                groupFactor <- sample_data(physeq)[[classlabel]]
                CT <- t(counts(dds, normalized = TRUE))
                n1 <- sum(grepl(levels(groupFactor)[1], rownames(CT)))
                n2 <- sum(grepl(levels(groupFactor)[2], rownames(CT)))
                res$Median_grp1 <- apply(CT[grepl(levels(groupFactor)[1], rownames(CT)),], 2, median)
                res$Median_grp2 <- apply(CT[grepl(levels(groupFactor)[2], rownames(CT)),], 2, median)
                res$Mean_grp1 <- apply(CT[grepl(levels(groupFactor)[1], rownames(CT)),], 2, mean)
                res$Mean_grp2 <- apply(CT[grepl(levels(groupFactor)[2], rownames(CT)),], 2, mean)
                # res$baseMeanSelf <- apply(CT, 2, mean) # exactly the same as baseMean!
                res$Zeros_grp1 <- apply(CT[grepl(levels(groupFactor)[1], rownames(CT)),], 2, function(cnts){sum(cnts == 0)})
                res$Zeros_grp2 <- apply(CT[grepl(levels(groupFactor)[2], rownames(CT)),], 2, function(cnts){sum(cnts == 0)})
                res$Present_grp1 <- n1-res$Zeros_grp1
                res$Present_grp2 <- n2-res$Zeros_grp2
                res$Sparsity_grp1 <- 100*(res$Zeros_grp1/n1)
                res$Sparsity_grp2 <- 100*(res$Zeros_grp2/n2)

                # -- add fisher exact test of presence differences (should be none in simulation) --
                Fisher <- t(sapply(1:nrow(res), FUN = function(i){
                        fisherMat <- matrix(c(res$Present_grp1[i], res$Zeros_grp1[i], res$Present_grp2[i],
                                              res$Zeros_grp2[i]), ncol = 2, dimnames = list(c("Present", "Absent"), c("grp1", "grp2")))
                        Test <- fisher.test(fisherMat)
                        cbind(Test$p.value, Test$estimate)
                }))
                
                fdrFisher <- p.adjust(Fisher[,1], method = "fdr")
                
                res$pFisher <- Fisher[,1]
                res$fdrFisher <- fdrFisher
                res$oddsRatioFisher <- Fisher[,2]
                
                res <- dplyr::select(res, teststat = stat, rawp = pvalue, adjp = padj, fdr, Median_grp1, Median_grp2, Mean_grp1,
                                     Mean_grp2, Present_grp1, Present_grp2, Zeros_grp1, Zeros_grp2, Sparsity_grp1, Sparsity_grp2, baseMean, 
                                     log2FoldChange, lfcSE, pFisher, fdrFisher, oddsRatioFisher)
                res <- res[order(res$rawp),]
                
                cbind(id = rownames(res), res)
        })
}


####################################
## CalcTbTMatrixesNEW:
###################################
# Reasoning: it is all about the treatment of zeros. 
# The method is based on direct comparisons between Taxa, thus avoiding compositionality
# I think that Values where one of the two taxa that are compared is missing should not be used to determine
# differential abundance (instead they should be separately used to check for differences in sparsity of a taxon)
# so if the host taxon is 0 in a sample, then the entire sample will basically be ignored
# Values where the host taxon or the other taxon is 0 should be set to NA, so that all 0 values come
# from cases log(x/y) where x = y but both x and y != 0
# NB: read also explanations within the code

## Input: 
# - simlist = a list with physeq objects (usually the simulations) with count data
## Output: 
# - list (length = number of simulations in simlist) in which each entry 
#   is again a list (length = number of taxa in the simulation) with the 
#   normalised and log transformed TbTMatrixes for each host taxon 

CalcTbTMatrixesNEW = function(simlist){
        
        TbTMatrixes_List <- lapply(simlist, FUN <- function(physeq){
                if(taxa_are_rows(physeq)){stop("Taxa should not be rows in the simulations")}
                CM <- t(as(otu_table(physeq), 'matrix')) # well now taxa are rows and samples are columns
                
                # == Take in-sample taxa by taxa ratios and take the log ==
                # NB: log(x/y) = log(x) - log(y)
                TbTMatrixes <- lapply(1:nrow(CM), function(i){apply(CM, 2, function(samp_cnts){log(samp_cnts[i]) - log(samp_cnts)})})
                # produces for each taxon (= host taxon) a TbTMatrix
                # NB: there are -Inf, Inf, and NaN values in the matrixes, specifically
                # 0/x = log(0) - log(x) = -Inf, x/0 = log(x) - log(0) = Inf; 0/0 = log(0) - log(0) = NaN!
                # Inf, -Inf, and NaN will be ignored for calculating the rowMeans (geometric means) in the
                # next step
                # ====
                
                # == dividing each row of a TbTMatrix by its geometric mean, and taking log ==
                # NB: because the ratios are already "logged" it is just subtracting the rowMeans, see NB2
                # NB2: remember (see keep Note: #statistics #work geometric mean): geometric_mean(x) with x = x1, ... xn = (x1*...*xn)^1/n = exp(mean(log(x))), i.e. exp((1/n)*(log(x1)+ ... + log(xn))). 
                # Therefore, log(geometric_mean(x)) = (1/n)*(log(x1)+ ... + log(xn))
                # Therefore: the rowMeans of the "logged" ratios are log(geometric_mean), and thus:
                # log(Ratio/geometric_mean) = log(Ratio) - rowMean!!
                # NB3: all Inf, -Inf, and NaN are set to NA, thus all 0 in the matrixes are from ratios x/y where x = y and x AND y != 0.
                TbTMatrixesGMLog <- lapply(TbTMatrixes, function(mat) {
                        mat[!is.finite(mat)] <- NA # puts all Inf, -Inf, and also NaN to NA!
                        newM <- mat-rowMeans(mat, na.rm = TRUE)
                        #newM[is.na(newM)] <- 0
                })
                names(TbTMatrixesGMLog) <- rownames(TbTMatrixesGMLog[[1]])
                TbTMatrixesGMLog
                # ====
        })
        
}


####################################
## resultsTbT: 
###################################
# the function calculates the sums for the two sample groups in the TbTMatrixes based on classlabel, and then
# adds count info from simlist ignoring zero counts
## Input: 
# - TbTMatrixes_list: The list with the list of TbTMatrixes for each simulation
# - simlist: The list of physeq objects that have been used to calculate the TBTMatrixes_list
# - classlabel: the name of the column in the simlist physeq objects with the factor that defines the two sample groups
#  usually in your simulations = sample_data(simlist[[1]])$postfix
## Output: 
# - list of result DF for each simulation ordered after TbT_groupSum
#   sumGrp1 and sumGrp2 can be used to see whether 
#   the Taxa is more abundant in grp1 (positive) or in grp2.

resultsTbT <- function(TbTMatrixes_list, simlist, classlabel = "postfix") {
        if(!identical(names(simlist),names(TbTMatrixes_list))){"names of simlist and TbTMatrixes_list do not fit"}
        
        lapply(names(TbTMatrixes_list), FUN = function(simName){
                Mat_list <- TbTMatrixes_list[[simName]]
                sim <- simlist[[simName]]
                groupFactor <- sample_data(sim)[[classlabel]]
                sumGrp1 <- sapply(Mat_list, function(mat){sum(mat[,grepl(levels(groupFactor)[1], colnames(mat))], na.rm = T)})
                sumGrp2 <- sapply(Mat_list, function(mat){sum(mat[,grepl(levels(groupFactor)[2], colnames(mat))], na.rm = T)})
                sumAll <- sapply(Mat_list, function(mat){sum(mat, na.rm = T)})
                DF <- data.frame(id = names(sumAll), sumGrp1 = sumGrp1, sumGrp2 = sumGrp2, sumAll = sumAll, TbT_groupSum = pmax(sumGrp1, sumGrp2))
                CT <- as(otu_table(sim), "matrix")
                CT[CT == 0] <- NA
                DF$Median_grp1 <- apply(CT[grepl(levels(groupFactor)[1], rownames(CT)),], 2, median, na.rm = T)
                DF$Median_grp2 <- apply(CT[grepl(levels(groupFactor)[2], rownames(CT)),], 2, median, na.rm = T)
                DF$Mean_grp1 <- apply(CT[grepl(levels(groupFactor)[1], rownames(CT)),], 2, mean, na.rm = T)
                DF$Mean_grp2 <- apply(CT[grepl(levels(groupFactor)[2], rownames(CT)),], 2, mean, na.rm = T)
                DF$Present_grp1 <- apply(CT[grepl(levels(groupFactor)[1], rownames(CT)),], 2, function(cnts){sum(!is.na(cnts))})
                DF$Present_grp2 <- apply(CT[grepl(levels(groupFactor)[2], rownames(CT)),], 2, function(cnts){sum(!is.na(cnts))})
                DF$Zeros_grp1 <- apply(CT[grepl(levels(groupFactor)[1], rownames(CT)),], 2, function(cnts){sum(is.na(cnts))})
                DF$Zeros_grp2 <- apply(CT[grepl(levels(groupFactor)[2], rownames(CT)),], 2, function(cnts){sum(is.na(cnts))}) 
                DF$Sparsity_grp1 <- 100*(DF$Zeros_grp1/sum(grepl(levels(groupFactor)[1], rownames(CT))))
                DF$Sparsity_grp2 <- 100*(DF$Zeros_grp2/sum(grepl(levels(groupFactor)[2], rownames(CT))))
                # -- add fisher exact test of presence differences (should be none in simulation) --
                Fisher <- t(sapply(1:ntaxa(sim), FUN = function(i){
                        fisherMat <- matrix(c(DF$Present_grp1[i], DF$Zeros_grp1[i], DF$Present_grp2[i],
                                              DF$Zeros_grp2[i]), ncol = 2, dimnames = list(c("Present", "Absent"), c("grp1", "grp2")))
                        Test <- fisher.test(fisherMat)
                        cbind(Test$p.value, Test$estimate)
                }))
                fdrFisher <- p.adjust(Fisher[,1], method = "fdr")
                DF$pFisher <- Fisher[,1]
                DF$fdrFisher <- fdrFisher
                DF$oddsRatioFisher <- Fisher[,2]
                DF[order(DF$TbT_groupSum, decreasing = T),]
        })
}



####################################
## resultsTbTwilcoxTest: 
###################################
# An attempt to analyse TbTMatrixes with wilcox.test, similar to wilcoxTestApply. Basically simply applying wilcox.test on the colSums or colMeans
# of the TbTMatrixes, i.e. the sum or the mean over all taxa for each sample
## Input: 
# - TbTMatrixes_list: The list with the list of TbTMatrixes for each simulation
# - simlist: The list of physeq objects that have been used to calculate the TBTMatrixes_list
# - classlabel: the name of the column in the simlist physeq objects with the factor that defines the two sample groups
#  usually in your simulations = sample_data(simlist[[1]])$postfix
# - colSums: if TRUE colSums are used, if FALSE colMeans
# - excludeZeros: as in wilcoxTestApply, probably here should always be TRUE, samples where the host taxon was not present are excluded from the wilcox.test.
## Output: 
# - list of result DF for each simulation ordered after the wilcox.test standardized test statistic
# - also here a fisher exact test is added to compare sparsity levels.

resultsTbTwilcoxTest <- function(TbTMatrixes_list, simlist, classlabel = "postfix", colSums = TRUE, excludeZeros = TRUE) {
        if(!identical(names(simlist),names(TbTMatrixes_list))){"names of simlist and TbTMatrixes_list do not fit"}
        lapply(names(TbTMatrixes_list), FUN = function(simName){
                Mat_list <- TbTMatrixes_list[[simName]]
                sim <- simlist[[simName]]
                groupFactor <- sample_data(sim)[[classlabel]]
                if(colSums){
                        measureMatrix <- sapply(Mat_list, colSums, na.rm = T)
                        # NB: samples in which host taxon was 0 are all NA in the TbTMatrix, colSums of all NA samples = 0!
                } else {
                        measureMatrix <- sapply(Mat_list, colMeans, na.rm = T)
                        measureMatrix[is.na(measureMatrix)] <- 0
                }
                
                res_mat <- apply(measureMatrix, 2, function(taxon_col_sum){
                        x <- taxon_col_sum[grepl(levels(groupFactor)[1], names(taxon_col_sum))]
                        Zeros_grp1 <- sum(x == 0) #NB: this assumes it never happens that a colSum is exactly 0 unless the host taxon was 0 and thus all ratios NA
                        Sparsity_grp1 <- 100*(Zeros_grp1/length(x))
                        Present_grp1 <- length(x)-Zeros_grp1
                        if(excludeZeros){
                                # in case there are only 0s in x
                                if(all(x == 0)){x[1] <- ceiling(mean(taxon_col_sum))+1} 
                                x <- x[x != 0]
                        }
                        Median_grp1 <- median(x, na.rm = T)
                        Mean_grp1 <- mean(x, na.rm = T)
                        y <- taxon_col_sum[grepl(levels(groupFactor)[2], names(taxon_col_sum))]
                        Zeros_grp2 <- sum(y == 0)
                        Present_grp2 <- length(y)-Zeros_grp2
                        Sparsity_grp2 <- 100*(Zeros_grp2/length(y))
                        if(excludeZeros){
                                if(all(y == 0)){y[1] <- ceiling(mean(taxon_col_sum))+1}
                                y <- y[y != 0]
                        }
                        Median_grp2 <- median(y, na.rm = T)
                        Mean_grp2 <- mean(y, na.rm = T)
                        wilcTest <- wilcox.test(x = x, y = y, alternative = "two", paired = F, exact = F)
                        pValue <- wilcTest$p.value
                        W <- wilcTest$statistic
                        # calculate standardized rank sum Wilcoxon statistics as in multtest::mt.minP
                        Ranks <- rank(c(x, y))
                        n1 <- length(x)
                        n2 <- length(y)
                        # Wx <- sum(Ranks[1:n1])-(n1*(n1+1)/2) # would be the same as W
                        standStat <- -1*((sum(Ranks[1:n1]) - n1*(n1+n2+1)/2)/sqrt(n1*n2*(n1+n2+1)/12))
                        # check that multtest::mt.minP would give the same statistic
                        # mati <- matrix(c(x,y), nrow = 1)
                        # grFac <- as.factor(c(rep("a",8), rep("b",8)))
                        # testStat <- multtest::mt.minP(mati, grFac, test = "wilcoxon")$teststat
                        # identical(standStat, testStat) # TRUE
                        
                        # -- add fisher exact test of presence differences (should be none in simulation) --
                        fisherMat <- matrix(c(Present_grp1, Zeros_grp1, Present_grp2, Zeros_grp2), 
                                            ncol = 2, dimnames = list(c("Present", "Absent"), c("grp1", "grp2")) )
                        Test <- fisher.test(fisherMat)
                        
                        c(teststat = standStat, rawp = pValue, Median_grp1 = Median_grp1, Median_grp2 = Median_grp2, 
                          Mean_grp1 = Mean_grp1, Mean_grp2 = Mean_grp2, n1 = n1, n2 = n2, Present_grp1 = Present_grp1, 
                          Present_grp2 = Present_grp2, Zeros_grp1 = Zeros_grp1, Zeros_grp2 = Zeros_grp2, 
                          Sparsity_grp1 = Sparsity_grp1, Sparsity_grp2 = Sparsity_grp2, W, 
                          pFisher = Test$p.value, Test$estimate)
                })
                
                res_mat <- t(res_mat)
                fdr <- p.adjust(res_mat[,"rawp"], method = "fdr")
                fdrFisher <- p.adjust(res_mat[,"pFisher"], method = "fdr")
                DF <- data.frame(id = rownames(res_mat), res_mat, adjp = fdr, fdr = fdr, fdrFisher = fdrFisher)
                DF <- dplyr::select(DF, 1:3, 19:20, 4:7, 10:15, 8:9, 16, 17, 21, 18)
                colnames(DF)[21] <- "oddsRatioFisher"
                DF <- dplyr::arrange(DF, desc(abs(teststat)))
        })
}


####################################
## pValueMatrixesTbTplusTile: 
###################################
# the function calculates ntaxa x ntaxa matrixes of p-values for each simulation.
# Based on the TbTMatrixes (probably ratios better than gm log ones) the p-value matrix
# indicates which taxa is enriched compared to which other taxa.
# in addition it adds a tile plot for each pValueMatrix
# wilcoxon.test is used
## Input: 
# - TbTMatrixes_list: The list with the list of TbTMatrixes for each simulation
# - simlist: The list of physeq objects that have been used to calculate the TBTMatrixes_list
# - classlabel: the name of the column in the simlist physeq objects with the factor that defines the two sample groups
#  usually in your simulations = sample_data(simlist[[1]])$postfix
## Output: 
# - list of pValMatrixes plus TileTr for each simulation. So for each simulation a list again with the matrix and the trellis.
# NB: I negated the p-values if host taxon was more abundant in grp2 compared to other taxon!!

pValueMatrixesTbTplusTile <- function(TbTMatrixes_list, simlist, classlabel = "postfix") {
        if(!identical(names(simlist),names(TbTMatrixes_list))){"names of simlist and TbTMatrixes_list do not fit"}
        pValueMat_Tr_List <- lapply(names(TbTMatrixes_list), FUN = function(simName){
                
                Mat_list <- TbTMatrixes_list[[simName]]
                sim <- simlist[[simName]]
                groupFactor <- sample_data(sim)[[classlabel]]
                
                pValMatrix <- sapply(Mat_list, function(mat){
                        apply(mat, 1, function(taxon_ratios){
                                x <- taxon_ratios[grepl(levels(groupFactor)[1], names(taxon_ratios))]
                                # NB: it is possible that x is completely NA (if taxon was only present when host taxon was not!), therefore
                                if(sum(!is.na(x)) == 0){x[1] <- 0}
                                y <- taxon_ratios[grepl(levels(groupFactor)[2], names(taxon_ratios))]
                                if(sum(!is.na(y)) == 0){y[1] <- 0}
                                pValue <- wilcox.test(x = x, y = y, alternative = "two", paired = F, exact = F)$p.value
                                # change sign of pValue, so the value is positive for a taxon in its row when it is more abundant in group1
                                Ranks <- rank(c(x[!is.na(x)], y[!is.na(y)]))
                                n1 <- length(x[!is.na(x)])
                                n2 <- length(y[!is.na(y)])
                                Wx <- sum(Ranks[1:n1])-(n1*(n1+1)/2)
                                Wy <- sum(Ranks[(n1+1):(n1+n2)])-(n2*(n2+1)/2)
                                if(Wx > Wy){pValue <- -1*pValue}
                                pValue
                        })
                })
                
                # -- add a tile plot of the pValMatrix --
                
                DF <- as.data.frame(pValMatrix)
                DF[is.na(DF)] <- 1
                DF$HostTaxon <- rownames(pValMatrix)
                DF <- tidyr::gather(DF, key = Taxon , value = pValue, - HostTaxon)
                DF$Taxon <- as.factor(DF$Taxon)
                DF$HostTaxon <- as.factor(DF$HostTaxon)
                LW <- rownames(pValMatrix)
                LW2 <- rev(LW)
                for(z in 1:nrow(pValMatrix)){
                        DF$Taxon <- relevel(DF$Taxon, ref = LW2[z])
                        DF$HostTaxon <- relevel(DF$HostTaxon, ref = LW[z])
                }
                # add color to the taxa names so you see up and down
                colyaxis <- vector(mode = "character", length = nrow(DF))
                colyaxis[] <- "black"
                colyaxis[grepl("TP-U", levels(DF$HostTaxon))] <- "#E69F00"
                colyaxis[grepl("TP-D", levels(DF$HostTaxon))] <- "#009E73"
                colxaxis <- vector(mode = "character", length = nrow(DF))
                colxaxis[] <- "black"
                colxaxis[grepl("TP-U", levels(DF$Taxon))] <- "#E69F00"
                colxaxis[grepl("TP-D", levels(DF$Taxon))] <- "#009E73"
                DF$Fill <- "not significant"
                DF$Fill[DF$pValue<0.05 & DF$pValue> 0] <- "up (p < 0.05)"
                DF$Fill[DF$pValue>-0.05 & DF$pValue< 0] <- "down (p < 0.05)"
                DF$Fill <- factor(DF$Fill, levels = c("up (p < 0.05)", "not significant", "down (p < 0.05)"))
                TileTr <- ggplot(DF, aes(x=Taxon, y=HostTaxon, fill=Fill))
                TileTr <- TileTr + 
                        geom_raster() + 
                        ggtitle(simName) +
                        scale_fill_manual(values = c("not significant" = "gray98", "up (p < 0.05)" = "#E69F00", "down (p < 0.05)" = "#009E73")) +
                        scale_x_discrete(position = "top") +
                        labs(x=NULL, y=NULL) +
                        theme_tufte(base_family="Helvetica") +
                        theme(plot.title=element_text(hjust=0)) +
                        theme(axis.ticks=element_blank()) +
                        theme(axis.text=element_text(size=7)) +
                        theme(legend.title=element_blank()) +
                        theme(legend.text=element_text(size=6)) +
                        theme(axis.text.y = element_text(colour = colyaxis),
                              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0, colour = colxaxis))
                list(pValMatrix = pValMatrix, TileTr = TileTr)
        })
        
        names(pValueMat_Tr_List) <- names(simlist)
        pValueMat_Tr_List
}






# ------------------------- deprecated function ---------------------------------


# ####################################
# ## CalcTbTMatrixesOLD:
# ###################################
# 
# ## Input: 
# # - simlist = a list with physeq objects (usually the simulations) with count data
# # - Pseudocount: by default = 0.5. Each 0 in the count data of a simulation is
# #   replaced by this Pseudocount. All other counts remain untouched. Necessary since
# #   0 would clash with ratios
# ## Output: 
# # - list (length = number of simulations in simlist) in which each entry 
# #   is again a list (length = number of OTUs in the simulation) with the 
# #   normalised and log transformed ObORatioMatrixes for each host OTU 
# 
# 
# CalcTbTMatrixesOLD = function(simlist, Pseudocount = .5){
#         TbTMatrixesList <- list()
#         for (i in 1:length(simlist)) {
#                 
#                 physeq <- simlist[[i]]
#                 CM <- t(as(otu_table(physeq), 'matrix'))
#                 
#                 # 1.) add Pseudocount ---------------------
#                 # (unique(sort(CM)))[2] # would give you the second smalles count
#                 CM[CM==0] <- Pseudocount
#                 
#                 # 2.) calculate TbyT Matrix for each Taxum ---------------------
#                 TbTMatrixes <- lapply(1:nrow(CM), function(i){apply(CM, 2, function(samp_cnts){samp_cnts[i]/samp_cnts})})
# 
#                 
#                 # 3.) Divide the TbTMatrixes by the geometric means over all samples and take the log: -------
#                 TbTMatrixesGMLog <- lapply(TbTMatrixes, function(x) {log(x)-rowMeans(log(x))})
#                 
#                 # # This was the short version of this commented longer version:
#                 # 
#                 # # 3.) Divide the TbTMatrixes by the geometric means over all samples.
#                 # 
#                 # TbTMatrixesGM <- lapply(1:length(TbTMatrixes), function(x) {
#                 #         GM <- apply(TbTMatrixes[[x]], 1, gm_own, zeros.count = TRUE) # because of the pseudocounts anyway no zeros
#                 #         Mat <- TbTMatrixes[[x]]/GM} )
#                 # 
#                 # # so NB that the products across the rows of these matrixes, i.e. over all samples is 1, e.g.
#                 # # apply(TbTMatrixesGM[[2]], 1, prod)
#                 # 
#                 # # 4.) take the log ---------------------
#                 # 
#                 # TbTMatrixesGMLog <- lapply(TbTMatrixesGM, log)
#                 
#                 # apply(TbTMatrixesGMLog[[2]], 1, sum) #practically 0 
# 
#                 names(TbTMatrixesGMLog) <- rownames(TbTMatrixesGMLog[[1]])
#                 TbTMatrixesList[[i]] <- TbTMatrixesGMLog
#         }
#         names(TbTMatrixesList) <- names(simlist)
#         TbTMatrixesList     
# }

      




