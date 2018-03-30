#######################################
### simulate_counts_withTPTaxa##
#################
## simulates a count table from a given template physeq where a defined number of true positive taxa will
# be more abundant either in grp1 or in grp2. The degree of abundance enrichment is defined by foldeffects. 
## requires: phyloseq, uses simulate_count_samples (NB for all physeqs generated taxa_are_rows = FALSE)
## Input: 
# - templatelist: a list with physeq template objects named with sampletypes (I usually use only one sampletype, 
# so only one template in a list)
# - simparams: character vector with the simulation parameters to use named in simparamslabels, specifically
#     - n: the median number of liberary sizes (TotalAmplicons) used to sample from for the simulated samples
#     - min and maxEffectSize: the fold effects will be sampled from this range
#     - J1 and J2: the number of samples that will be generated for grp1 and grp2 respectively
# - nTP: the number of taxa that will be true positives
# - NB: uses more inputs that should be defined before in the code
# 
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
                
                # - sample LibrarySizes (TotalAmplicons) for the simulated samples in each group -
                scaledSums = round(n*(sample_sums(templatei)/median(sample_sums(templatei))))
                # in case of n = median(sample_sums(physeqi)) = default, scaledSums = sample_sums(templatei)
                n1   = sample(scaledSums, size = Ji1, replace = TRUE)
                n2   = sample(scaledSums, size = Ji2, replace = TRUE)
                # --
                
                # -  Set the TP taxa -
                TPTaxa = sample(taxa_names(templatei), nTP, replace=FALSE) # NB: will throw an error if nTP bigger ntaxa(templatei)
                # randomly sample how many TruePositives will be more abundant in grp1:
                upGr1 <- sample(nTP, 1)
                # set which TP will be more abundant in Gr1 and which less abundant in Grp1 (more abundant in grp2)
                TPTaxa1 <- sample(TPTaxa, upGr1)
                TPTaxa2 <- setdiff(TPTaxa, TPTaxa1)
                
                wh.TP1 = taxa_names(templatei) %in% TPTaxa1
                wh.TP2 = taxa_names(templatei) %in% TPTaxa2
                
                # sample the foldeffects for the different TP taxa so you can put the FE into the TP taxa names
                # NB: I need if statement to check whether minEffectSize < maxEffectSize because sample samples from 1:x if x 
                # is a single numeric, and thus 3:3 = 3 would produce values between 1 and 3 instead of all 3!!
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
                # --
                
                # -  Simulate the data separately for grp1 and grp2 using simulate_count_samples function -
                sim1 = simulate_count_samples(postfix = paste0(sampletypei, ";grp1"), physeq = templatei, NoSamples = Ji1, SampleSizes = n1,
                                              TPNames = newname1, FEffects = effectsizes1)
                sim2 = simulate_count_samples(postfix = paste0(sampletypei, ";grp2"), physeq = templatei, NoSamples = Ji2, SampleSizes = n2,
                                              TPNames = newname2, FEffects = effectsizes2)
                # --
                
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
# NoSamples samples will have the same length, or n of length NoSamples to define length of each
# sample) 
# TPNames: the names of the taxa that will be upregulated (over abundant)
# FEffects: the TPTaxa will be upregulated by the given FEffects
## Output: a phyloseq object with the simulated count table with some very simple sample data
# to understand the simulation parameters, and the tax_table from physeq

simulate_count_samples = function(postfix="sim", physeq, NoSamples, SampleSizes, TPNames, FEffects){
        
        # - generate the urn matrix to draw the counts from later, including the 0 for the sparsity -
        pi = taxa_sums(physeq)
        
        # add the effect at the TPTaxa
        pi[which(names(pi) %in% TPNames)] <- pi[which(names(pi) %in% TPNames)]*FEffects
        
        # new (22.04.2017): correct pi with NoSamples/nsamples(physeq) 
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
        
        # changing pi into a matrix to deal with sparsity problem 
        # New 23.04.2017: Remember old version produced data with very low sparsity
        piMat <- matrix(rep(pi, NoSamples), ncol = NoSamples)
        # each column will be used as an urn
        # Determine the sparsity of each taxon 
        if(taxa_are_rows(physeq)) { physeq <- t(physeq) }
        TaxaSpars <- colSums(otu_table(physeq) == 0)/nrow(otu_table(physeq))
        NoZerosPerTaxa <- floor(TaxaSpars*NoSamples)
        # put zeros randomly into the pimatrix
        ZeroIndexList <- lapply(NoZerosPerTaxa, FUN = function(x){sample(NoSamples,x)}) # tells you the Sample IDs in which the Taxon will be zero
        piMat <- t(sapply(1:nrow(piMat), function(e){piMat[e,ZeroIndexList[[e]]] = 0;
        piMat[e,]}))
        # NB: plot(colSums(piMat)) you see that the colsums now differ, so a taxon can be more likely to become abundant in one sample than in another,
        # But note, it is random, and will be the same for the samples of the other group, so it should not significantly affect the TP status of a taxon

        # the next line is for the special case that NoSamples is = 1
        if(ncol(piMat) != NoSamples && nrow(piMat) == NoSamples){piMat <- t(piMat)}
        if(ncol(piMat) != NoSamples){ stop("NoSamples unequal to the number of columns in piMat") }
        # result: piMat with as many zeros (roughly) as the template physeq
        # --
        
        # - pick the counts based on piMat using sparsity_subsample -  
        simct = sapply(1:NoSamples, function(i){sparsity_subsample(piMat, SampleSizes, i)})
        # --
        
        # - transform count matrix into physeq object -
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
# ProbMatrix: columns define the probabilities with which the sample will be generated
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


