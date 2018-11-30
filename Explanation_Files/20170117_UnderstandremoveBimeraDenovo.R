# Some kind of summary here:
# - they claim that thanks to the precision of ASVs it is much easier to detect bimera based on ASV than with OTUs. Note however, it seems that they
# completely ignore multimera and restrict themselves to bimera
# - for a short description of the algorithm see help isBimeraDenovo and later isBimeraDenovoTable
# - I describe here the differences between Pooled PerSample and Consensus based on an example seqtab with 6 samples and 807 ASV but any seqtab should do
# my conclusion is:
#       - I think pooled is a bit weird because it removes most ASV and can potentially remove ASVs even though the parents were not present in that sample
#       - so I vote for per-sample or consensus which is per sample but when an ASV is considered Bimera in a sufficient number of samples it is considered
#       - bimera in all samples. I think that makes sense and it is also the default. 



# - load the seqtab (NB: if this path does not exist find whatever seqtab) -
datapath <- "/Users/jvb740/MarieCurie_Work/Project_Normalization/DK_Healthy_Normalization_20181023_Reads/Dada_Analysis_filteredReads42500/Dada_Data"
load(file.path(datapath, "DenoisedData.RData"))
# --
dim(seqtab) # six samples rows, and 807 ASV
# the interesting thing about this seqtab is that Samples 1 3 4 and 2 5 6 are from the same Person (even stool sample)
apply(seqtab, 1, function(row){sum(row > 0)}) # not a single sample has more than 300 of these 807 ASV that are found in all smaples
min(seqtab[seqtab > 0]) # yes there are in fact two singletons in this dataset to my surprise
sum(seqtab == 1)

# remove bimeras with the three different methods
seqtab.nochim.Pooled <- removeBimeraDenovo(seqtab, method="pooled", verbose = T)
seqtab.nochim.PerSample <- removeBimeraDenovo(seqtab, method="per-sample", verbose = T)
seqtab.nochim.Consensus <- removeBimeraDenovo(seqtab, method="consensus", verbose = T)

NoBimerasPooled <- ncol(seqtab) - ncol(seqtab.nochim.Pooled)
NoBimerasPerSample <- ncol(seqtab) - ncol(seqtab.nochim.PerSample)
NoBimerasConsensus <- ncol(seqtab) - ncol(seqtab.nochim.Consensus)


OverviewDF <- data.frame(Method = c("pooled", "per-sample", "consensus"), NoBimeras =  c(NoBimerasPooled, NoBimerasPerSample, NoBimerasConsensus))
# so you see that pooled removes most ASV, per-sample the least and consensus is in between
# let's try to understand the methods one by one and go a bit more into detail with what is removed



# -- "pooled": just run isBimeraDenovo on colSums(seqtab) --

pooledSeqtab <- matrix(colSums(seqtab), ncol = ncol(seqtab), dimnames = list("", colnames(seqtab)))

identical(sum(isBimeraDenovo(pooledSeqtab)), NoBimerasPooled)
# TRUE
# ----



# -- "per-sample": run isBimeraDenovo on each sample, i.e. row in seqtab --
BimeraPerSample <- t(apply(seqtab, 1, isBimeraDenovo))


# it also checks for sequences with abundance = 0 and puts them all to Bimera
rowSums(BimeraPerSample)
# so to get the actual bimeras found in each sample you have to count only those that were not zero abundant
seqtabLogical <- seqtab > 0
NoBimeraPerSample <- sapply(1:6, function(i){sum(BimeraPerSample[i,] & seqtabLogical[i,])})
# quickly prove that the isBimeraDenovo result is the same for a sample if you ignored the zeros from the beginning
NoBimeraPerSample2 <- sapply(1:6, function(i) {sum(isBimeraDenovo(seqtab[i, seqtab[i,] > 0]))})
identical(NoBimeraPerSample, NoBimeraPerSample2)
# ----



# -- compare Pooled directly to PerSample to understand the situation a bit better --
# In Pooled you sum up the abundances of the ASV in all Samples and then you run isBemeraDenovo. There is the option that a low abundant ASV from a sample 
# is considered bimera based on parents that were only found in another sample. This could be considered weird. 
# on the other hand, there is also the option that an ASV is low abundant in one sample but more abundant in another and does not become bimera.
# let's compare, i.e. for each sample find overlap, ASV that did become Bimera only in pooled and only in per Sample
NoASVs <- apply(seqtab, 1, function(row){sum(row > 0)})

BimeraPooled <- isBimeraDenovo(pooledSeqtab)
NoBimeraPooled <- apply(seqtab, 1, function(row){
        Test <- row > 0
        sum(Test & BimeraPooled)
})
NoBimeraPerSample <- NoASVs - apply(seqtab.nochim.PerSample, 1, function(row){sum(row > 0)})
BimeraPerSample <- t(apply(seqtab, 1, isBimeraDenovo))
seqtabLogical <- seqtab > 0
Overlap <- sapply(1:6, function(i){sum(seqtabLogical[i,] & BimeraPerSample[i,] & BimeraPooled)})
OnlyPooled <- sapply(1:6, function(i){sum(seqtabLogical[i,] & !BimeraPerSample[i,] & BimeraPooled)})
OnlyPerSample <- sapply(1:6, function(i){sum(seqtabLogical[i,] & BimeraPerSample[i,] & !BimeraPooled)})

SumDF <- data.frame(Sample = rownames(seqtab), NoASV = NoASVs, NoBimeraPooled = NoBimeraPooled, NoBimeraPerSample = NoBimeraPerSample,
                    Overlap = Overlap, OnlyPooled = OnlyPooled, OnlyPerSample = OnlyPerSample)
# As SumDF shows you there are indeed very few ASV that are only considered Bimera in the PerSample option, but clearly more are only considered Bimera
# in the Pooled method. I clearly think that the pooled method is a bit weird considering that bimera should be mainly a PCR problem that is
# sample specific
# ----



# -- "consensus": run isBimeraDenovo on each sample, then check if bimera in "majority, default = 0.9" of samples in which it was present
# see: isBimeraDenovoTable, real code is simply bim <- isBimeraDenovoTable(seqtab)

bimeraConsensus <- isBimeraDenovoTable(seqtab)
identical(sum(bimeraConsensus), NoBimerasConsensus) # yes finds same amount of bimera

# try to recapitulate the algorithm
# first calculate the BimeraPerSample matrix
BimeraPerSample <- t(apply(seqtab, 1, isBimeraDenovo))

BimeraPerSample[seqtab == 0] <- NA
# NB: rowSums(BimeraPerSample, na.rm = T) is again NoBimeraPerSample as in SumDF
NumberOfTrues <- colSums(BimeraPerSample, na.rm = T)
NumberOfFalses <- colSums(!BimeraPerSample, na.rm = T)
# see ignoreNNegatives in isBimeraDenovo is by default 1, so reduce Number of FALSES by 1
NumberOfFalses[NumberOfFalses > 0] <- NumberOfFalses[NumberOfFalses > 0] - 1
Total <- NumberOfTrues + NumberOfFalses
Fraction <- NumberOfTrues/Total
identical(NoBimerasConsensus, sum(Fraction >= 0.9, na.rm = T)) # sometimes TRUE BUT not in my case here
# so I had to look further into this but the procedure seems at least relatively clear so I stop here for now
# NB: there are quite some NaN in the Fraction. Maybe that has an effect. 

# --


