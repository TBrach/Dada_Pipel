# THE NEWER ATTMEPT (also see older below almost the same)

# get mergers from going throught the DanFunD data
Df <- mergers[[1]]
ddF <- dd_F[[1]]
drpF <- drp_F[[1]]
ddR <- dd_R[[1]]
drpR <- drp_R[[1]]

# Df tells you that denoised F 1 fitted to denoised R 1 with abundance 1640.
# To get to the abundance you simply go backwards
# First from denoised to uniques using ddF$map and ddR$map
UniquesFittingDenoisedF1 <- which(ddF$map == 1)
UniquesFittingDenoisedR1 <- which(ddR$map == 1)
# drpF$uniques[UniquesFittingDenoisedF1] # gives you the unique sequences

# Second: you find all the reads that in the end were assigned to the denoised sequence using drpF$map, and drpR$map 

ReadsFittingDenoisedF1 <- which(drpF$map %in% UniquesFittingDenoisedF1)
ReadsFittingDenoisedR1 <- which(drpR$map %in% UniquesFittingDenoisedR1)

# NB:
length(ReadsFittingDenoisedF1) >= length(UniquesFittingDenoisedF1) 

# NOW to get the abundance, you realise it gets crucial that FW and RV reads are in the same order (same ID),
# you know the FW reads that were assigned to the denoised FW sequence, same for the reverse,
# you want to know the intersection of these IDs.
length(intersect(ReadsFittingDenoisedF1, ReadsFittingDenoisedR1)) # 1640

# NB: THIS explains why a denoised RV sequence can be merged with several FW sequences!
# the abundances are just counted for the read IDs that actually fit together (from same cluster)

# THE OLDER ATTEMPT

derepF = derepFastq(system.file("extdata", "sam1F.fastq.gz", package="dada2"))
derepR = derepFastq(system.file("extdata", "sam1R.fastq.gz", package="dada2"))
dadaF <- dada(derepF, err=tperr1, errorEstimationFunction=loessErrfun, selfConsist=TRUE)
dadaR <- dada(derepR, err=tperr1, errorEstimationFunction=loessErrfun, selfConsist=TRUE)
merger <- mergePairs(dadaF, derepF, dadaR, derepR)
merger2 <- mergePairs(dadaF, derepF, dadaR, derepR, returnRejects = T)

# merger tells you that denoised F 1 fitted to denoised R 1, how to get to the abundance of 584?

UniquesFittingDenoisedF1 <- which(dadaF$map == 1)
UniquesFittingDenoisedR1 <- which(dadaR$map == 1)

ReadsFittingUniquesDenoisedF1 <- which(derepF$map %in% UniquesFittingDenoisedF1)
ReadsFittingUniquesDenoisedR1 <- which(derepR$map %in% UniquesFittingDenoisedR1)

length(intersect(ReadsFittingUniquesDenoisedF1, ReadsFittingUniquesDenoisedR1))
# exactly the 584

# this explains why you cannot get from Denoised_F and bimera_F to Unique Amplicons

F_DenoisedWOChimera <- ReadSummary$Denoised_F - ReadSummary$bimera_F
R_DenoisedWOChimera <- ReadSummary$Denoised_R - ReadSummary$bimera_R


## Name suggestions:
AllReads
FilteredReads
MergedReads
MergedReadsWOBimera
UniqueSequences_F
DenoisedSequences_F
BimeraSequences_F
UniqueAmplicons
UniqueAmpliconsWOBimera/Richness/Observed