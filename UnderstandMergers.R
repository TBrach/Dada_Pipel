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