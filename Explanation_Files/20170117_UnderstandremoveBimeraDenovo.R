# you have a seqtab (on which no bimera removal has been taken place)

seqtab.nochim.Pooled <- removeBimeraDenovo(seqtab, method="pooled", verbose = T)
seqtab.nochim.PerSample <- removeBimeraDenovo(seqtab, method="per-sample", verbose = T)
seqtab.nochim.Consensus <- removeBimeraDenovo(seqtab, method="consensus", verbose = T)

NoBimerasPooled <- ncol(seqtab) - ncol(seqtab.nochim.Pooled)
NoBimerasPerSample <- ncol(seqtab) - ncol(seqtab.nochim.PerSample)
NoBimerasConsensus <- ncol(seqtab) - ncol(seqtab.nochim.Consensus)

# -- "pooled": just run isBimeraDenovo on colSums(seqtab) --

pooledSeqtab <- matrix(colSums(seqtab), ncol = ncol(seqtab), dimnames = list("", colnames(seqtab)))

identical(sum(isBimeraDenovo(pooledSeqtab)), NoBimerasPooled)
# True

# ----

# -- "per-sample": run isBimeraDenovo on each sample, i.e. row in seqtab --

mat <- t(apply(seqtab, 1, isBimeraDenovo))
# weirdly it also checks for sequences with abundance = 0
seqtab.PS <- seqtab
seqtab.PS[mat] <- 0
identical(NoBimerasPerSample, sum(colSums(seqtab.PS) == 0)) # TRUE

# ----

# -- "consensus": run isBimeraDenovo on each sample, then check if bimera in "majority, default = 0.9" of samples leaving one sample out
# see: isBimeraDenovoTable, real code is simply bim <- isBimeraDenovoTable(seqtab)

# exclude 
mat[seqtab == 0] <- NA
NumberOfTrues <- colSums(mat, na.rm = T)
NumberOfFalses <- colSums(!mat, na.rm = T)
# see ignoreNNegatives in isBimeraDenovo is by default 1, so reduce Number of FALSES by 1
NumberOfFalses[NumberOfFalses > 0] <- NumberOfFalses[NumberOfFalses > 0] - 1
Total <- NumberOfTrues + NumberOfFalses
Fraction <- NumberOfTrues/Total
identical(NoBimerasConsensus, sum(Fraction >= 0.9, na.rm = T)) # sometimes TRUE

# --



# # extra fiddling
# NoBimerasConsensus - NoBimerasPerSample
# ExtraOutByConsensus <- which(colnames(seqtab) %in% colnames(seqtab.nochim.PerSample)[!(colnames(seqtab.nochim.PerSample) %in% colnames(seqtab.nochim.Consensus))])
# MatOuts <- mat[,ExtraOutByConsensus]
# MatOuts2 <- matB[,ExtraOutByConsensus]
# SeqOuts <- seqtab[,ExtraOutByConsensus]

