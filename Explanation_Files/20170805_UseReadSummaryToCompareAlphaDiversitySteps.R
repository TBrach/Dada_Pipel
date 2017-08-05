samdf <- sample_data(ps)
samdf <- as(samdf, "data.frame")
samdf$Sample <- rownames(samdf)

RS1 <- merge(ReadSummary, samdf, by = "Sample")


# Filtered Reads the ones that the dereplication and denoising is done on
Tr <- ggplot(RS1, aes(x = Group, y = FilteredReads, col = Group))
Tr <- Tr + geom_boxplot() +
        geom_jitter()
Tr
pairwise.t.test(x = RS1$FilteredReads, g = RS1[[group_var]], alternative = "two", p.adjust.method = "none", var.equal = F, pool.sd = F)


# Richness
Tr <- ggplot(RS1, aes(x = Group, y = UniqueAmpliconsWOBimera, col = Group))
Tr <- Tr + geom_boxplot() +
        geom_jitter()
Tr
pairwise.t.test(x = RS1$UniqueAmpliconsWOBimera, g = RS1[[group_var]], alternative = "two", p.adjust.method = "none", var.equal = F, pool.sd = F)


# Richness in seqtab (i.e. without bimera removal)
Tr <- ggplot(RS1, aes(x = Group, y = Unique_Amplicons, col = Group))
Tr <- Tr + geom_boxplot() +
        geom_jitter()
Tr
pairwise.t.test(x = RS1$Unique_Amplicons, g = RS1[[group_var]], alternative = "two", p.adjust.method = "none", var.equal = F, pool.sd = F)

# so quite some effect comes from bimera (could try other bimera removal options)

# Already more unique sequences_F ?
Tr <- ggplot(RS1, aes(x = Group, y = UniqueSequences_F, col = Group))
Tr <- Tr + geom_boxplot() +
        geom_jitter()
Tr
pairwise.t.test(x = RS1$UniqueSequences_F, g = RS1[[group_var]], alternative = "two", p.adjust.method = "none", var.equal = F, pool.sd = F)


# Already more unique sequences_R ?
Tr <- ggplot(RS1, aes(x = Group, y = UniqueSequences_R, col = Group))
Tr <- Tr + geom_boxplot() +
        geom_jitter()
Tr
pairwise.t.test(x = RS1$UniqueSequences_R, g = RS1[[group_var]], alternative = "two", p.adjust.method = "none", var.equal = F, pool.sd = F)

# so at least the old samples got more unique sequences out of the same number of filtered reads

# Already more denoised sequences_F ?
Tr <- ggplot(RS1, aes(x = Group, y = DenoisedSequences_F, col = Group))
Tr <- Tr + geom_boxplot() +
        geom_jitter()
Tr
pairwise.t.test(x = RS1$DenoisedSequences_F, g = RS1[[group_var]], alternative = "two", p.adjust.method = "none", var.equal = F, pool.sd = F)


# Already more denoised sequences_R ?
Tr <- ggplot(RS1, aes(x = Group, y = DenoisedSequences_R, col = Group))
Tr <- Tr + geom_boxplot() +
        geom_jitter()
Tr
pairwise.t.test(x = RS1$DenoisedSequences_R, g = RS1[[group_var]], alternative = "two", p.adjust.method = "none", var.equal = F, pool.sd = F)


# but more sequences are actually considered error, in the denoised sequences there is little difference. But then more merged, and less removed by bimeras

# More bimeras anywhere

RS1$RemovedSVs <- RS1$Unique_Amplicons - RS1$UniqueAmpliconsWOBimera

# Already more denoised sequences_R ?
Tr <- ggplot(RS1, aes(x = Group, y = RemovedSVs, col = Group))
Tr <- Tr + geom_boxplot() +
        geom_jitter()
Tr
pairwise.t.test(x = RS1$RemovedSVs, g = RS1[[group_var]], alternative = "two", p.adjust.method = "none", var.equal = F, pool.sd = F)



