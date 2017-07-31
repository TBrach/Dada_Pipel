# from: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4955027/
# on page 8

savepath <- "/Users/jvb740/MarieCurie_Work/MouseProject/ResultsAndProtocols/ManiAging_Results/16S_Sequencing/2017-07-13_DK_age_ManiAging/Dada2_Analysis/Dada_phylogenetic_tree"

source("https://bioconductor.org/biocLite.R")
biocLite("DECIPHER")
install.packages("phangorn")

library(DECIPHER); packageVersion("DECIPHER")
library(phangorn); packageVersion("phangorn")

seqs <- getSequences(seqtab.nochim) # simple character string
names(seqs) <- seqs

# ---- multiple alignment of SVs using "DECIPHER" (takes only 10 seconds) ----

alignment <- DECIPHER::AlignSeqs(DNAStringSet(seqs), anchor = NA)

# --------

# ---- use alignment to construct phylogenetic tree (phangorn package) -----

# -- First construct a neighbor-joining tree --

phang.align <- phangorn::phyDat(as(alignment, "matrix"), type = "DNA")
dm <- phangorn::dist.ml(phang.align) # computes pairwise distances

treeNJ <- phangorn::NJ(dm)

# NB: They warn: tip order != sequence order
# however, I got a TRUE here
# all.equal(treeNJ$tip.label, seqs, check.attributes = FALSE)

fit <- phangorn::pml(treeNJ, data = phang.align) # warns: negative edges length changed to 0!

# ----

# -- use the neighbor-joining tree as a start point to fit a GTR+G+I (generalized time-reversible with Gamma rate variation) maximum 
# likelihood tree --

fitGTR <- update(fit, k = 4, inv = 0.2)
fitGTR <- phangorn::optim.pml(fitGTR, model = "GTR", optInv = TRUE, optGamma = TRUE,
                              rearrangement = "stochastic", control = pml.control(trace = 0))

# ----

# --------





