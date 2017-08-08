# interesting points learned

# - elongated PCOAs based on what the principle coordinates explain
# - plot_tree (ok)
# - check if library size explains some axes on your ordination plots
# - double principal coordinates analysis (DPCoA) is a phylogenetic ordination method that provides a biplot representation for both samples and taxonomic categories
# - biplots from DPCoA and weighted unifrac can be used to interpret the axes

# load a phyloseq object

datapath <- "/Users/jvb740/MarieCurie_Work/MouseProject/ResultsAndProtocols/ManiAging_Results/16S_Sequencing/2017-07-13_DK_age_ManiAging/Dada2_Analysis/phyloseq_analysis"

ps <- readRDS(file.path(datapath, "physeq.rds"))

functionpath <- "/Users/jvb740/MarieCurie_Work/BackgroundKnowledge/16S_Learning/Dada_Pipel/Functions/"
source(file.path(functionpath, "Dada_TaxonomyFunctions.R"))
source(file.path(functionpath, "Dada_PlotFunctions.R"))

# -- Taxonomic Filtering --

# - remove taxa not assigned to phylum level -
# assumption: taxa that could not be assigned to phylum level are most likely sequence artifacts that don't exist in nature

table(tax_table(ps)[, "Phylum"], exclude = NULL) # without exclude NAs would not be counted

ps0 <- subset_taxa(ps, !is.na(Phylum))

# --

# - filter Phyla by prevalence - 

prevdf <- data.frame(SV_ID = 1:ntaxa(ps0), 
                         total_abundance = taxa_sums(ps0),
                         prevalence = colSums(as(otu_table(ps0), "matrix") != 0),
                         sparsity = colSums(as(otu_table(ps0), "matrix") == 0), 
                         mean_abundance_nonzero = apply(as(otu_table(ps0), "matrix"), 2, function(x){mean(x[x > 0])}),
                         median_abundance_nonzero = apply(as(otu_table(ps0), "matrix"), 2, function(x){median(x[x > 0])}))
prevdf <- cbind(prevdf, tax_table(ps0))


summarise(group_by(prevdf, Phylum), mean_prevalence = mean(prevalence), total_prevalence = sum(prevalence), n_taxa = n())

# so I would not use subset_taxa here to remove an entire phylum, so
prevdf1 <- prevdf

# --

# - Filter individual Taxa by prevalence -

# motivation: don't pay attention to taxa that are only found in very rare samples (unless you have a reason to)
# again I use my plot function here, generates the same plot/ similar plot, maybe indeed you should go to percentages

filterList <- plot_abundance_prev_filter(physeq = ps0, prevalence = 5, taxa_sums_quantile = 100)

filterList[[2]] # basically his plot


prevalenceThreshold <- 0.05 * nsamples(ps0)

keepTaxa <- rownames(prevdf1)[(prevdf1$prevalence >= prevalenceThreshold)]

ps2 <- prune_taxa(keepTaxa, ps0)

# I would say filter_taxa is smoother here

# ps2 <- filter_taxa(ps0, function(x){sum(x > 0) >= prevalenceThreshold}, prune = TRUE)

# --
# ----

# - Agglomerate taxa -

# motivation: ideally you would agglomerate based on known functional redundancy (using merge_taxa()).
# in principle you hope here that taxonomic agglomeration also reflects functional redundancy

length(get_taxa_unique(ps2, taxonomic.rank = "Genus")) # This counts NA in
ps3 <- tax_glom(ps2, "Genus", NArm = TRUE) # here all taxa with Genus = NA removed!

# <-- Do tax_glom yourself to understand it -->

# shows you that tax_glom just sums up, I did not get on which SV becomes the new name, but the name also becomes meaningless and should
# be replaced by genus name (or taxa name)

# CT <- as(otu_table(ps2), "matrix")
# TT <- tax_table(ps2)
# TT <- as.data.frame(TT)
# # kick Genus = NA out
# GenusNA <- is.na(TT[,"Genus"])
# 
# CT <- CT[,!GenusNA]
# LookUpDF <- data.frame(SV = rownames(TT)[!GenusNA], Genus = TT[!GenusNA,"Genus"])
# CT <- as.data.frame(CT)
# CT$Sample <- rownames(CT)
# CT_l <- gather(CT, key = SV, value = abundance, -Sample)
# 
# CT_l <- merge(CT_l, LookUpDF, by = "SV")
# 
# CT_l <- group_by(CT_l, Sample, Genus)
# CT_Sum <- summarise(CT_l, abundance = sum(abundance))
# CT_p3 <- tidyr::spread(CT_Sum, key = Genus, value = abundance)
# CT_p3 <- as.data.frame(CT_p3)
# rownames(CT_p3) <- CT_p3$Sample
# CT_p3 <- CT_p3[-1]
# 
# CT_realp3 <- as(otu_table(ps3), "matrix")
# LookUpDF2 <- LookUpDF[LookUpDF$SV %in% colnames(CT_realp3), ]
# # give yours the same order
# CT_p3 <- CT_p3[, match(LookUpDF2$Genus, colnames(CT_p3))]
# 
# all.equal(unname(CT_p3), unname(as.data.frame(CT_realp3))) #TRUE
# # How about the names?
# # does it always keep the one with max abundance?
# CT_l <- group_by(CT_l, Genus)
# CT_SumName <- summarise(CT_l, abundance = max(abundance), SV = SV[which.max(abundance)])
# CT_SumName$SV %in% colnames(CT_realp3)
# <-- -->

# if you only have a tree but no taxonomy:
# this here is very similar to OTU clustering he claims:

h1 <- 0.4
ps4 <- tip_glom(ps2, h = h1)

# Figure 4
p2tree <- plot_tree(ps2, method = "treeonly", ladderize = "left", title = "Before Agglomeration")
p3tree <- plot_tree(ps3, method = "treeonly", ladderize = "left", title = "by genus")
p4tree <- plot_tree(ps4, method = "treeonly", ladderize = "left", title = "by height")
grid.arrange(p2tree, p3tree, p4tree, nrow = 1)

p3tree <- plot_tree(ps3, method = "sampledodge", ladderize = "left", title = "by genus")


# -- Abundance value transformation  --
# - going to relative abundances -

# remember the psmelt function to from a phyloseq object to a plottable data.frame
# Note also the log scaling
plot_abundance = function(physeq, title = "",
                          Facet = "Order", Color = "Phylum"){
        # Arbitrary subset, based on Phylum, for plotting
        p1f = subset_taxa(physeq, Phylum %in% c("Firmicutes"))
        mphyseq = psmelt(p1f)
        mphyseq <- subset(mphyseq, Abundance > 0)
        ggplot(data = mphyseq, mapping = aes_string(x = "Group",y = "Abundance")) +
                geom_violin(fill = NA) +
                geom_point(aes_string(col = Color, fill = Color), size = 1, alpha = 0.3,
                           position = position_jitter(width = 0.3)) +
                facet_wrap(facets = Facet) + scale_y_log10() +
                theme(legend.position="none")
}

ps3ra <- transform_sample_counts(ps3, function(x) {x/sum(x)})


plotBefore = plot_abundance(ps3, Color = "Sample")
plotAfter = plot_abundance(ps3ra, Color = "Sample")
# Combine each plot into one graphic.
grid.arrange(nrow = 2, plotBefore, plotAfter)

psOrd = subset_taxa(ps3ra, Order == "Lactobacillales")
plot_abundance(psOrd, Facet = "Genus", Color = NULL)

# --
# ----

# - Preprocessing -
# I add the sample_sums
sample_data(ps)$sample_sums <- sample_sums(ps)

qplot(log10(sample_sums(ps)), binwidth = 0.01)

# <-- NB: here they go back to counts and just log transform with pseudocount 1 as a "variance stabilization transformation"
pslog <- transform_sample_counts(ps, function(x) log(1 + x))
#sample_data(pslog)$age_binned <- cut(sample_data(pslog)$age,
 #                                    breaks = c(0, 100, 200, 400))
out.wuf.log <- ordinate(pslog, method = "MDS", distance = "wunifrac")
evals <- out.wuf.log$values$Eigenvalues
plot_ordination(pslog, out.wuf.log, color = "Group") +
        labs(col = "Group") +
        coord_fixed(sqrt(evals[2] / evals[1]))

rel_abund <- t(apply(otu_table(ps), 1, function(x) x / sum(x)))
qplot(rel_abund[, 12], geom = "histogram") +
        xlab("Relative abundance")
# I do not have any drastic outliers with relative abundance above 90% as they had
# Nowhere they back up their claim that their two outliers have the lowest diversity

# --


# - Different ordination projections -
# NB there was some order mix up in their manuscript

out.bc.log <- ordinate(pslog, method = "MDS", distance = "bray")
evals <-  out.bc.log$values$Eigenvalues
plot_ordination(pslog, out.bc.log, shape = "Group", color = "sample_sums") +
        coord_fixed(sqrt(evals[2] / evals[1])) +
        labs(col = "Binned Age")
# so note: this is including all taxa, just log transformed


# Figure 11
out.dpcoa.log <- ordinate(pslog, method = "DPCoA")
evals <- out.dpcoa.log$eig
plot_ordination(pslog, out.dpcoa.log, color = "Group",
                shape = "Group") +
        coord_fixed(sqrt(evals[2] / evals[1])) # +
        #labs(col = "Binned Age", shape = "Litter")

# Figure 12 interpret axes biplot
plot_ordination(pslog, out.dpcoa.log, type = "species", color = "Phylum") + coord_fixed(sqrt(evals[2] / evals[1]))

# Fig 13
out.wuf.log <- ordinate(pslog, method = "PCoA", distance ="wunifrac")
evals <- out.wuf.log$values$Eigenvalues 
plot_ordination(pslog, out.wuf.log, color = "sample_sums", shape = "Group") +
        coord_fixed(sqrt(evals[2] / evals[1])) 

plot_ordination(pslog, out.wuf.log, color = "Group", shape = "Group") +
        coord_fixed(sqrt(evals[2] / evals[1])) 

# use biplot again to interpret the axes like for DPCoA now for weighted unifrac
plot_ordination(pslog, out.wuf.log, type = "species", color = "Phylum") + coord_fixed(sqrt(evals[2] / evals[1]))
# indeed also in my data DPCoA was better in giving a taxonomic explanation of the axes. 

