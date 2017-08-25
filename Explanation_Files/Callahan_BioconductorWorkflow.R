# packages
.cran_packages <- c("knitr", "phyloseqGraphTest", "phyloseq", "shiny",
                    "miniUI", "caret", "pls", "e1071", "ggplot2", "randomForest",
                    "vegan", "plyr", "dplyr", "ggrepel", "nlme",
                    "reshape2","devtools", "PMA", "structSSI", "ade4",
                    "igraph", "ggnetwork", "intergraph", "scales")
.github_packages <- c("jfukuyama/phyloseqGraphTest")
.bioc_packages <- c("phyloseq", "genefilter", "impute")
# Install CRAN packages (if not already installed)
.inst <- .cran_packages %in% installed.packages()
if (any(!.inst)){
        install.packages(.cran_packages[!.inst],repos = "http://cran.rstudio.com/")
}
.inst <- .github_packages %in% installed.packages()
if (any(!.inst)){
        devtools::install_github(.github_packages[!.inst])
}
.inst <- .bioc_packages %in% installed.packages()
if (any(!.inst)){
        source("http://bioconductor.org/biocLite.R")
        biocLite(.bioc_packages[!.inst])
}


# RUN library for all packages of generalized phyloseq analysis
# in addition
library(ade4)
library(reshape2)
library(ggrepel)
library(caret)
library(igraph)
library(ggnetwork)
library(phyloseqGraphTest)
library(gridExtra)
library(phyloseqGraphTest)
library(nlme)
library(DESeq2)
library(structSSI)
library(plyr)

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
        labs(col = "sample_sumns")
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


# -- PCA on ranks --
# when no transformation helps to bring data to normality, it can be safer to ignore abundances altogether and go for ranks

abund <- otu_table(pslog)
# NB: data still has a 0 where there were zeros, because log(1 + 0) = 0
abund_ranks <- t(apply(abund, 1, rank))
# each SV within each sample, NB: ties were not broken, so you have for example 325.5 in
# The issue he talks about here: If many microbes have very low abundance, there could be large differences in ranks that are actually small in abundance.
# So he introduces with the following transformation a rank threshold, all microbes below this rank are set to 1, the ranks of the others are shifted down

# he uses 329 as a threshold, after seeing the plot below corresponding to figure 15 I realized that I had to use a higher threshold to get to a similar plot
abund_ranks <- abund_ranks - 700
abund_ranks[abund_ranks < 1] <- 1

# - Figure 15 -
#  since the code for figure 15 does not seem to be given (yes see blow), let's try to recapitulate
abund5 <- as.data.frame(abund[1:5,])
abund_ranks5 <- as.data.frame(abund_ranks[1:5,])
colnames(abund5) <- paste("SV", 1:869, sep = "_")
colnames(abund_ranks5) <- paste("SV", 1:869, sep = "_")
abund5 <- as.data.frame(t(abund5))
abund5$SV <- rownames(abund5)
abund_ranks5 <- as.data.frame(t(abund_ranks5))
abund_ranks5$SV <- rownames(abund_ranks5)

abund5 <- gather(abund5, key = Sample, value = logCount, -SV)
abund_ranks5 <- gather(abund_ranks5, key = Sample, value = rank, -SV)

abund5 <- merge(abund5, abund_ranks5, by = c("Sample", "SV"))

Tr <- ggplot(abund5, aes(x = logCount, y = rank, col = Sample))
Tr <- Tr +
        geom_point()

# here his code found at the end of page 21
# abund2 <- as.data.frame(abund)
#melt_abund <- melt(abund, value.name = "abund")
#melt_abund_ranks <- melt(abund_ranks, value.name = "rank")
#joined <- left_join(melt_abund, melt_abund_ranks)
# melt_abund2 <- melt(abund2, value.name = "abund")
abund_df <- melt(abund, value.name = "abund") %>%
        left_join(melt(abund_ranks, value.name = "rank"))
colnames(abund_df) <- c("sample", "seq", "abund", "rank")

sample_ix <- sample(1:nrow(abund_df), 8)
ggplot(abund_df %>%
               filter(sample %in% abund_df$sample[sample_ix])) +
        geom_point(aes(x = abund, y = rank, col = sample),
                   position = position_jitter(width = 0.02), size = .7) +
        labs(x = "Abundance", y = "Thresholded rank") +
        scale_color_brewer(palette = "Set2")
# --


# now the actual PCA from his code
ranks_pca <- dudi.pca(abund_ranks, scannf = F, nf = 3)
row_scores <- data.frame(li = ranks_pca$li,
                         Sample = rownames(abund_ranks))
col_scores <- data.frame(co = ranks_pca$co,
                         seq = colnames(abund_ranks))
tax <- tax_table(ps)@.Data %>%
        data.frame(stringsAsFactors = FALSE)
tax$seq <- rownames(tax)
main_orders <- c("Clostridiales", "Bacteroidales", "Lactobacillales",
                 "Coriobacteriales")
tax$Order[!(tax$Order %in% main_orders)] <- "Other"
tax$Order <- factor(tax$Order, levels = c(main_orders, "Other"))
tax$otu_id <- seq_len(ncol(otu_table(ps)))

sampleDF <- sample_data(pslog)
sampleDF$Sample <- rownames(sampleDF)
row_scores <- row_scores %>%
        left_join(sampleDF)

col_scores <- col_scores %>%
        left_join(tax)


# - Fig 16 - 
evals_prop <- 100 * (ranks_pca$eig / sum(ranks_pca$eig))
Tr <- ggplot() +
        # geom_point(data = row_scores, aes(x = li.Axis1, y = li.Axis2, shape = Group)) +
        geom_point(data = row_scores, aes(x = li.Axis1, y = li.Axis2), shape = 2) +
        geom_point(data = col_scores, aes(x = 25 * co.Comp1, y = 25 * co.Comp2, col = Order),
                   size = .3, alpha = 0.6) +
        scale_color_brewer(palette = "Set2") +
        facet_grid(~ Group) +
        guides(col = guide_legend(override.aes = list(size = 3))) +
        labs(x = sprintf("Axis1 [%s%% variance]", round(evals_prop[1], 2)),
             y = sprintf("Axis2 [%s%% variance]", round(evals_prop[2], 2))) +
        coord_fixed(sqrt(ranks_pca$eig[2] / ranks_pca$eig[1])) +
        theme(panel.border = element_rect(color = "#787878", fill = alpha("white", 0)))

# --
# I see the point in showing a biplot again, i.e. not only showing the samples but also which order contributed to the axes, but this facet grid is weird

Tr <- ggplot() +
        geom_point(data = col_scores, aes(x = 25 * co.Comp1, y = 25 * co.Comp2, col = Order),
                   size = 1, alpha = 0.6) +
        # geom_point(data = row_scores, aes(x = li.Axis1, y = li.Axis2, shape = Group)) +
        geom_point(data = row_scores, aes(x = li.Axis1, y = li.Axis2, shape = Group), size = 1.5) +
        scale_color_brewer(palette = "Set2") +
        # facet_grid(~ Group) +
        guides(col = guide_legend(override.aes = list(size = 3))) +
        labs(x = sprintf("Axis1 [%s%% variance]", round(evals_prop[1], 2)),
             y = sprintf("Axis2 [%s%% variance]", round(evals_prop[2], 2))) +
        coord_fixed(sqrt(ranks_pca$eig[2] / ranks_pca$eig[1])) +
        theme(panel.border = element_rect(color = "#787878", fill = alpha("white", 0))) +
        theme_bw()

# ----


# -- Canonical correspondance --
# biplots where both species signatures and environmental characteristics go in

ps_ccpna <- ordinate(pslog, "CCA", formula = pslog ~ Group + Sample.Integrity)

ps_scores <- vegan::scores(ps_ccpna)
sites <- data.frame(ps_scores$sites)
sites$Sample <- rownames(sites)
sampleDF <- sample_data(pslog)
sampleDF$Sample <- rownames(sampleDF)
sites <- sites %>%
        left_join(sampleDF)
species <- data.frame(ps_scores$species)
species$otu_id <- seq_along(colnames(otu_table(ps)))
species <- species %>%
        left_join(tax) # here joined using this otu_id

# - Fig 17 -

evals_prop <- 100 * ps_ccpna$CCA$eig[1:2] / sum(ps_ccpna$CA$eig)
Tr <- ggplot() +
        geom_point(data = sites, aes(x = CCA1, y = CCA2), shape = 2, alpha = 0.5) +
        geom_point(data = species, aes(x = CCA1, y = CCA2, col = Order), size = 0.5) +
        geom_text_repel(data = species %>% filter(CCA2 < -0.80106),
                        aes(x = CCA1, y = CCA2, label = otu_id),
                        size = 1.5, segment.size = 0.1) +
        facet_grid(. ~ Group) +
        guides(col = guide_legend(override.aes = list(size = 3))) +
        labs(x = sprintf("Axis1 [%s%% variance]", round(evals_prop[1], 2)),
             y = sprintf("Axis2 [%s%% variance]", round(evals_prop[2], 2))) +
        scale_color_brewer(palette = "Set2") +
        coord_fixed(sqrt(ps_ccpna$CCA$eig[2] / ps_ccpna$CCA$eig[1])*0.33) +
        theme(panel.border = element_rect(color = "#787878", fill = alpha("white", 0)))


# well again, I do not like the faceting, and clearly not the labelling of some species that are outliers

Tr <- ggplot() +
        geom_point(data = species, aes(x = CCA1, y = CCA2, col = Order), size = 1) +
        geom_point(data = sites, aes(x = CCA1, y = CCA2, shape = Group), size = 2, alpha = 0.8) +
        # geom_text_repel(data = species %>% filter(CCA2 < -0.80106),
        #                 aes(x = CCA1, y = CCA2, label = otu_id),
        #                 size = 1.5, segment.size = 0.1) +
        # facet_grid(. ~ Group) +
        guides(col = guide_legend(override.aes = list(size = 3))) +
        labs(x = sprintf("Axis1 [%s%% variance]", round(evals_prop[1], 2)),
             y = sprintf("Axis2 [%s%% variance]", round(evals_prop[2], 2))) +
        scale_color_brewer(palette = "Set2") +
        coord_fixed(sqrt(ps_ccpna$CCA$eig[2] / ps_ccpna$CCA$eig[1])*0.33) +
        theme(panel.border = element_rect(color = "#787878", fill = alpha("white", 0))) +
        theme_bw()

# run it two times and note that sample integrity hardly changed anything.

# --

# - Fig 18 -
# basically the same just using Sample.Integrity instead of Group

Tr <- ggplot() +
        geom_point(data = species, aes(x = CCA1, y = CCA2, col = Order), size = 0.5) +
        geom_point(data = sites, aes(x = CCA1, y = CCA2), shape = 2, alpha = 0.5) +
        # geom_text_repel(data = species %>% filter(CCA2 < -0.80106),
        #                 aes(x = CCA1, y = CCA2, label = otu_id),
        #                 size = 1.5, segment.size = 0.1) +
        facet_grid(. ~ Sample.Integrity) +
        guides(col = guide_legend(override.aes = list(size = 3))) + # makes the legend dots big even though they are small, smart!
        labs(x = sprintf("Axis1 [%s%% variance]", round(evals_prop[1], 2)),
             y = sprintf("Axis2 [%s%% variance]", round(evals_prop[2], 2))) +
        scale_color_brewer(palette = "Set2") +
        coord_fixed(sqrt(ps_ccpna$CCA$eig[2] / ps_ccpna$CCA$eig[1])*0.45 ) +
        theme(panel.border = element_rect(color = "#787878", fill = alpha("white", 0)))

# again, I like much more the non faceted version

Tr <- ggplot() +
        geom_point(data = species, aes(x = CCA1, y = CCA2, col = Order), size = 1) +
        geom_point(data = sites, aes(x = CCA1, y = CCA2, shape = Sample.Integrity), size = 2, alpha = 0.8) +
        # geom_text_repel(data = species %>% filter(CCA2 < -0.80106),
        #                 aes(x = CCA1, y = CCA2, label = otu_id),
        #                 size = 1.5, segment.size = 0.1) +
        # facet_grid(. ~ Sample.Integrity) +
        guides(col = guide_legend(override.aes = list(size = 3))) + # makes the legend dots big even though they are small, smart!
        labs(x = sprintf("Axis1 [%s%% variance]", round(evals_prop[1], 2)),
             y = sprintf("Axis2 [%s%% variance]", round(evals_prop[2], 2))) +
        scale_color_brewer(palette = "Set2") +
        coord_fixed(sqrt(ps_ccpna$CCA$eig[2] / ps_ccpna$CCA$eig[1])*0.45 ) +
        theme(panel.border = element_rect(color = "#787878", fill = alpha("white", 0))) +
        theme_bw()

# --
# ----


# -- Supervised Learning --

# supervised meaning training data set with known outcome
# caret package wraps many prediction algorithms available in R
# Here: can we predict age from microbiome composition?
# ask yourself if this is at all interesting

# - Partial Least Squares (PLS) -

# split into training and test data sets

# setup_example(c("phyloseq", "ggplot2", "caret", "plyr", "dplyr"))
#sample_data(pslog)$age2 <- cut(sample_data(pslog)$age, c(0, 100, 400))
sample_data(pslog)$age2 <- sample_data(pslog)$Group
dataMatrix <- data.frame(age = sample_data(pslog)$age2, otu_table(pslog))
# take 8 mice at random to be the training set, and the remaining 4 the test set
# since I have not done mice longitudinal we can only take a random selection of samples
#trainingMice <- sample(unique(sample_data(pslog)$host_subject_id), size = 8)
#inTrain <- which(sample_data(pslog)$host_subject_id %in% trainingMice)
inTrain <- sample(1:nrow(dataMatrix), 40)
training <- dataMatrix[inTrain,]
testing <- dataMatrix[-inTrain,]
# train the model
plsFit <- train(age ~ ., data = training,
                method = "pls", preProc = "center")
# test the prediction on the test data
plsClasses <- predict(plsFit, newdata = testing)
table(plsClasses, testing$age)

# that is doing a pretty good job in predicting definitely for the young samples
# I change order a bit and go directly to the PLS biplot that should help to interpret the results

# = Fig 19 = 

pls_biplot <- list("loadings" = loadings(plsFit$finalModel),
                   "scores" = scores(plsFit$finalModel))
class(pls_biplot$scores) <- "matrix"
pls_biplot$scores <- data.frame(sample_data(pslog)[inTrain, ],
                                pls_biplot$scores)
tax <- tax_table(ps)@.Data %>%
        data.frame(stringsAsFactors = FALSE)
main_orders <- c("Clostridiales", "Bacteroidales", "Lactobacillales",
                 "Coriobacteriales", "Verrucomicrobiales")
tax$Order[!(tax$Order %in% main_orders)] <- "Other"
tax$Order <- factor(tax$Order, levels = c(main_orders, "Other"))
class(pls_biplot$loadings) <- "matrix"
pls_biplot$loadings <- data.frame(tax, pls_biplot$loadings)

Tr <- ggplot() +
        geom_point(data = pls_biplot$scores,
                   aes(x = Comp.1, y = Comp.2), shape = 2) +
        geom_point(data = pls_biplot$loadings,
                   aes(x = 25 * Comp.1, y = 25 * Comp.2, col = Order),
                   size = 0.3, alpha = 0.6) +
        scale_color_brewer(palette = "Set2") +
        labs(x = "Axis1", y = "Axis2") +
        guides(col = guide_legend(override.aes = list(size = 3))) +
        facet_grid( ~ age2) +
        theme(panel.border = element_rect(color = "#787878", fill = alpha("white", 0)))

# again, I like much more the non faceted version

Tr <- ggplot() +
        geom_point(data = pls_biplot$loadings,
                   aes(x = 25 * Comp.1, y = 25 * Comp.2, col = Order),
                   size = 0.3, alpha = 0.6) +
        geom_point(data = pls_biplot$scores,
                   aes(x = Comp.1, y = Comp.2, shape = Group), size = 2, alpha = 0.8) +
        labs(x = "Axis1", y = "Axis2") +
        guides(col = guide_legend(override.aes = list(size = 3))) +
        theme(panel.border = element_rect(color = "#787878", fill = alpha("white", 0))) +
        theme_bw()        

# based on the biplot you almost ask why not all predictions were correct       

# ==

# --

# - random forest (rf) -
# NB: it is slower than the partial least squares method
# but it performes better in distinguishing MiddleAged from old 

rfFit <- train(age ~ ., data = training, method = "rf",
               preProc = "center", proximity = TRUE)
rfClasses <- predict(rfFit, newdata = testing)
table(rfClasses, testing$age)

# == Fig 20 ==
# a random forest proximity plot, a distance is calculated between samples based on how frequently samples occur in the same tree partition
#  in the rf bootstrapping procedure.
# The PCoA gives a glimpse into the otherwise complex rf classification mechanism

rf_prox <- cmdscale(1 - rfFit$finalModel$proximity) %>%
        data.frame(sample_data(pslog)[inTrain, ])
Tr <- ggplot(rf_prox) +
        geom_point(aes(x = X1, y = X2, col = Group),
                   size = 2, alpha = 0.6) +
        scale_color_manual(values = c("#A66EB8", "#238DB5", "#748B4F")) +
        guides(col = guide_legend(override.aes = list(size = 3))) +
        labs(col = "", x = "Axis1", y = "Axis2")

# could be maybe used to find the samples that are hard to classify
# ====
# == Fig 21 ==

# FIND the microbe with the most influence in the random forest prediction

impSVname <- as.vector(tax_table(ps)[which.max(importance(rfFit$finalModel)), c("Family", "Genus")])
impSVname <- paste(impSVname, collapse = " ")
# NB: some are just as important:
plot(importance(rfFit$finalModel))

impSV <- as.vector(otu_table(pslog)[,which.max(importance(rfFit$finalModel))])
maxImpDF <- data.frame(sample_data(pslog), abund = impSV)
Tr <- ggplot(maxImpDF) + geom_histogram(aes(x = abund)) +
        facet_grid(Group ~ .) +
        labs(x = paste("Abundance of ", impSVname), y = "Number of samples")

# this is nice how different it is, it is only present in the young mice
# ====
# --
# ----

# -- Graph-based visualization and testing --
# - Creating and plotting graphs -
# Phyloseq has functionality for creating graphs based on thresholding a distance matrix, and the resulting networks can
# be plotting using the ggnetwork
# This package overloads the ggplot syntax, so you can use the function ggplot on an
# igraph object and add geom_edges and geom_nodes geoms to plot the network.

# here it would be better to have more variables to illustrate potential grouping, such as litter as in their example

# == Fig 22 ==

net <- make_network(ps, max.dist=0.6) # used by default jaccard distance
sampledata <- data.frame(sample_data(ps))
# V(net)$id <- sampledata[names(V(net)), "host_subject_id"]
V(net)$id <- as.character(sampledata[names(V(net)), "Group"])
# V(net)$litter <- sampledata[names(V(net)), "family_relationship"]
V(net)$integrity <- sampledata[names(V(net)), "Sample.Integrity"]
ggplot(net, aes(x = x, y = y, xend = xend, yend = yend), layout = "fruchtermanreingold") +
        geom_edges(color = "darkgray") +
        geom_nodes(aes(color = id, shape = integrity)) +
        theme(axis.text = element_blank(), axis.title = element_blank(),
              legend.key.height = unit(0.5,"line")) +
        guides(col = guide_legend(override.aes = list(size = .25)))

# ==
# --

# - Graph-based two sample tests -
# introduced by Friedman and Rafsky as a generalization of the Wald-Wlfowity runs test. 
# They proposed the use of a minimum spanning tree (MST) based on the distances between the samples, and
# then counting the number of edges on the tree that were between samples in different groups
# No matter what graph we build between the samples, we can approximate a null distribution
# by permuting the labels of the nodes of the graph.

# ++ Minimum Spanning Tree (MST) ++
# We first perform a test using an MST with Jaccard dissimilarity
# the compare family_relationship (read), I compare influence of Sample integrity, do samples of different integrity come from the same distribution
# so just permutate the integrity labels

# gt <- graph_perm_test(ps, "family_relationship", grouping = "host_subject_id",
#                      distance = "jaccard", type = "mst")
# NB: I could do no grouping because: all values of sampletype must be the same within each level of grouping
# like in their case within a host_subject_id all samples of course were from the same litter
gt <- graph_perm_test(ps, sampletype = "Sample.Integrity", distance = "jaccard", type = "mst")
gt$pval

# == Fig 23 ==
plotNet1 = plot_test_network(gt) + theme(legend.text = element_text(size = 8),
                                       legend.title = element_text(size = 9))
plotPerm1=plot_permutations(gt)
grid.arrange(ncol = 2, plotNet1, plotPerm1)

# do the same for the groups
gt2 <- graph_perm_test(ps, sampletype = "Group", distance = "jaccard", type = "mst")
gt2$pval

# == Fig 23 ==
plotNet1 = plot_test_network(gt2) + theme(legend.text = element_text(size = 8),
                                         legend.title = element_text(size = 9))
plotPerm1=plot_permutations(gt2)
grid.arrange(ncol = 2, plotNet1, plotPerm1)
# ====
# ++++

# ++ Nearest neigbors ++
# The k-nearest neighbors graph is obtained by putting an edge between two samples whenever one of them is in the set
# of k-nearest neighbors of the other. We see from Figure 24 that if a pair of samples has an edge between them in the
# nearest neighbor graph, they are overwhelmingly likely to be in the same litter.

gt <- graph_perm_test(ps, sampletype = "Sample.Integrity", distance = "jaccard", type = "knn", knn = 1)
gt$pval

# == Fig 24 ==
plotNet2 = plot_test_network(gt) + theme(legend.text = element_text(size = 8),
                                         legend.title = element_text(size = 9))
plotPerm2=plot_permutations(gt)
grid.arrange(ncol = 2, plotNet2, plotPerm2)

# again also for the groups
gt <- graph_perm_test(ps, sampletype = "Group", distance = "jaccard", type = "knn", knn = 1)
gt$pval

plotNet2 = plot_test_network(gt) + theme(legend.text = element_text(size = 8),
                                         legend.title = element_text(size = 9))
plotPerm2=plot_permutations(gt)
grid.arrange(ncol = 2, plotNet2, plotPerm2)

# Not sure if all this should actually always just be for two groups, but well, no complains about three here

# ====

# same test for bray curtis and two-nearest neighbors
gt <- graph_perm_test(ps, sampletype = "Sample.Integrity", distance = "bray", type = "knn", knn = 2)
gt$pval

plotNet2 = plot_test_network(gt) + theme(legend.text = element_text(size = 8),
                                         legend.title = element_text(size = 9))
plotPerm2=plot_permutations(gt)
grid.arrange(ncol = 2, plotNet2, plotPerm2)

# with Group

gt <- graph_perm_test(ps, sampletype = "Group", distance = "bray", type = "knn", knn = 2)
gt$pval

plotNet2 = plot_test_network(gt) + theme(legend.text = element_text(size = 8),
                                         legend.title = element_text(size = 9))
plotPerm2=plot_permutations(gt)
grid.arrange(ncol = 2, plotNet2, plotPerm2)


# SO I WONDER: IS THIS BETTER THAN JUST ADONIS ON PCOA, in the end it semms all comes down to the dissimilarity measures

# ++++

# ++ Distance threshold ++

# Another way of making a graph between samples is to threshold the distance matrix, this is called a geometric graph27.
# The testing function lets the user supply an absolute distance threshold; alternatively, it can find a distance threshold
# such that there are a prespecified number of edges in the graph.
# Below we make a graph to have twice as many edges as samples
# Heuristically, the graph we obtain isn’t as good,
# because there are many singletons. This reduces power, and so if the thresholded graph has this many singletons it is
# better to either modify the threshold or consider a MST or k-nearest neighbors graph.

# == Fig 25 ==
gt <- graph_perm_test(ps, "Group",
                      distance = "bray", type = "threshold.nedges", nedges = 120, # I have 61 samples
                      keep.isolates = FALSE)
gt$pval
plotNet3 <- plot_test_network(gt) + theme(legend.text = element_text(size = 8),
                                        legend.title = element_text(size = 9))
plotPerm3 <- plot_permutations(gt)
grid.arrange(ncol = 2, plotNet3, plotPerm3)

# what happens with Sample.Integrity

gt <- graph_perm_test(ps, "Sample.Integrity",
                      distance = "bray", type = "threshold.nedges", nedges = 120, # I have 61 samples
                      keep.isolates = FALSE)
gt$pval
plotNet3 <- plot_test_network(gt) + theme(legend.text = element_text(size = 8),
                                          legend.title = element_text(size = 9))
plotPerm3 <- plot_permutations(gt)
grid.arrange(ncol = 2, plotNet3, plotPerm3)

# ====
# ++++

# ++ Linear modeling ++
# It is often of interest to evaluate the degree to which microbial community diversity reflects characteristics of the environment
# from which it was sampled. Unlike ordination, the purpose of this analysis is not to develop a representation of
# many bacteria with respect to sample characteristics; rather, it is to describe how a single measure of overall community
#structure (In particular, it need not be limited to diversity – defining univariate measures of community stability is also
 #          common, for example.) is associated with sample characteristics
# This is a somewhat simpler statistical goal, and can be addressed through linear modeling, for which there are a range of approaches in R. 

# As an example, we will used a mixed-effects model to study the relationship between mouse microbial community diversity and 
# the age and litter (me age/group and Sample.Integrity) variables that have been our focus so far.
# This choice was motivated by the observation that younger mice have noticeably lower Shannon diversities, 
# but that different mice have different baseline diversities. The mixed-effects model is 
# a starting point for formalizing this observation.

# compute shannon diversity associated with each sample and join it with sample annotation

ps_alpha_div <- estimate_richness(ps, split = TRUE, measure = "Shannon")
ps_alpha_div$SampleID <- rownames(ps_alpha_div) %>%
        as.factor()
ps_samp <- sample_data(ps) %>%
        # unclass() %>%
        data.frame() 
ps_samp$SampleID <- rownames(ps_samp) %>% as.factor()

ps_samp <- left_join(ps_samp, ps_alpha_div, by = "SampleID")
ps_samp <- reshape2::melt(ps_samp, measure.vars = "Shannon",
             variable.name = "diversity_measure",
             value.name = "alpha_diversity")

# to being able to recapitulate this mixed-effects model, I have to make up some host_subject_id data.
# their data is longitudinal, i.e. the same mice have been analysed at different ages.
# we had 24 Young mice, 25 Old, and 12 middle age, I pretend it had been 4 mice
set.seed(1332)
ps_samp$host_subject_id <- ""
ps_samp$host_subject_id[ps_samp$Group == "Young"] <- paste("M_", sample(4, size = 24, replace = T), sep = "")
ps_samp$host_subject_id[ps_samp$Group == "Old"] <- paste("M_", sample(4, size = 25, replace = T), sep = "")
ps_samp$host_subject_id[ps_samp$Group == "MiddleAged"] <- paste("M_", sample(4, size = 12, replace = T), sep = "")


# reorder facets from lowest to highest diversity
diversity_means <- ps_samp %>%
        group_by(host_subject_id) %>%
        summarise(mean_div = mean(alpha_diversity)) %>%
        arrange(mean_div)
ps_samp$host_subject_id <- factor(ps_samp$host_subject_id,
                                  diversity_means$host_subject_id)

# now the mixed effect model
alpha_div_model <- lme(fixed = alpha_diversity ~ Group, data = ps_samp,
                       random = ~ 1 | host_subject_id)

# == Fig 26 ==

new_data <- expand.grid(host_subject_id = levels(ps_samp$host_subject_id),
                        Group = levels(ps_samp$Group))
# I had to change to character so that the predict function works
new_data$host_subject_id <- as.character(new_data$host_subject_id)
new_data$Group <- as.character(new_data$Group)
new_data$pred <- predict(alpha_div_model, newdata = new_data)
# I had to play a bit with this prediction, so I'm not fully sure of it
# Therefore compare as an alternative 
new_data_pred <- cbind(ps_samp[, c("SampleID", "host_subject_id", "Group")], predict(alpha_div_model))
new_data_pred <- new_data_pred[!duplicated(new_data_pred[,c("host_subject_id", "Group")]), ]
new_data_pred <- new_data_pred[, -1]
colnames(new_data_pred)[3] <- "pred"

new_data$host_subject_id <- factor(new_data$host_subject_id, diversity_means$host_subject_id)
new_data_pred$host_subject_id <- factor(new_data_pred$host_subject_id, diversity_means$host_subject_id)
new_data$Group <- factor(new_data$Group, levels = levels(ps_samp$Group), ordered = TRUE)
new_data_pred$Group <- factor(new_data_pred$Group, levels = levels(ps_samp$Group), ordered = TRUE)

# I think the key thing for this figure are the prediction intervals
# would have also been nice to look and describe
summary(alpha_div_model)
X <- model.matrix(eval(eval(alpha_div_model$call$fixed)[-2]),
                  new_data[-ncol(new_data)])
pred_var_fixed <- diag(X %*% alpha_div_model$varFix %*% t(X))
new_data$pred_var <- pred_var_fixed + alpha_div_model$sigma ^ 2

X2 <- model.matrix(eval(eval(alpha_div_model$call$fixed)[-2]),
                   new_data[-ncol(new_data_pred)])

pred_var_fixed2 <- diag(X2 %*% alpha_div_model$varFix %*% t(X2))
new_data_pred$pred_var <- pred_var_fixed2 + alpha_div_model$sigma ^ 2


Tr <- ggplot(ps_samp %>% left_join(new_data)) +
        geom_errorbar(aes(x = Group, ymin = pred - 2 * sqrt(pred_var), # so assuming pred_var is variance, the error bars are ~95% intervals, because 2 SD away
                          ymax = pred + 2 * sqrt(pred_var)),
                      col = "#858585", size = .1) +
        geom_point(aes(x = Group, y = alpha_diversity,
                       col = Sample.Integrity), size = 1.5) +
        facet_wrap(~host_subject_id) +
        scale_y_continuous(limits = c(3, 5), breaks = seq(0, 5, .5)) +
        scale_color_brewer(palette = "Set2") +
        labs(x = "Group", y = "Shannon Diversity", color = "Sample.Integrity") +
        guides(col = guide_legend(override.aes = list(size = 4))) +
        theme(panel.border = element_rect(color = "#787878", fill = alpha("white", 0)),
              axis.text.x = element_text(angle = -90, size = 6),
              axis.text.y = element_text(size = 6))

Tr1 <- ggplot(ps_samp %>% left_join(new_data_pred)) +
        geom_errorbar(aes(x = Group, ymin = pred - 2 * sqrt(pred_var), # so assuming pred_var is variance, the error bars are ~95% intervals, because 2 SD away
                          ymax = pred + 2 * sqrt(pred_var)),
                      col = "#858585", size = .1) +
        geom_point(aes(x = Group, y = alpha_diversity,
                       col = Sample.Integrity), size = 1.5) +
        facet_wrap(~host_subject_id) +
        scale_y_continuous(limits = c(3, 5), breaks = seq(0, 5, .5)) +
        scale_color_brewer(palette = "Set2") +
        labs(x = "Group", y = "Shannon Diversity", color = "Sample.Integrity") +
        guides(col = guide_legend(override.aes = list(size = 4))) +
        theme(panel.border = element_rect(color = "#787878", fill = alpha("white", 0)),
              axis.text.x = element_text(angle = -90, size = 6),
              axis.text.y = element_text(size = 6))

# so if you want to use it, you have to figure out the correct prediction method because you see that the error bars are different

# ====

# ++++

# ++ Hierarchical multiple testing ++ and DESeq2

# Hierarchical from higher taxonomies down to reduce multiple testing issues:
# it is likely that if one Ruminococcus species is strongly associated
# with age, then others are as well. To integrate this information29,30, proposed a hierarchical testing procedure, where
# taxonomic groups are only tested if higher levels are found to be be associated. In the case where many related species
# have a slight signal, this pooling of information can increase power.

# First the DESEQ count normalisation, it seems they stress again variance stabilizing transformation, which can be used instead of log
# transformation. 
# If you compare with log transform (plus pseudocount): 
# One difference is that, after accounting for size factors, the
# histogram of row sums for DESeq is more spread out in the lower values, refer to Figure 27. This is the motivation of
# using such a transformation, although for high abundance counts, it is equivalent to the log, for lower and mid range
# abundances it does not crush the data and yields more powerful results

# ps_dds <- phyloseq_to_deseq2(ps, ~ Group + Sample.Integrity)
sample_data(ps)$Group <- factor(sample_data(ps)$Group, ordered = FALSE)
# NB: the factor is not allowed to be ordered because otherwise DESeq2 thinks they are all the same value!
ps_dds <- phyloseq_to_deseq2(ps, ~ Group)

Test <- varianceStabilizingTransformation(ps_dds, blind = TRUE, fitType = "parametric") # not sure why this is run here

ps_dds <- estimateSizeFactors(ps_dds)
ps_dds <- estimateDispersions(ps_dds)
abund <- getVarianceStabilizedData(ps_dds)

# We use the structSSI to perform the hierarchical testing33. For more convenient printing, we first shorten the names of
# each microbe.

# short_names <- substr(rownames(abund), 1, 5) %>%
  #      make.names(unique = TRUE)
short_names <- paste("SV_", 1:nrow(abund), sep = "")

rownames(abund) <- short_names

# Unlike standard multiple hypothesis testing, the hierarchical testing procedure needs univariate tests for each higherlevel
# taxonomic group, not just every bacteria. A helper function, treePValues, is available for this; it expects an
# edgelist encoding parent-child relationships, with the first row specifying the root node.
# NOT sure if I had already only univariate


el <- phy_tree(pslog)$edge # why here from pslog??
el0 <- el
el0 <- el0[nrow(el):1, ]
el_names <- c(short_names, seq_len(phy_tree(pslog)$Nnode))
el[, 1] <- el_names[el0[, 1]]
el[, 2] <- el_names[as.numeric(el0[, 2])]
unadj_p <- treePValues(el, abund, sample_data(pslog)$Group)
# so this was for the DESeq data, let's try with the pslog
abund_log <- t(as(otu_table(pslog), "matrix"))
rownames(abund_log) <- short_names
unadj_p_log <- treePValues(el, abund_log, sample_data(pslog)$Group)


# We can now correct p-value using the hierarchical testing procedure. The test results are guaranteed to control several
# variants of FDR control, but at different levels; we defer details to 29,30,33.

hfdr_res <- hFDR.adjust(unadj_p, el, .75)
summary(hfdr_res)

hfdr_res_log <- hFDR.adjust(unadj_p_log, el, .75)
summary(hfdr_res_log)

# so not exactly the same but clearly related
# ok, clearly have to think about it, but looks nice as if edges are in, the ones without SV, but is it significant between old and young or old and middle age and sso on


# == Fig 27 ==

# because I get negative sums for the variance stabilized data, I wondered why theirs was in such a similar range to the pslog
# in general, I'm not sure if this variance stabilization should be done at all for testing differential abundance
abund_sums <- rbind(data.frame(sum = colSums(abund),
                               sample = colnames(abund),
                               type = "DESeq2"),
                    data.frame(sum = sample_sums(pslog),
                               sample = sample_names(pslog),
                               type = "log(1 + x)"))
Tr <- ggplot(abund_sums) +
        geom_histogram(aes(x = sum), binwidth = 20) +
        facet_grid(type ~ .) +
        xlab("Total abundance within sample")

# ====

# == Fig 28 ==

plot(hfdr_res, height = 12000) # opens in new browser window
plot(hfdr_res_log, height = 12000)

# you can see that the tree is the same, but different edges have not been tested

# The plot opens in a new browser – a static screenshot of a subtree is displayed in Figure 28. Nodes are shaded
# according to p-values, from blue to orange, representing the strongest to weakest associations. Grey nodes were
# never tested, to focus power on more promising subtrees. Scanning the full tree, it becomes clear that the association
# between age group and bacterial abundance is present in only a few isolated taxonomic groups, but that it is quite strong
# in those groups. To give context to these results, we can retrieve the taxonomic identity of the rejected hypotheses.

OTU <- t(counts(ps_dds, normalized = TRUE))





