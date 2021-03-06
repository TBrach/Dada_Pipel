# --
#######################################
### FUNCTION: check_phyla_distribution#
#######################################
# NB: throws error if "Phylum" is not in colnames(tax_table(physeq))
# outputs a data.frame, summarising how many taxa each phylum has, and


check_phyla_distribution <- function(physeq) {
        
        
        if (phyloseq::taxa_are_rows(physeq)) {
                physeq <- t(physeq)
        }
        
        df_ab_prev <- data.frame(Taxon_No = 1:ntaxa(physeq), 
                                 total_counts = taxa_sums(physeq),
                                 prevalence = colSums(as(otu_table(physeq), "matrix") != 0))
        
        
        df_ab_prev <- cbind(df_ab_prev, tax_table(physeq))
        
        PhylaDistribution <- dplyr::summarise(group_by(df_ab_prev, Phylum), 
                                              taxa = n(), 
                                              PC_of_taxa = round(100*taxa/ntaxa(ps),1),
                                              PC_of_counts = round(100*sum(total_counts)/sum(otu_table(physeq)), 1),
                                              PC_of_prevalence = round(100*sum(prevalence)/sum(otu_table(physeq) != 0), 1),
                                              mean_taxa_sum = round(mean(total_counts)),
                                              median_taxa_sum = round(median(total_counts)),
                                              mean_prevalence_in_PC = round(100*mean(prevalence)/nsamples(ps), 1)) %>% 
                arrange(desc(PC_of_counts), desc(taxa), desc(PC_of_prevalence))
        
        PhylaDistribution
        
}
# --

# --
#######################################
### FUNCTION: get_assignemnt_distribution
#######################################


get_assignemnt_distribution <- function(physeq){
        
        taxa <- tax_table(physeq)
        
        countNA <- apply(taxa, 2, function(x){sum(is.na(x))})
        countNonNA <- apply(taxa, 2, function(x){sum(!is.na(x))})
        # ambiguous <- apply(taxa, 2, function(x){sum(grepl("/", x))})
        total <- countNA + countNonNA
        assignment_distribution <- data.frame(assigned = countNonNA, 
                                              # assigned_unamb = countNonNA - ambiguous,
                                              total = total,
                                              PC_assigned = round(100*countNonNA/total, 1))
                                              # PC_assigned_unamb = round(100*(countNonNA - ambiguous)/total, 1))
        
}
# --



# --
#######################################
### FUNCTION: check_assignment_vs_abundance
#######################################
# checks if there is a trend for better assignment for more abundant SVs

check_assignment_vs_abundance <- function(physeq, abundanceQuantiles = seq(0, 90, by = 10)){
        
        if (phyloseq::taxa_are_rows(physeq)) {
                physeq <- t(physeq)
        }
        
        taxa = tax_table(physeq) 
        seqtab = as(otu_table(physeq), "matrix")
        
        
        total_counts <- colSums(seqtab)
        abQuantiles <- quantile(total_counts, probs = abundanceQuantiles/100)
        remaining_SVs <- lapply(abQuantiles, function(quant) {
                total_counts >= quant
        })
        No_ASVs <- sapply(remaining_SVs, sum)
        filtered.taxas <- lapply(remaining_SVs, function(indexes){
                taxa[indexes,]
        })
        
        assignment_distributions <- lapply(filtered.taxas, get_assignemnt_distribution)
        
        PC_assigned <- sapply(assignment_distributions, function(distri){distri[["PC_assigned"]]})
        
        rownames(PC_assigned) <- colnames(taxa)
        colnames(PC_assigned) <- paste("Ab_", abundanceQuantiles, "_", round(abQuantiles), "_", No_ASVs, sep = "")
        
        PC_assigned <- as.data.frame(PC_assigned)
        PC_assigned$Level <- rownames(PC_assigned)
        df_long <- gather(PC_assigned, key = Ab, value = PC, -Level)
        df_long <- separate(df_long, col = Ab, into = c("Type", "Quant", "abQuant", "No_ASVs"), sep = "_")
        df_long$Level <- factor(df_long$Level, levels = colnames(taxa), ordered = TRUE)
        df_long$Quant <- factor(df_long$Quant, levels = as.character(abundanceQuantiles), ordered = TRUE)
        
        tr <- ggplot(df_long, aes(x = Quant, y = PC, col = Level))
        tr <- tr + 
                geom_point(size = 2) +
                scale_colour_manual("", values = cbPalette[2:8]) +
                scale_x_discrete(breaks = abundanceQuantiles, labels = paste(abundanceQuantiles, " (", No_ASVs, ")", sep = "")) +
                ylab("percentage of taxa assigned") +
                xlab("total counts quantile (No of remaining taxa)") +
                theme_bw() +
                theme(axis.text.x = element_text(angle=90, vjust=0.5))
        
        list(PC_assigned, tr)
        
}
# --




# --
#######################################
### FUNCTION: check_assignment_vs_prevalence
#######################################
# checks if there is a trend for better assignment for more prealent SVs

check_assignment_vs_prevalence <- function(physeq, prevalences = seq(0, 90, by = 10)){
        
        if (phyloseq::taxa_are_rows(physeq)) {
                physeq <- t(physeq)
        }
        
        taxa = tax_table(physeq) 
        seqtab = as(otu_table(physeq), "matrix")
        
        prev <- colSums(seqtab != 0)
        remaining_SVs <- lapply(prevalences, function(preva) {
                prev > (preva/100)*nrow(seqtab)
        })
        No_ASVs <- sapply(remaining_SVs, sum)
        filtered.taxas <- lapply(remaining_SVs, function(indexes){
                taxa[indexes,]
        })
        
        assignment_distributions <- lapply(filtered.taxas, get_assignemnt_distribution)
        
        PC_assigned <- sapply(assignment_distributions, function(distri){distri[["PC_assigned"]]})
        
        rownames(PC_assigned) <- colnames(taxa)
        colnames(PC_assigned) <- paste("Prev_", prevalences, "_", No_ASVs, sep = "")
        
        PC_assigned <- as.data.frame(PC_assigned)
        PC_assigned$Level <- rownames(PC_assigned)
        df_long <- gather(PC_assigned, key = Prev, value = PC, -Level)
        df_long <- separate(df_long, col = Prev, into = c("Type", "Prevalence", "No_ASVs"), sep = "_")
        df_long$Level <- factor(df_long$Level, levels = colnames(taxa), ordered = TRUE)
        df_long$Prevalence <- factor(df_long$Prevalence, levels = as.character(prevalences), ordered = TRUE)
        
        
        tr <- ggplot(df_long, aes(x = Prevalence, y = PC, col = Level))
        tr <- tr + 
                geom_point(size = 2) +
                scale_colour_manual("", values = cbPalette[2:8]) +
                scale_x_discrete(breaks = prevalences, labels = paste(prevalences, " (", No_ASVs, ")", sep = "")) +
                ylab("percentage of taxa assigned") +
                xlab("prevalence (No of remaining taxa)") +
                theme_bw() +
                theme(axis.text.x = element_text(angle=90, vjust=0.5))
        
        list(PC_assigned, tr)
        
}
# --



# --
#######################################
#### plot_correlations_abundance_prev_sparsity
#######################################
# physeq
# col = here a string of a taxonomic level

plot_correlations_abundance_prev_sparsity <- function(physeq, col = NULL){ 
        
        
        if (phyloseq::taxa_are_rows(physeq)) {
                physeq <- t(physeq)
        }
        
        df_ab_prev <- data.frame(Taxon_No = 1:ntaxa(physeq), 
                                 total_counts = taxa_sums(physeq),
                                 prevalence = colSums(as(otu_table(physeq), "matrix") != 0),
                                 sparsity = colSums(as(otu_table(physeq), "matrix") == 0),
                                 mean_count_nonzero = apply(as(otu_table(physeq), "matrix"), 2, function(x){mean(x[x > 0])}),
                                 median_count_nonzero = apply(as(otu_table(physeq), "matrix"), 2, function(x){median(x[x > 0])})
                                 )
        
        
        df_ab_prev <- cbind(df_ab_prev, tax_table(physeq))
        
        nsamples <- nsamples(physeq)
        
        if (!is.null(col)) {
                
                CountOrder <- dplyr::group_by_(df_ab_prev, col) %>% dplyr::summarise(total_count_sum = sum(total_counts)) %>% dplyr::arrange(desc(total_count_sum))
                # - make sure that also NAs are plotted, i.e. where coloring taxonomic rank is NA -
                CountOrder[[col]] <- as.character(CountOrder[[col]])
                CountOrder[[col]][is.na(CountOrder[[col]])] <- "NA"
                df_ab_prev[[col]] <- as.character(df_ab_prev[[col]])
                df_ab_prev[[col]][is.na(df_ab_prev[[col]])] <- "NA"
                df_ab_prev[[col]] <- factor(df_ab_prev[[col]], levels = CountOrder[[col]], ordered = TRUE)
                
                if (nrow(CountOrder) <= 15) {
                        custom_colors <- QuantColors15[1:nrow(CountOrder)]
                        
                } else {
                        custom_colors <- viridis(nrow(CountOrder))
                        
                }
                
                names(custom_colors) <- levels(df_ab_prev[[col]])
        }
        
        
        
        Tr_ab <- ggplot(df_ab_prev, aes(x = Taxon_No, y = total_counts))
        if (is.null(col)) {
                Tr_ab <- Tr_ab + geom_point(col = cbPalette[2], size = 2, alpha = 0.95) 
        } else {
                Tr_ab <- Tr_ab + geom_point(aes_string(col = col), size = 2, alpha = 0.95)
                Tr_ab <- Tr_ab + scale_color_manual("", values = custom_colors)
                
        }
        Tr_ab <- Tr_ab +
                ylab("total counts of taxon (= taxa_sums())") +
                theme_bw()
        
        
        Tr_prev <- ggplot(df_ab_prev, aes(x = Taxon_No, y = prevalence))
        if (is.null(col)) {
                Tr_prev <- Tr_prev + geom_point(col = cbPalette[2], size = 2, alpha = 0.95) 
        } else {
                Tr_prev <- Tr_prev + geom_point(aes_string(col = col), size = 2, alpha = 0.95)
                Tr_prev <- Tr_prev + scale_color_manual(values = custom_colors)
        }
        Tr_prev <- Tr_prev +
                theme_bw()
        
        
        
        fit_prev_log10 <- lm(formula = prevalence ~ log10(total_counts), data = df_ab_prev)
        pval_prev_log10 <- lmp(fit_prev_log10)
    
        fit_prev_mean <- lm(formula = prevalence ~ mean_count_nonzero, data = df_ab_prev)
        pval_prev_mean <- lmp(fit_prev_mean)
        fit_prev_mean_log10 <- lm(formula = prevalence ~ log10(mean_count_nonzero), data = df_ab_prev)
        pval_prev_mean_log10 <- lmp(fit_prev_mean_log10)
        fit_prev_median <- lm(formula = prevalence ~ median_count_nonzero, data = df_ab_prev)
        pval_prev_median <- lmp(fit_prev_median)
        fit_prev_median_log10 <- lm(formula = prevalence ~ log10(median_count_nonzero), data = df_ab_prev)
        pval_prev_median_log10 <- lmp(fit_prev_median_log10)
        
        
        Tr_prev_vs_log10_ab <- ggplot(df_ab_prev, aes(x = total_counts, y = prevalence))
        if (is.null(col)) {
                Tr_prev_vs_log10_ab <- Tr_prev_vs_log10_ab + geom_point(col = cbPalette[2], size = 2, alpha = 0.95) 
        } else {
                Tr_prev_vs_log10_ab <- Tr_prev_vs_log10_ab + geom_point(aes_string(col = col), size = 2, alpha = 0.95)
                Tr_prev_vs_log10_ab <- Tr_prev_vs_log10_ab + scale_color_manual("", values = custom_colors)
        }
        Tr_prev_vs_log10_ab <- Tr_prev_vs_log10_ab +
                scale_x_log10() +
                scale_y_continuous(limits = c(-1, max(df_ab_prev$prevalence) + 5)) +
                geom_smooth(method = "lm") +
                # annotate("text", x = max(df_ab_prev$abundance)/5, y = 5, label = paste("R.square of lm: ", as.character(round(summary(fit_prev)$r.squared,4), sep = ""))) +
                xlab("total counts of taxon (= taxa_sums())") +
                ggtitle(paste("lm fit: p.val: ", format(pval_prev_log10, digits = 4), "R.square: ", as.character(round(summary(fit_prev_log10)$r.squared,4), sep = ""))) +
                theme_bw(12)
        
        
        
        Tr_prev_vs_log10_meanab <- ggplot(df_ab_prev, aes(x = mean_count_nonzero, y = prevalence))
        if (is.null(col)) {
                Tr_prev_vs_log10_meanab <- Tr_prev_vs_log10_meanab + geom_point(col = cbPalette[2], size = 2, alpha = 0.95) 
        } else {
                Tr_prev_vs_log10_meanab <- Tr_prev_vs_log10_meanab + geom_point(aes_string(col = col), size = 2, alpha = 0.95)
                Tr_prev_vs_log10_meanab <- Tr_prev_vs_log10_meanab + scale_color_manual("", values = custom_colors)
        }
        Tr_prev_vs_log10_meanab <- Tr_prev_vs_log10_meanab +
                scale_x_log10() +
                geom_smooth(method = "lm") +
                scale_y_continuous(limits = c(-1, max(df_ab_prev$prevalence) + 5)) +
                xlab("mean count of taxon in non-zero samples") +
                ggtitle(paste("lm fit: p.val: ", format(pval_prev_mean_log10, digits = 4), "R.square: ", as.character(round(summary(fit_prev_mean_log10)$r.squared,4), sep = ""))) +
                theme_bw(12)
        
        Tr_prev_vs_log10_medianab <- ggplot(df_ab_prev, aes(x = median_count_nonzero, y = prevalence))
        if (is.null(col)) {
                Tr_prev_vs_log10_medianab <- Tr_prev_vs_log10_medianab + geom_point(col = cbPalette[2], size = 2, alpha = 0.95) 
        } else {
                Tr_prev_vs_log10_medianab <- Tr_prev_vs_log10_medianab + geom_point(aes_string(col = col), size = 2, alpha = 0.95)
                Tr_prev_vs_log10_medianab <- Tr_prev_vs_log10_medianab + scale_color_manual("", values = custom_colors)
        }
        Tr_prev_vs_log10_medianab <- Tr_prev_vs_log10_medianab +
                scale_x_log10() +
                geom_smooth(method = "lm") +
                scale_y_continuous(limits = c(-1, max(df_ab_prev$prevalence) + 5)) +
                xlab("median count of taxon in non-zero samples") +
                ggtitle(paste("lm fit: p.val: ", format(pval_prev_median_log10, digits = 4), "R.square: ", as.character(round(summary(fit_prev_median_log10)$r.squared,4), sep = ""))) +
                theme_bw(12)
        
        
        
        
        fitlist <- list(fit_prev_log10 = fit_prev_log10, fit_prev_mean = fit_prev_mean, 
                        fit_prev_mean_log10 = fit_prev_mean_log10, fit_prev_median = fit_prev_median,
                        fit_prev_median_log10 = fit_prev_median_log10)
        
        out <- list(Tr_ab = Tr_ab, 
                    Tr_prev = Tr_prev, 
                    Tr_prev_vs_log10_ab = Tr_prev_vs_log10_ab, 
                    Tr_prev_vs_log10_meanab = Tr_prev_vs_log10_meanab,
                    Tr_prev_vs_log10_medianab = Tr_prev_vs_log10_medianab,
                    fitlist = fitlist)
        
}
# --



# --
#######################################
#### plot_ab_pev_distributions
#######################################


plot_ab_pev_distributions <- function(physeq, prevalence = 5) {
        
        if (phyloseq::taxa_are_rows(physeq)) {
                physeq <- t(physeq)
        }
        
        seqtab <- as(otu_table(physeq), "matrix")
        
        FinalNumbersSeq <- data.frame(Sequence = colnames(seqtab), InNumberSamples = colSums(seqtab != 0), TotalCounts = colSums(seqtab))
        # a look at it shows you that seqtab is ordered after total counts
        FinalNumbersSeq <- group_by(FinalNumbersSeq, InNumberSamples)
        CountDistribution <- dplyr::summarise(FinalNumbersSeq, No_ASVs = n(), TotalCounts = sum(TotalCounts))
        CountDistribution$CumSumUnique <- rev(cumsum(rev(CountDistribution$No_ASVs)))
        CountDistribution$CumPerCUnique <- rev(cumsum(rev(CountDistribution$No_ASVs/ncol(seqtab))))
        CountDistribution$CumSumTotal <- rev(cumsum(rev(CountDistribution$TotalCounts)))
        CountDistribution$CumPerCTotal <- rev(cumsum(rev(CountDistribution$TotalCounts/sum(colSums(seqtab)))))
        
        PCValue <- ceiling((prevalence/100)*dim(seqtab)[1]) # tells you in how many samples a SV must be present to meet the prevalence 
        
        # Diff <- CountDistribution$InNumberSamples - PCValue
        # index <- which.max(Diff[Diff<0]) + which.min(Diff[Diff>=0])
        index <- which(CountDistribution$InNumberSamples >= PCValue)[1]
        PCKeptAtPCValue <- CountDistribution$CumPerCTotal[index]
        SVskeptAtPCValue <- CountDistribution$CumPerCUnique[index]
        
        # The number of samples the SVs are present in
        Tr <- ggplot(CountDistribution, aes(x = InNumberSamples, y = No_ASVs))
        Tr <- Tr + geom_point(col = "#E69F00", size = 3) +
                xlab("prevalence") +
                ylab("number of taxa") +
                theme_bw() +
                theme(panel.grid.minor = element_blank(),
                      panel.grid.major.y = element_blank(),
                      panel.grid.major.x = element_line(color = "#999999", size = .15))
        
        
        In1Index <- which(CountDistribution$InNumberSamples == 1)
        if (length(In1Index) != 0) {
                Tr <- Tr + ggtitle(paste(CountDistribution$No_ASVs[In1Index], " of ", CountDistribution$CumSumUnique[1], " taxa (", round(100*CountDistribution$No_ASVs[In1Index]/CountDistribution$CumSumUnique[1], 1), " %)", " were only found in 1 sample", sep = ""))
        } 
        
        
        # Cumulative Percentage of SVs
        Tr1 <- ggplot(CountDistribution, aes(x = InNumberSamples, y = CumPerCUnique))
        Tr1 <- Tr1 + geom_point(col = "#E69F00", size = 3) +
                xlab("prevalence") +
                ylab("cumulative percentage of taxa") +
                theme_bw() +
                theme(panel.grid.minor = element_blank(),
                      panel.grid.major.y = element_blank(),
                      panel.grid.major.x = element_line(color = "#999999", size = .15))
        
        Tr1 <- Tr1 + 
                geom_hline(yintercept = SVskeptAtPCValue, lty =  "dashed") +
                geom_vline(xintercept = (prevalence/100)*dim(seqtab)[1], lty = 'dashed') +
                ggtitle(paste("prevalence ", prevalence, " % = ", round((prevalence/100)*dim(seqtab)[1],1), "; ", CountDistribution$CumSumUnique[index], " of ", CountDistribution$CumSumUnique[1], 
                              " taxa (", round(100*CountDistribution$CumSumUnique[index]/CountDistribution$CumSumUnique[1], 1), " %) have higher prevalence", sep = ""))
        
        
        Tr2 <- ggplot(CountDistribution, aes(x = InNumberSamples, y = TotalCounts))
        Tr2 <- Tr2 + geom_point(col = "#E69F00", size = 3) +
                xlab("prevalence") +
                ylab("total counts of taxa with given prevalence") +
                theme_bw() +
                theme(panel.grid.minor = element_blank(),
                      panel.grid.major.y = element_blank(),
                      panel.grid.major.x = element_line(color = "#999999", size = .15))
        
        Tr2 <- Tr2 + ggtitle(paste(CountDistribution$CumSumTotal[1] - CountDistribution$CumSumTotal[index], " of ",
                                   CountDistribution$CumSumTotal[1], " (", round(100*(CountDistribution$CumSumTotal[1] - CountDistribution$CumSumTotal[index])/CountDistribution$CumSumTotal[1], 2),
                                   " %) counts are from taxa present in less than ", round((prevalence/100)*dim(seqtab)[1],1), " samples.", sep = ""))
        
        Tr3 <- ggplot(CountDistribution, aes(x = InNumberSamples, y = CumPerCTotal))
        Tr3 <- Tr3 + geom_point(col = "#E69F00", size = 3) +
                xlab("prevalence") +
                ylab("cumulative percentage of counts") +
                theme_bw() +
                theme(panel.grid.minor = element_blank(),
                      panel.grid.major.y = element_line(color = "#999999", size = .15),
                      panel.grid.major.x = element_line(color = "#999999", size = .15))
        
        Tr3 <- Tr3 + 
                geom_hline(yintercept = PCKeptAtPCValue, lty =  "dashed") +
                geom_vline(xintercept = (prevalence/100)*dim(seqtab)[1], lty = 'dashed') +
                ggtitle(paste("prevalence ", prevalence, " % = ", round((prevalence/100)*dim(seqtab)[1],1), "; ", CountDistribution$CumSumTotal[index], " of ", CountDistribution$CumSumTotal[1], 
                              " counts (", round(100*CountDistribution$CumSumTotal[index]/CountDistribution$CumSumTotal[1], 1), " %) would remain", sep = ""))
        
        
        TrList <- list(Tr, Tr1, Tr2, Tr3, CountDistribution)
        
        TrList
        
}
# --



# --
#######################################
### boxplot_sampleSums
#######################################
# Output:
# a list of boxplots for the alpha diversity measures given (including resid plots), the plots indicate the p-values between all levels in the group variable

# NB: makes use of ggpubr::stat_compare_means that allows you to adjust more settings, e.g. plotting symbols instead of numbers (change label to "p.signif")

boxplot_sampleSums <- function(physeq, group_var, color_levels, shape, test = "t.test", p_adjust_method = "fdr",
                               symnum.args = list(cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")),
                               hide.ns = FALSE){
        
        DF <- cbind(sample_data(physeq), total_count = sample_sums(physeq))
        
        Tr <-  ggplot(DF, aes_string(x = group_var, y = "total_count", color = group_var)) 
        Tr <- Tr + 
                geom_boxplot(outlier.color = NA) +
                geom_jitter(aes_string(shape = shape), width = .2, height = 0, alpha = 0.65) +
                xlab("") +
                ylab("sample_sums()") +
                scale_color_manual("", values = color_levels) +
                theme_bw()
        if(is.null(shape)){
                Tr <- Tr + theme(legend.position = "none")
        }
        
        # - since you might have more than two levels in each plot you need to set the comparisons argument in stat_compare_means -
        group_fac <- factor(sample_data(physeq)[[group_var]])
        fac_levels <- levels(group_fac)
        
        comparisonList <- get_unique_facLevel_combinations(fac_levels)
        
        Tr <- Tr + ggpubr::stat_compare_means(comparisons = comparisonList, label = "p.signif", method = test, hide.ns = hide.ns)
        
        formulaa <- as.formula(paste("total_count ~", group_var, sep = " "))
        
        pVals <- compare_means(formula = formulaa, data = DF, method = test, p.adjust.method = p_adjust_method, symnum.args = symnum.args)
        
        list(pVals = pVals, Tr = Tr)
}
# --
