cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#######################################
### lmp (to get p_value from lm fit)
#######################################
# needed to get the p-value from-linear fit objects (from stackoverflow)
lmp <- function (modelobject) {
        if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
        f <- summary(modelobject)$fstatistic
        p <- pf(f[1],f[2],f[3],lower.tail=F)
        attributes(p) <- NULL
        return(p)
}


#######################################
### FUNCTION: get_assignemnt_distribution
#######################################

# Function to determine how many sequences could not be assigned using a given minBoot value
## Input
# taxa: the output of the assignTaxonomy (dada2) command. A matrix with nrow = number of sequences and ncol usually 6, the taxonomic levels. When 
# the minBoot criterion was not fulfilled at a taxonomic level, assignTaxonomy assigns an NA. 
## Output
# assignment_distribution: data frame that states for each taxonomic level in taxa how many SVs have been assigned 

get_assignemnt_distribution <- function(taxa){
        
        countNA <- apply(taxa, 2, function(x){sum(is.na(x))})
        countNonNA <- apply(taxa, 2, function(x){sum(!is.na(x))})
        ambiguous <- apply(taxa, 2, function(x){sum(grepl("/", x))})
        total <- countNA + countNonNA
        assignment_distribution <- data.frame(assigned = countNonNA, 
                                              assigned_unamb = countNonNA - ambiguous,
                                              total = total,
                                              PC_assigned = round(100*countNonNA/total, 1),
                                              PC_assigned_unamb = round(100*(countNonNA - ambiguous)/total, 1))
        
}



#######################################
### FUNCTION: check_assignment_vs_abundance
#######################################
# checks if there is a trend for better assignment for more abundant SVs

check_assignment_vs_abundance <- function(taxa, seqtab, abundanceQuantiles = seq(0, 90, by = 10)){
        
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
                scale_x_discrete(breaks = abundanceQuantiles, labels = paste(abundanceQuantiles, " (", round(abQuantiles), ", ", No_ASVs, ")", sep = "")) +
                ylab("percentage of ASVs assigned") +
                xlab("total counts quantile (total counts filter, No ASVs remaining)") +
                theme_bw(12) +
                theme(axis.text.x = element_text(angle=90, vjust=0.5))
        
        list(PC_assigned, tr)
        
}


#######################################
### FUNCTION: check_assignment_vs_prevalence
#######################################
# checks if there is a trend for better assignment for more prealent SVs

check_assignment_vs_prevalence <- function(taxa, seqtab, prevalences = seq(0, 90, by = 10)){
        
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
                ylab("percentage of ASVs assigned") +
                xlab("prevalence (No ASVs remaining)") +
                theme_bw(12) +
                theme(axis.text.x = element_text(angle=90, vjust=0.5))
        
        list(PC_assigned, tr)
        
}



#######################################
#### plot_correlations_abundance_prev_sparsity
#######################################
# df_ab_prev: data frame with "ASV_ID", "total_counts_of_ASV", "prevalence", "sparsity", "mean_count_nonzero",
# "median_count_nonzero"
# outputs a list with different plots and fits

plot_correlations_abundance_prev_sparsity <- function(df_ab_prev, col = NULL){ # NB: you could color by Phylum for example
        
        nsamples <- df_ab_prev$prevalence[1] + df_ab_prev$sparsity[1]
        
        Tr_ab <- ggplot(df_ab_prev, aes(x = ASV_ID, y = total_counts_of_ASV))
        if (is.null(col)) {
                Tr_ab <- Tr_ab + geom_point(col = cbPalette[2], size = 2, alpha = 0.7) 
        } else {
                Tr_ab <- Tr_ab + geom_point(aes_string(col = col), size = 2, alpha = 0.7)
        }
        Tr_ab <- Tr_ab +
                ylab("total counts of ASV (= taxa_sums())") +
                theme_bw(12)
        
        
        Tr_prev <- ggplot(df_ab_prev, aes(x = ASV_ID, y = prevalence))
        if (is.null(col)) {
                Tr_prev <- Tr_prev + geom_point(col = cbPalette[2], size = 2, alpha = 0.7) 
        } else {
                Tr_prev <- Tr_prev + geom_point(aes_string(col = col), size = 2, alpha = 0.7)
        }
        Tr_prev <- Tr_prev +
                theme_bw(12)
        
        
        # - associations of total abundance and log10(total_counts_of_ASV) to sparsity/prevalence -
        # NB: turned out: log10(total_counts_of_ASV) is far better correlated with prevalence/sparsity, so I stick with the log fits
        # Also NB: since prevalence = constant - sparsity, a fit prevalence ~ abundance is equal to fit sparsity ~ abundance just with
        # opposite coefficient signs
        
        # fit_spar <- lm(formula = sparsity ~ total_counts_of_ASV, data = df_ab_prev)
        # pval_spar <- lmp(fit_spar)
        # fit_prev <- lm(formula = prevalence ~ total_counts_of_ASV, data = df_ab_prev)
        # pval_prev <- lmp(fit_prev)
        fit_prev_log10 <- lm(formula = prevalence ~ log10(total_counts_of_ASV), data = df_ab_prev)
        pval_prev_log10 <- lmp(fit_prev_log10)
        # fit_spar_log10 <- lm(formula = sparsity ~ log10(total_counts_of_ASV), data = df_ab_prev)
        # pval_spar_log10 <- lmp(fit_spar_log10)
        # identical(pval_prev_log10, pval_spar_log10) # TRUE
        
        # it comes natural that total abundance and prevalence/sparsity are correlated, but how about mean abundance in non zero samples
        # here I prepare all combinations, but stick to prevalence for now (comment sparsity out) since correlations are the same
        
        
        # fit_spar_mean <- lm(formula = sparsity ~ mean_count_nonzero, data = df_ab_prev)
        # pval_spar_mean <- lmp(fit_spar_mean)
        # fit_spar_mean_log10 <- lm(formula = sparsity ~ log10(mean_count_nonzero), data = df_ab_prev)
        # pval_spar_mean_log10 <- lmp(fit_spar_mean_log10)
        # fit_spar_median <- lm(formula = sparsity ~ median_count_nonzero, data = df_ab_prev)
        # pval_spar_median <- lmp(fit_spar_median)
        # fit_spar_median_log10 <- lm(formula = sparsity ~ log10(median_count_nonzero), data = df_ab_prev)
        # pval_spar_median_log10 <- lmp(fit_spar_median_log10)
        
        fit_prev_mean <- lm(formula = prevalence ~ mean_count_nonzero, data = df_ab_prev)
        pval_prev_mean <- lmp(fit_prev_mean)
        fit_prev_mean_log10 <- lm(formula = prevalence ~ log10(mean_count_nonzero), data = df_ab_prev)
        pval_prev_mean_log10 <- lmp(fit_prev_mean_log10)
        fit_prev_median <- lm(formula = prevalence ~ median_count_nonzero, data = df_ab_prev)
        pval_prev_median <- lmp(fit_prev_median)
        fit_prev_median_log10 <- lm(formula = prevalence ~ log10(median_count_nonzero), data = df_ab_prev)
        pval_prev_median_log10 <- lmp(fit_prev_median_log10)
        
        
        # Tr_ab_vs_prev <- ggplot(df_ab_prev, aes(x = total_counts_of_ASV, y = prevalence))
        # Tr_ab_vs_prev <- Tr_ab_vs_prev +
        #         geom_point(col = cbPalette[2], size = 2, alpha = 0.7) +
        #         geom_smooth(method = "lm") +
        #         scale_y_continuous(limits = c(-1, max(df_ab_prev$prevalence) + 5)) +
        #         # annotate("text", x = max(df_ab_prev$abundance)/5, y = 5, label = paste("R.square of lm: ", as.character(round(summary(fit_prev)$r.squared,4), sep = ""))) +
        #         xlab("total abundance (taxa_sums())") +
        #         ggtitle(paste("lm fit: p.val: ", format(pval_prev, digits = 4), "R.square: ", as.character(round(summary(fit_prev)$r.squared,4), sep = ""))) +
        # theme_bw(12)
        # Tr_ab_vs_prev_75Q <- Tr_ab_vs_prev + coord_cartesian(xlim = c(-5, quantile(df_ab_prev$abundance, .75)))
        
        Tr_prev_vs_log10_ab <- ggplot(df_ab_prev, aes(x = total_counts_of_ASV, y = prevalence))
        if (is.null(col)) {
                Tr_prev_vs_log10_ab <- Tr_prev_vs_log10_ab + geom_point(col = cbPalette[2], size = 2, alpha = 0.7) 
        } else {
                Tr_prev_vs_log10_ab <- Tr_prev_vs_log10_ab + geom_point(aes_string(col = col), size = 2, alpha = 0.7)
        }
        Tr_prev_vs_log10_ab <- Tr_prev_vs_log10_ab +
                scale_x_log10() +
                scale_y_continuous(limits = c(-1, max(df_ab_prev$prevalence) + 5)) +
                geom_smooth(method = "lm") +
                # annotate("text", x = max(df_ab_prev$abundance)/5, y = 5, label = paste("R.square of lm: ", as.character(round(summary(fit_prev)$r.squared,4), sep = ""))) +
                xlab("total counts of ASV (= taxa_sums())") +
                ggtitle(paste("lm fit: p.val: ", format(pval_prev_log10, digits = 4), "R.square: ", as.character(round(summary(fit_prev_log10)$r.squared,4), sep = ""))) +
                theme_bw(12)
        
        # NB: you could color or facet by phylum
        
        # Tr_prev_vs_log10_ab <- ggplot(df_ab_prev, aes(x = total_counts_of_ASV, y = prevalence))
        # Tr_prev_vs_log10_ab <- Tr_prev_vs_log10_ab +
        #         geom_point(aes(col = Phylum), size = 2, alpha = 0.7) +
        #         scale_x_log10() +
        #         geom_smooth(method = "lm") +
        #         # annotate("text", x = max(df_ab_prev$abundance)/5, y = 5, label = paste("R.square of lm: ", as.character(round(summary(fit_prev)$r.squared,4), sep = ""))) +
        #         xlab("total abundance (taxa_sums())") +
        #         ggtitle(paste("lm fit: p.val: ", format(pval_prev_log10, digits = 4), "R.square: ", as.character(round(summary(fit_prev_log10)$r.squared,4), sep = ""))) +
        #         theme_bw(12)
        
        
        # Tr_spar_vs_meanab <- ggplot(df_ab_prev, aes(x = mean_count_nonzero, y = sparsity))
        # Tr_spar_vs_meanab <- Tr_spar_vs_meanab +
        #         geom_point(col = cbPalette[2], size = 2, alpha = 0.7) +
        #         geom_smooth(method = "lm") +
        #         scale_y_continuous(limits = c(-1, nsamples + 5)) +
        #         xlab("SVs mean abundance in non-zero samples") +
        #         ggtitle(paste("lm fit: p.val: ", format(pval_spar_mean, digits = 4), "R.square: ", as.character(round(summary(fit_spar_mean)$r.squared,4), sep = ""))) +
        #         theme_bw(12)
        # Tr_spar_vs_meanab_75Q <- Tr_spar_vs_meanab + coord_cartesian(xlim = c(-5, quantile(df_ab_prev$mean_count_nonzero, .75)))
        # 
        
        # Tr_spar_vs_log10_meanab <- ggplot(df_ab_prev, aes(x = mean_count_nonzero, y = sparsity))
        # Tr_spar_vs_log10_meanab <- Tr_spar_vs_log10_meanab +
        #         geom_point(col = cbPalette[2], size = 2, alpha = 0.7) +
        #         scale_x_log10() +
        #         geom_smooth(method = "lm") +
        #         xlab("SVs mean abundance in non-zero samples") +
        #         ggtitle(paste("lm fit: p.val: ", format(pval_spar_mean_log10, digits = 4), "R.square: ", as.character(round(summary(fit_spar_mean_log10)$r.squared,4), sep = ""))) +
        #         theme_bw(12)
        
        Tr_prev_vs_log10_meanab <- ggplot(df_ab_prev, aes(x = mean_count_nonzero, y = prevalence))
        if (is.null(col)) {
                Tr_prev_vs_log10_meanab <- Tr_prev_vs_log10_meanab + geom_point(col = cbPalette[2], size = 2, alpha = 0.7) 
        } else {
                Tr_prev_vs_log10_meanab <- Tr_prev_vs_log10_meanab + geom_point(aes_string(col = col), size = 2, alpha = 0.7)
        }
        Tr_prev_vs_log10_meanab <- Tr_prev_vs_log10_meanab +
                scale_x_log10() +
                geom_smooth(method = "lm") +
                scale_y_continuous(limits = c(-1, max(df_ab_prev$prevalence) + 5)) +
                xlab("mean count of ASV in non-zero samples") +
                ggtitle(paste("lm fit: p.val: ", format(pval_prev_mean_log10, digits = 4), "R.square: ", as.character(round(summary(fit_prev_mean_log10)$r.squared,4), sep = ""))) +
                theme_bw(12)
        
        Tr_prev_vs_log10_medianab <- ggplot(df_ab_prev, aes(x = median_count_nonzero, y = prevalence))
        if (is.null(col)) {
                Tr_prev_vs_log10_medianab <- Tr_prev_vs_log10_medianab + geom_point(col = cbPalette[2], size = 2, alpha = 0.7) 
        } else {
                Tr_prev_vs_log10_medianab <- Tr_prev_vs_log10_medianab + geom_point(aes_string(col = col), size = 2, alpha = 0.7)
        }
        Tr_prev_vs_log10_medianab <- Tr_prev_vs_log10_medianab +
                scale_x_log10() +
                geom_smooth(method = "lm") +
                scale_y_continuous(limits = c(-1, max(df_ab_prev$prevalence) + 5)) +
                xlab("median count of ASV in non-zero samples") +
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


#######################################
#### plotSVdistributions
#######################################
## SORRY FOR THE TERMINOLOGY CONFUSION< from now on counts is used instead of amplicons in the plots
## REQUIRES dplyr
## Input:
# Seqtab: seqtab from dada2 wrapper
# prevalence: a percentage of samples (bw 0 and 100)
## Output: 
# TrList: List of four Trelis objects, plus distribution data.frame


plotSVdistributions <- function(seqtab, prevalence = 10) {
        
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
        
        return(TrList)
        
}