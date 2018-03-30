#######################################
### gm_own: calculate geometric mean
#######################################
# see commented below, this function comes from 
# <http://stackoverflow.com/questions/2602583/geometric-mean-is-there-a-built-in>

## Input:
# x numeric vector
# na.rm: if FALSE you get NA as soon as an NA is in your data, if TRUE the NA get basically treated as 0 (but NOTE these zeros always count also when zero.count = FALSE)
# zeros.count: This is IMPORTANT, if TRUE 0s count, so if x contains a 0 the GM will be lower than when zero.count = FALSE, in 
# which case the gm is calculated for the same vector in which all 0 have been removed.
## Output:
# the geometric mean of x
gm_own = function(x, na.rm=FALSE, zeros.count = TRUE){
        if(any(x < 0, na.rm = TRUE)){
                return(NaN)
        }
        if(zeros.count){
                exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
        } else {
                exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x[x > 0]))
        }
}

#######################################
### adj_LS
#######################################
# implementation of DESeq2 library size similar to estimateSizeFactorsForMatrix
# in difference to adjust_LS (obsolete) ignore.zero.ratios is always TRUE (so 0 ratios are always ignored)

## Input:
# physeq
# zeros.count: zeros.count of gm_own, if TRUE 0 will be considered when calculating the geometric mean
# if FALSE not and thus the geometric means will be bigger (see gm_own)
# percentile: the percentile to select the size factor SF of a sample based on its count ratios to the reference sample. In
# DESeq percentile = 50, i.e. stats::median is used. 
# plots: if TRUE SFs and plots will be given in an output list, otherwise just the adjusted physeq is returned


adj_LS <- function(physeq, zeros.count = FALSE, percentile = 50, plots = FALSE)  {
        
        # ---- Step 1: Calculate Geometric mean for each taxa over all samples ------
        # NB: these GM is basically the reference sample
        if(taxa_are_rows(physeq)){
                GM <- apply(otu_table(physeq), 1, gm_own, zeros.count = zeros.count)   
        } else {
                GM <- apply(otu_table(physeq), 2, gm_own, zeros.count = zeros.count) 
        }
        
        # ---- Step 2: Calculate Size factors --------
        
        # NB: x/y = exp(log(x) - log(y))
        
        if (taxa_are_rows(physeq)) {
                SFs <- apply(as(otu_table(physeq), "matrix"), 2, function(sample_cnts){exp(quantile((log(sample_cnts)-log(GM))[sample_cnts > 0], probs = percentile/100, na.rm = T))})
                SFs <- SFs/exp(mean(log(SFs)))
        } else {
                SFs <- apply(as(otu_table(physeq), "matrix"), 1, function(sample_cnts){exp(quantile((log(sample_cnts)-log(GM))[sample_cnts > 0], probs = percentile/100, na.rm = T))})
                SFs <- SFs/exp(mean(log(SFs)))
        } 
        
        
        if(min(SFs) == 0) { warning("in at least one sample the Size Factor was 0!") }
        
        
        # --- 3: calculate the new counts and put into a physeq object
        
        if(taxa_are_rows(physeq)){
                if (!identical(colnames(otu_table(physeq)), names(SFs))) {stop("names SFs do not fit to physeq")}
                Mat <- sweep(otu_table(physeq),2,SFs, "/")
                phynew <- physeq
                otu_table(phynew) <- otu_table(Mat, taxa_are_rows = TRUE)
        } else {
                if (!identical(rownames(otu_table(physeq)), names(SFs))) {stop("names SFs do not fit to physeq")}
                Mat <- sweep(otu_table(physeq),1,SFs, "/") 
                phynew <- physeq
                otu_table(phynew) <- otu_table(Mat, taxa_are_rows = FALSE)
        }
        
        if (plots){
                
                cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
                
                
                # compare calculated SFs to library sizes of the samples
                if(!identical(names(SFs), names(sample_sums(physeq)))){warning("Probably some Mix Up CHECK!")}
                comp <- data.frame(Sample = names(SFs), TotAmps = sample_sums(physeq), SFs = SFs)
                comp$SFsNormed <- comp$SFs/median(comp$SFs)
                # comp$TotAmpsNormed <- comp$TotAmps/mean(comp$TotAmps)
                comp$TotAmpsNormed <- comp$TotAmps/median(comp$TotAmps)
                comp <- dplyr::arrange(comp, desc(TotAmpsNormed))
                comp$Sample <- as.factor(comp$Sample)
                LevelsWant <- as.character(comp$Sample)
                for(i in 1:length(LevelsWant)){
                        comp$Sample <- relevel(comp$Sample, LevelsWant[i])
                }
                
                comp <- comp[c(1,4,5)]
                names(comp)[c(2,3)] <- c("SizeFactor_DESeq", "TotalAmplicons_relAb")
                comp <- tidyr::gather(comp, key = Corrector, value = NormedValue, -Sample)
                Tr <- ggplot(comp, aes(x = Sample, y = NormedValue, color = Corrector))
                Tr <- Tr + geom_point(size = 2) +
                        xlab("") +
                        ylab("correction value (normalized to median)") +
                        ggtitle(paste("Median SF: ", round(median(SFs),3), " Median TotAmp: ", round(median(sample_sums(physeq)),3), sep = "")) +
                        scale_color_manual(values = cbPalette[c(4,2)]) +
                        theme_bw() +
                        theme(panel.grid.minor = element_blank(),
                              panel.grid.major.y = element_blank(),
                              panel.grid.major.x = element_line(color = "#999999", size = .15),
                              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
                              legend.title = element_blank())
                
                # -- 3c: save Histograms of the SFs and sample_sums
                histo2 <- function(x, xtitle, gtitle) {
                        x <- data.frame(x = x)
                        Tr <- ggplot(x, aes(x = x))
                        Tr <- Tr + geom_histogram(binwidth = diff(range(x))/60, col = "black", fill = "#E69F00") +
                                geom_rug() +
                                geom_vline(xintercept = median(x$x), col = "#009E73", size = 1) +
                                ylab("No Samples") + 
                                xlab(xtitle) +
                                ggtitle(gtitle) +
                                theme_bw() + 
                                theme(panel.grid.minor = element_blank(),
                                      panel.grid.major.y = element_blank(),
                                      panel.grid.major.x = element_line(color = "#999999", size = .15))
                }
                
                Tr2 <- histo2(SFs, xtitle = "Size Factors", gtitle = "Size Factors a la DESeq")
                Tr3 <- histo2(sample_sums(physeq), xtitle = "Total Amplicons", gtitle = "Size Factors a la relative abundance")
                
                List <- list(Physeq = phynew, SFs = SFs, SizeFactorCompare = Tr, SFHistos = list(Tr2, Tr3), RefSample = GM)
                
                
        } else {
                
                list(physeq = phynew, SFs = SFs)
                
        }
}


#######################################
#### plot_abundance_prev_filter
#######################################
# df_ab_prev: data frame with "SV_ID", "total_counts_of_ASV", "prevalence", "sparsity", "mean_count_nonzero",
# "median_count_nonzero" plus tax_table

plot_abundance_prev_filter <- function(physeq, prevalence, taxa_sums_quantile){ 
        
        df_ab_prev <- data.frame(SV_ID = 1:ntaxa(physeq), 
                                 total_counts_of_ASV = taxa_sums(physeq),
                                 prevalence = colSums(as(otu_table(physeq), "matrix") != 0),
                                 sparsity = colSums(as(otu_table(physeq), "matrix") == 0), 
                                 mean_count_nonzero = apply(as(otu_table(physeq), "matrix"), 2, function(x){mean(x[x > 0])}),
                                 median_count_nonzero = apply(as(otu_table(physeq), "matrix"), 2, function(x){median(x[x > 0])}))
        
        df_ab_prev <- cbind(df_ab_prev, tax_table(physeq))
        
        prev_thresh <- (prevalence/100)*nsamples(physeq)
        abund_thresh <- quantile(taxa_sums(physeq), probs = taxa_sums_quantile/100)
        
        df_ab_prev_filt <- dplyr::filter(df_ab_prev, prevalence > prev_thresh | total_counts_of_ASV > abund_thresh)
        
        no_samples <- df_ab_prev$prevalence[1] + df_ab_prev$sparsity[1]
        shade_df <- data.frame(total_counts_of_ASV = 0, prevalence = 0)
        
        
        Tr_prev_vs_log10_ab <- ggplot(df_ab_prev, aes(x = total_counts_of_ASV, y = prevalence))
        Tr_prev_vs_log10_ab <- Tr_prev_vs_log10_ab +
                geom_point(col = cbPalette[2], size = 2, alpha = 0.7) +
                scale_x_log10() +
                geom_rect(data = shade_df, xmin = -Inf, xmax = log10(abund_thresh), ymin = -Inf, ymax = prev_thresh, fill = "#660033", alpha = 0.4) +
                geom_hline(yintercept = (prevalence/100)*nsamples(physeq), col = cbPalette[1], lty = "dashed") +
                geom_vline(xintercept = quantile(taxa_sums(physeq), probs = taxa_sums_quantile/100), col = cbPalette[1], lty = "dashed") +
                xlab("total counts of ASVs (taxa_sums())") + 
                theme_bw(12) +
                ggtitle(paste(nrow(df_ab_prev_filt), " of ", nrow(df_ab_prev), " ASVs (", round(100*nrow(df_ab_prev_filt)/nrow(df_ab_prev), 1),
                              " %) and ", round(sum(df_ab_prev_filt$total_counts_of_ASV)), " of ", round(sum(df_ab_prev$total_counts_of_ASV)), " counts (",
                              round((sum(df_ab_prev_filt$total_counts_of_ASV)/sum(df_ab_prev$total_counts_of_ASV))*100, 1), " %) remain", sep = ""))
        
        
        Tr_prev_vs_log10_ab_col <- ggplot(df_ab_prev, aes(x = total_counts_of_ASV, y = prevalence))
        Tr_prev_vs_log10_ab_col <- Tr_prev_vs_log10_ab_col +
                geom_point(aes(col = Phylum), size = 2, alpha = 0.7) +
                scale_x_log10() +
                geom_rect(data = shade_df, xmin = -Inf, xmax = log10(abund_thresh), ymin = -Inf, ymax = prev_thresh, fill = "#660033", alpha = 0.4) +
                geom_hline(yintercept = (prevalence/100)*nsamples(physeq), col = cbPalette[1], lty = "dashed") +
                geom_vline(xintercept = quantile(taxa_sums(physeq), probs = taxa_sums_quantile/100), col = cbPalette[1], lty = "dashed") +
                xlab("total counts of ASVs (taxa_sums())") +
                facet_wrap(~Phylum) +
                theme_bw(12) +
                theme(legend.position = "none")
        
        
        # phylum_df <- df_ab_prev[, c("Phylum", "total_counts_of_ASV", "prevalence")]
        # phylum_df <- group_by(phylum_df, Phylum)
        # phylum_df <- dplyr::summarise(phylum_df, SVs = n(), abundance = round(sum(total_counts_of_ASV)))
        # phylum_df_filt <- df_ab_prev_filt[, c("Phylum", "total_counts_of_ASV", "prevalence")]
        # phylum_df_filt <- group_by(phylum_df_filt, Phylum)
        # phylum_df_filt <- dplyr::summarise(phylum_df_filt, SVs = n(), abundance = round(sum(total_counts_of_ASV)))
        # phylum_df_summary <- merge(phylum_df, phylum_df_filt, by = "Phylum")
        # colnames(phylum_df_summary) <- c("Phylum", "SVs_before", "abundance_before", "SVs_after", 'abundance_after')
        # phylum_df_summary <- mutate(phylum_df_summary, SV_r_PC = round(100*SVs_after/SVs_before, 1), abundance_r_PC = round(100*abundance_after/abundance_before, 1),
        #                             SV_PC = round(100*SVs_after/sum(SVs_after), 1), abundance_PC = round(100*abundance_after/sum(abundance_after), 1))
        
        Before <- dplyr::summarise(group_by(df_ab_prev, Phylum), ASVs_bef = n(), PC_ASV_bef = round(100*n()/nrow(df_ab_prev),1), PC_total_pre_bef = round(100*sum(prevalence)/sum(df_ab_prev$prevalence), 1),
                                   PC_total_ab_bef = round(100*sum(total_counts_of_ASV)/sum(df_ab_prev$total_counts_of_ASV), 1), 
                                   mean_pre_bef_inPC = round(100*mean(prevalence)/no_samples, 1),
                                   mean_tot_ab_bef = round(mean(total_counts_of_ASV)), 
                                   mean_mean_ab_nonzero_bef = round(mean(mean_count_nonzero)),
                                   med_med_ab_nonzero_bef = round(median(median_count_nonzero)),
                                   total_ab_bef = sum(total_counts_of_ASV))
        
        After <- dplyr::summarise(group_by(df_ab_prev_filt, Phylum), ASVs_aft = n(), PC_ASV_aft = round(100*n()/nrow(df_ab_prev_filt),1), PC_total_pre_aft = round(100*sum(prevalence)/sum(df_ab_prev_filt$prevalence), 1),
                                  PC_total_ab_aft = round(100*sum(total_counts_of_ASV)/sum(df_ab_prev_filt$total_counts_of_ASV), 1), 
                                  mean_pre_aft_inPC = round(100*mean(prevalence)/no_samples, 1),
                                  mean_tot_ab_aft = round(mean(total_counts_of_ASV)),
                                  mean_mean_ab_nonzero_aft = round(mean(mean_count_nonzero)),
                                  med_med_ab_nonzero_aft = round(median(median_count_nonzero)),
                                  total_ab_aft = sum(total_counts_of_ASV))
        
        Before_total <- dplyr::summarise(df_ab_prev, ASVs_bef = n(), PC_ASV_bef = round(100*n()/nrow(df_ab_prev),1), PC_total_pre_bef = round(100*sum(prevalence)/sum(df_ab_prev$prevalence), 1),
                                         PC_total_ab_bef = round(100*sum(total_counts_of_ASV)/sum(df_ab_prev$total_counts_of_ASV), 1), 
                                         mean_pre_bef_inPC = round(100*mean(prevalence)/no_samples, 1),
                                         mean_tot_ab_bef = round(mean(total_counts_of_ASV)), 
                                         mean_mean_ab_nonzero_bef = round(mean(mean_count_nonzero)),
                                         med_med_ab_nonzero_bef = round(median(median_count_nonzero)),
                                         total_ab_bef = sum(total_counts_of_ASV))
        
        Before_total <- data.frame(Phylum = "Total", Before_total)
        
        After_total <- dplyr::summarise(df_ab_prev_filt, ASVs_aft = n(), PC_ASV_aft = round(100*n()/nrow(df_ab_prev_filt),1), PC_total_pre_aft = round(100*sum(prevalence)/sum(df_ab_prev_filt$prevalence), 1),
                                        PC_total_ab_aft = round(100*sum(total_counts_of_ASV)/sum(df_ab_prev_filt$total_counts_of_ASV), 1), 
                                        mean_pre_aft_inPC = round(100*mean(prevalence)/no_samples, 1),
                                        mean_tot_ab_aft = round(mean(total_counts_of_ASV)),
                                        mean_mean_ab_nonzero_aft = round(mean(mean_count_nonzero)),
                                        med_med_ab_nonzero_aft = round(median(median_count_nonzero)),
                                        total_ab_aft = sum(total_counts_of_ASV))
        
        After_total <- data.frame(Phylum = "Total", After_total)
        
        Before <- rbind(Before, Before_total)
        
        After <- rbind(After, After_total)
        
        
        Merged <- merge(Before, After, by = "Phylum", all = TRUE, sort = FALSE)
        
        Merged <- dplyr::mutate(Merged, PC_ASV_rem = round(100*ASVs_aft/ASVs_bef, 1), PC_ab_rem = round(100*total_ab_aft/total_ab_bef, 1))
        
        Merged <- dplyr::select(Merged, 1, 20, 21, 2, 11, 3, 12, 4, 13, 5, 14, 6, 15, 7, 16, 8, 17, 9, 18, 10, 19)
        
        # in case some phyla have been kicked out: you need to reorder like it was "Before" to have total at the last place
        Merged[is.na(Merged)] <- 0 # Phyla that have been removed are NA
        colnames(Merged)[20:21] <- c("tot_ab_bef", "tot_ab_aft")
        Merged[, 20] <- round(Merged[, 20])
        Merged[, 21] <- round(Merged[, 21])
        
        #Merged$tot_ab_bef <- round(Before$total_ab_bef)
        #Merged$tot_ab_aft <- round(After$total_ab_aft)
        
        out <- list(Tr_prev_vs_log10_ab = Tr_prev_vs_log10_ab,
                    Tr_prev_vs_log10_ab_col = Tr_prev_vs_log10_ab_col,
                    phylum_df_summary = Merged)
        
}
