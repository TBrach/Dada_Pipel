#######################################
### calc_alphadiv_plusLmResids
#######################################
# Output:
# outlist: 
#    [[1]]: DF of alpha_diversity values, plus residuals of the values to a linear fit against total counts (sample_sums)
#    [[2]]: list of the fit objects derived from lm of alpha diversity values against sample_sums()

# NB: is based on estimate_richness function from phyloseq
# You could easily calculate Shannon, Chao1, Observed self, see alphaDiversityMeasures.Rmd
# estimate_richness itself uses functions from the vegan package

calc_alphadiv_plusLmResids <- function(physeq, measures = c("Observed", "Shannon")) {
        DF_alpha <- suppressWarnings(phyloseq::estimate_richness(physeq, measures = measures))
        
        DF_alpha$Total <- sample_sums(physeq)
        
        
        # because linear fits of alpha diversity measures to total counts are often highly significant, I add the residuals of these
        # linear fits. 

        fitlist <- list()
        ncol_df_alpha <- ncol(DF_alpha)
        
        for (i in 1:length(measures)){
                fit <- lm(DF_alpha[,measures[i]] ~ DF_alpha[,"Total"])
                fitlist[[i]] <- fit
                DF_alpha[, ncol_df_alpha + i] <- residuals(fit)
                colnames(DF_alpha)[ncol_df_alpha + i] <- paste(measures[i], "resids", sep = "_")
        }
        
        names(fitlist) <- measures
        
        
        DF_alpha <- data.frame(Sample = rownames(DF_alpha), DF_alpha, sample_data(physeq))
        
        outlist <- list(DF_alpha = DF_alpha, fitlist = fitlist)
}


#######################################
### calc_pVals_alphdiv
#######################################
# Output:
# a DF with the p-values of all pairwise comparisons between the levels in group_var for the different alpha diversity measures
# NB: uses ggpubr: compare_means

calc_pVals_alphdiv <- function(DF_alpha, measures, group, test = "t.test", 
                               symnum.args = list(cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")),
                               p.adjust.method = "BH"){
        
        y_columns <- which(colnames(DF_alpha) %in% c(measures, paste0(measures, "_resids", sep = "")))
        ttestList <- list()
        for (i in 1:length(y_columns)){
                comparison <- as.formula(paste(colnames(DF_alpha)[y_columns[i]], " ~ ", group, sep = ""))
                ttestList[[i]] <- ggpubr::compare_means(formula = comparison, data = DF_alpha, method = test, p.adjust.method = p.adjust.method, symnum.args = symnum.args)
        }
        alpha_div_ttests <- do.call("rbind", ttestList) %>% arrange(.y.)
}


#######################################
### boxplots_alphdiv
#######################################
# Output:
# a list of boxplots for the alpha diversity measures given (including resid plots), the plots indicate the p-values between all levels in the group variable

# NB: makes use of ggpubr::stat_compare_means that allows you to adjust more settings, e.g. plotting symbols instead of numbers (change label to "p.signif")

boxplots_alphdiv <- function(DF_alpha, measures, group, shape, color_levels, test = "t.test"){
        
        boxplotList <- list()
        group_var_levels <- levels(DF_alpha[[group]])
        
        y_columns <- which(colnames(DF_alpha) %in% c(measures, paste0(measures, "_resids", sep = "")))
        
        # in case you have more than two levels in group_var_levels
        fac_levels_num <- setNames(seq_along(group_var_levels), group_var_levels) 
        i_s <- outer(fac_levels_num, fac_levels_num, function(ivec, jvec){
                sapply(seq_along(jvec), function(x){
                        i <- jvec[x]
                })
        })
        j_s <- outer(fac_levels_num, fac_levels_num, function(ivec, jvec){
                sapply(seq_along(ivec), function(x){
                        j <- ivec[x]
                })
        })
        i_s <- i_s[lower.tri(i_s)]
        j_s <- j_s[lower.tri(j_s)]
        
        comparisonList <- list()
        for (k in seq_along(i_s)){
                comparisonList[[k]] <- c(group_var_levels[i_s[k]], group_var_levels[j_s[k]])
        }
        
        for (i in 1:length(y_columns)){
                aes_map <- aes_string(x = group, y = colnames(DF_alpha)[y_columns[i]], color = group)
                Tr <-  ggplot(DF_alpha, aes_map) + 
                        geom_boxplot(na.rm = TRUE, outlier.color = NA) +
                        geom_jitter(aes_string(shape = shape), width = .2, height = 0, alpha = 0.65) +
                        xlab("") +
                        scale_color_manual("", values = color_levels) +
                        theme_bw()
                Tr <- Tr + stat_compare_means(comparisons = comparisonList, method = "t.test", label = "p.format") # "p.signif"
                boxplotList[[i]] <- Tr
                names(boxplotList)[i] <- colnames(DF_alpha)[y_columns[i]]
        }
        boxplotList
}

#######################################
### lmPlots_alphdiv
#######################################
# Output:
# a list of linear fit plots for the alpha diversity measures given against Total = sample_sums()

# NB: requires the lmp function!

lmPlots_alphdiv <- function(DF_alpha, lm_fitlist, measures, group, shape, color_levels, test = "t.test"){
        
        lm_plotList <- list()
        for (i in 1:length(measures)){
                aes_map <- aes_string(x = "Total", y = measures[i], color = group, shape = shape)
                Tr <-  ggplot(DF_alpha, aes_map) + 
                        geom_point(na.rm = TRUE) +
                        xlab("total counts (sample_sums())") +
                        scale_color_manual("", values = color_levels) +
                        theme_bw()
                fit <- lm_fitlist[[measures[i]]]
                adjR2 <- round(summary(fit)$adj.r.squared,3)
                if (adjR2 != 0){
                        pfit <- lmp(fit)
                } else {
                        pfit <- NA
                }
                Tr <- Tr + geom_smooth(aes_string(x = "Total", y = measures[i]), method = "lm", se = TRUE, inherit.aes = F)
                
                Tr <- Tr + ggtitle(paste("lm fit: p-value: ", format(pfit, digits = 4), ", adj.r2: ", adjR2, sep = ""))
                lm_plotList[[i]] <- Tr
                names(lm_plotList)[i] <- measures[i]
                
        }
        lm_plotList
}



