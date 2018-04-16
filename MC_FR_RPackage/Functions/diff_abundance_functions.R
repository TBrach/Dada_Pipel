pal <- function(col, border = "light gray", ...){
        n <- length(col)
        plot(0, 0, type="n", xlim = c(0, 1), ylim = c(0, 1),
             axes = FALSE, xlab = "", ylab = "", ...)
        rect(0:(n-1)/n, 0, 1:n/n, 1, col = col, border = border)
}


# tol1qualitative=c("#4477AA")
# tol2qualitative=c("#4477AA", "#CC6677")
# tol3qualitative=c("#4477AA", "#DDCC77", "#CC6677")
# tol4qualitative=c("#4477AA", "#117733", "#DDCC77", "#CC6677")
# tol5qualitative=c("#332288", "#88CCEE", "#117733", "#DDCC77", "#CC6677")
# tol6qualitative=c("#332288", "#88CCEE", "#117733", "#DDCC77", "#CC6677","#AA4499")
# tol7qualitative=c("#332288", "#88CCEE", "#44AA99", "#117733", "#DDCC77", "#CC6677","#AA4499")
# tol8qualitative=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677","#AA4499")
# tol9qualitative=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677", "#882255", "#AA4499")
# tol10qualitative=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#882255", "#AA4499")
# tol11qualitative=c("#332288", "#6699CC", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#882255", "#AA4499")
# tol12qualitative=c("#332288", "#6699CC", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#AA4466", "#882255", "#AA4499")
# tol14rainbow=c("#882E72", "#B178A6", "#D6C1DE", "#1965B0", "#5289C7", "#7BAFDE", "#4EB265", "#90C987", "#CAE0AB", "#F7EE55", "#F6C141", "#F1932D", "#E8601C", "#DC050C")
tol15rainbow=c("#114477", "#4477AA", "#77AADD", "#117755", "#44AA88", "#99CCBB", "#777711", "#AAAA44", "#DDDD77", "#771111", "#AA4444", "#DD7777", "#771144", "#AA4477", "#DD77AA")
#tol18rainbow=c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")
# ...and finally, the Paul Tol 21-color salute
#tol21rainbow= c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")

# make your own colorscheme to well distinguish 15 quantitiative colors
QuantColors15 <- tol15rainbow[c(1, 12, 4, 7, 13, 2, 11, 5, 8, 14, 3, 10, 6, 9, 15)]

# pal(QuantColors15)

#######################################
### get_unique_facLevel_combinations
#######################################
# often needed to set the comparisons attribute in stat_compare_means from ggplot, it excludes comparison with each other and also removes 2 to 3 vs 3 to 2

get_unique_facLevel_combinations <- function(fac_levels){
        
        fac_levels_num <- setNames(seq_along(fac_levels), fac_levels) 
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
                comparisonList[[k]] <- c(fac_levels[i_s[k]], fac_levels[j_s[k]])
        }
        
        comparisonList
        
}


#######################################
### get_overview_of_physeq##
#################
# The function came to live while doing simulations: it became clear that the count variation (both across samples and taxa) is
# usually higher in real data than in simulated data. 
# Input: just physeq
# Output: a DF summarising features of the abundance table. 

get_overview_of_physeq <- function(physeq){
        
        if(taxa_are_rows(physeq)) { physeq <- t(physeq) }
        # I prefer taxa_are_rows = FALSE so rows (= observations = samples), and colums = variables = taxa
        
        CT <- as(otu_table(physeq), "matrix") # taxa are columns
        CT_NA <- CT
        CT_NA[CT_NA == 0] <- NA
        # also get a relative abundance Table
        CT_RA <- CT/rowSums(CT)
        CT_RA_NA <- CT_RA
        CT_RA_NA[CT_RA_NA == 0] <- NA
        
        
        Overview <- data.frame(NoSamples = nsamples(physeq),
                               NoTaxa = ntaxa(physeq),
                               MedianSampleSum = round(median(sample_sums(physeq))),
                               MedianTaxaSum = round(median(taxa_sums(physeq))),
                               Sparsity = round(100*(sum(otu_table(physeq) == 0)/(ntaxa(physeq)*nsamples(physeq))), 3),
                               MaxCount = max(CT),
                               MedianTaxaSD = round(median(apply(CT, 2, sd), na.rm = TRUE), 3),
                               MedTaxaSDNoZ = round(median(apply(CT_NA, 2, sd, na.rm = TRUE), na.rm = TRUE), 3),
                               
                               # NB: SD correlates with taxa_sums (= colSums(CT))
                               # plot(colSums(CT), apply(CT, 2, sd))
                               # plot(colSums(CT_NA, na.rm = TRUE), apply(CT_NA, 2, sd, na.rm = TRUE))
                               # in principle you could correct by taxa_sums
                               # MedianTaxaSD_Cor = round(median(apply(CT*(max(colSums(CT))/colSums(CT)), 2, sd), na.rm = TRUE), 3),
                               # just keep in mind it is confounded
                               
                               # get same taxa variation for relative abundance table
                               MedianTaxaSD_RA = round(median(apply(CT_RA, 2, sd), na.rm = TRUE), 6),
                               MedTaxaSDNoZ_RA = round(median(apply(CT_RA_NA, 2, sd, na.rm = TRUE), na.rm = TRUE), 6),
                               
                               MedianSampleSD = round(median(apply(CT, 1, sd), na.rm = TRUE), 3),
                               MedianSampleSDNoZ = round(median(apply(CT_NA, 1, sd, na.rm = TRUE), na.rm = TRUE), 3),
                               # makes in principle only sense on relative abundance
                               MedianSampleSD_RA = round(median(apply(CT_RA, 1, sd), na.rm = TRUE), 6),
                               MedianSampleSDNoZ_RA = round(median(apply(CT_RA_NA, 1, sd, na.rm = TRUE), na.rm = TRUE), 3)
        )
        
}



#######################################
### make_heat_map_physeq##
#################
## Input:
# physeq object
# group_var: the name of the group_fac column in sample_data used to order the samples in the heat map
# max_abundance_for_color: if null the 90% quantile of count/relative abundance in the data is used. all counts above this value will be
# shown yellow in the heat map
# tax_order: character vector of the original taxon names in the order you want them shown. if NULL >> tax_order = taxa_names(physeq)
# tax_names: the names that will be used for the taxons, if Null Taxon_1, Taxon_2 and so on will be used. NB: renaming of course after
# ordering. 
# color_sample_names: if TRUE and if you have less than 7 levels, the sample names will be colored using cbPalette
# gradient_steps: the steps the blue to green to yellow will be distributed in the viridis gradient: 4 numbers ranging from 1e-14 to 1,
# see default, you might wanna try c(0.25, 0.5, 0.75, 1) as well
## Output:
# the heat map trellis

make_heat_map_physeq <- function(physeq, group_var, color_levels, max_abundance_for_color = NULL, tax_order = NULL,
                                 tax_names = NULL, color_sample_names = TRUE, gradient_steps = c(0.15, 0.3, 0.45, 1)){
        
        gradient_steps <- c(0, 1e-14, gradient_steps)
        
        if(taxa_are_rows(physeq)) { physeq <- t(physeq) }
        
        if (!is.factor(sample_data(physeq)[[group_var]])) {sample_data(physeq)[[group_var]] <- as.factor(sample_data(physeq)[[group_var]])}
        
        
        CT <- as(otu_table(physeq), "matrix") # taxa are columns
        DF_CT <- as.data.frame(t(CT))
        
        if (is.null(tax_order)){tax_order <- taxa_names(physeq)}
        
        if (length(tax_order) != nrow(DF_CT) || !all(tax_order %in% rownames(DF_CT))) {stop("the given tax_order must contain all taxa_names(physeq). If you want to subset,
                                                                                            do it before with prune_taxa")}
        DF_CT <- DF_CT[tax_order, ]
        
        if (is.null(tax_names)) {tax_names <- paste("Taxon_", 1:nrow(DF_CT), sep = "")}
        
        if (length(tax_names) != nrow(DF_CT) ) {stop("the given tax_names must be a character vector of length = ntaxa(physeq)")}
        
        
        rownames(DF_CT) <- make.unique(tax_names)
        DF_CT$Taxa <- rownames(DF_CT)
        DF_CT <- tidyr::gather(DF_CT, key = Sample , value = Count, -Taxa)
        DF_CT$Taxa <- factor(DF_CT$Taxa, levels = rev(unique(DF_CT$Taxa)), ordered = TRUE)
        LookUpDF <- data.frame(Sample = sample_names(physeq), Group = sample_data(physeq)[[group_var]])
        LookUpDF <- LookUpDF[order(match(LookUpDF$Group, levels(LookUpDF$Group))), ]
        DF_CT$Sample <- factor(DF_CT$Sample, levels = LookUpDF$Sample, ordered = TRUE)
        
        if (is.null(max_abundance_for_color)) {
                max_abundance_for_color <- quantile(DF_CT$Count, .9)
        }
        
        if (max_abundance_for_color == 0) {max_abundance_for_color <- max(DF_CT$Count)}
        
        colxaxis <- color_levels[LookUpDF$Group]
        
        
        hmTr <- ggplot(DF_CT, aes(x = Sample, y = Taxa, fill = Count))
        hmTr <- hmTr + 
                geom_raster() + 
                scale_fill_gradientn("", limits = c(0, max_abundance_for_color), colors = c("red", viridis(5)), values = gradient_steps, oob = squish) +
                scale_x_discrete(position = "top") +
                #coord_equal() +
                labs(x=NULL, y=NULL) +
                theme_tufte(base_family = "Helvetica") +
                theme(axis.ticks=element_blank(),
                      axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0,
                                                 colour = colxaxis))
        hmTr
}



#######################################
### make_heat_map_physeq_levels##
#################
# same heat map plots as from make_heat_map_physeq but for each combination of levels in your physeq
## Input:
# physeq object
# group_var: the name of the group_fac column in sample_data used to order the samples in the heat map
# max_abundance_for_color: if null the 90% percentile count/relative abundance in the data is used. all counts above this value will be
# shown yellow in the heat map
# tax_order: character vector of the original taxon names in the order you want them shown. if NULL >> tax_order = taxa_names(physeq)
# tax_names: the names that will be used for the taxons, if Null Taxon_1, Taxon_2 and so on will be used. NB: renaming of course after
# ordering. 
# color_sample_names: if TRUE and if you have less than 7 levels, the sample names will be colored using cbPalette
# gradient_steps: the steps the blue to green to yellow will be distributed in the viridis gradient: 4 numbers ranging from 1e-14 to 1,
# see default, you might wanna try c(0.25, 0.5, 0.75, 1) as well
# Output: list of heat maps for each level combination

make_heat_map_physeq_levels <- function(physeq, group_var, color_levels, max_abundance_for_color = NULL, tax_order = NULL,
                                        tax_names = NULL, color_sample_names = TRUE, gradient_steps = c(0.15, 0.3, 0.45, 1)){
        
        gradient_steps <- c(0, 1e-14, gradient_steps)
        
        if(taxa_are_rows(physeq)) { physeq <- t(physeq) }
        
        if (!is.factor(sample_data(physeq)[[group_var]])) {sample_data(physeq)[[group_var]] <- as.factor(sample_data(physeq)[[group_var]])}
        
        CT <- as(otu_table(physeq), "matrix") # taxa are columns
        DF_CT <- as.data.frame(t(CT))
        
        if (is.null(tax_order)){tax_order <- taxa_names(physeq)}
        
        if (length(tax_order) != nrow(DF_CT) || !all(tax_order %in% rownames(DF_CT))) {stop("the given tax_order must contain all taxa_names(physeq). If you want to subset,
                                                                                            do it before with prune_taxa")}
        DF_CT <- DF_CT[tax_order, ]
        
        if (is.null(tax_names)) {tax_names <- paste("Taxon_", 1:nrow(DF_CT), sep = "")}
        
        if (length(tax_names) != nrow(DF_CT) ) {stop("the given tax_names must be a character vector of length = ntaxa(physeq)")}
        
        tax_names[is.na(tax_names)] <- "NA"
        rownames(DF_CT) <- make.unique(tax_names)
        DF_CT$Taxa <- rownames(DF_CT)
        DF_CT <- tidyr::gather(DF_CT, key = Sample , value = Count, -Taxa)
        DF_CT$Taxa <- factor(DF_CT$Taxa, levels = rev(unique(DF_CT$Taxa)), ordered = TRUE)
        LookUpDF <- data.frame(Sample = sample_names(physeq), Group = sample_data(physeq)[[group_var]])
        LookUpDF <- LookUpDF[order(match(LookUpDF$Group, levels(LookUpDF$Group))), ]
        DF_CT$Sample <- factor(DF_CT$Sample, levels = LookUpDF$Sample, ordered = TRUE)
        DF_CT$Group <- LookUpDF$Group[match(DF_CT$Sample, LookUpDF$Sample)]
        
        
        group_fac <- factor(sample_data(physeq)[[group_var]])
        
        fac_levels <- levels(group_fac)
        
        # -- get the level combis --
        fac_levels_num <- setNames(seq_along(fac_levels), fac_levels) 
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
        # ----
        
        plot_list <- vector("list", length = length(i_s))
        
        
        
        for (k in seq_along(i_s)){
                i <- i_s[k]
                j <- j_s[k]
                
                group_fac_current <- droplevels(group_fac[as.numeric(group_fac) %in% c(i, j)])
                DF_CT_current <- filter(DF_CT, Group %in% group_fac_current)
                DF_CT_current$Group <- factor(DF_CT_current$Group, levels = c(fac_levels[i], fac_levels[j]), ordered = TRUE)
                LookUpDF_current <- LookUpDF[LookUpDF$Sample %in% DF_CT_current$Sample, ]
                DF_CT_current$Sample <- factor(DF_CT_current$Sample, levels = LookUpDF_current$Sample, ordered = TRUE)
                
                if (is.null(max_abundance_for_color)) {
                        max_abundance_for_color_current <- quantile(DF_CT_current$Count, .9)
                } else {
                        max_abundance_for_color_current <- max_abundance_for_color
                }
                
                if (max_abundance_for_color_current == 0) {max_abundance_for_color_current = max(DF_CT_current$Count)}
                
                # Color the sample names based on color_levels
                colxaxis <- color_levels[LookUpDF_current$Group]
                
                
                hmTr <- ggplot(DF_CT_current, aes(x = Sample, y = Taxa, fill = Count))
                hmTr <- hmTr + 
                        geom_raster() + 
                        scale_fill_gradientn("", limits = c(0, max_abundance_for_color_current), colors = c("red", viridis(5)), values = gradient_steps, oob = squish) +
                        scale_x_discrete(position = "top") +
                        #coord_equal() +
                        labs(x=NULL, y=NULL) +
                        theme_tufte(base_family = "Helvetica") +
                        theme(axis.ticks=element_blank(),
                              axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0,
                                                         colour = colxaxis))
                
                plot_list[[k]] <- hmTr
                names(plot_list)[k] <- paste(fac_levels[i], "_vs_", fac_levels[j], sep = "")
        }
        plot_list
}



#######################################
### plot_toptaxa_boxAndviolin##
#################
# generates box and violin plots of the taxons in physeq in the order given by tax_order and named by tax_names
# if tax_order and tax_names are NULL the order and names of physeq will be used
# facet_cols determines ncol in facet_wrap
# color_levels must be a named vector of colors for all levels found in group_var
# if ttestp = "yes" ggpubr is used to add significance asterisk from simple t.tests

# OUTPUT:
# generates for each level combination of the grouping factor (defined by group_var) seven plots, so for each combi a list of 8 plots,
# specifically: boxplot, boxplot faceted, violin plot, violin plot faceted, and the same again for y axis log10 (NB: need to make new plots
# there because that needs pseydocounts for looking good)

plot_toptaxa_boxAndviolin <- function(physeq, group_var, tax_order = NULL, tax_names = NULL, facet_cols = 5, color_levels, ttestp = "yes"){
        
        if(taxa_are_rows(physeq)) { physeq <- t(physeq) }
        # I prefer taxa_are_rows = FALSE so rows (= observations = samples), and colums = variables = taxa
        
        if (!is.factor(sample_data(physeq)[[group_var]])) {sample_data(physeq)[[group_var]] <- as.factor(sample_data(physeq)[[group_var]])}
        
        CT <- as(otu_table(physeq), "matrix") # taxa are columns
        DF_CT <- as.data.frame(t(CT))
        
        if (is.null(tax_order)){tax_order <- taxa_names(physeq)}
        
        if (length(tax_order) != nrow(DF_CT) || !all(tax_order %in% rownames(DF_CT))) {stop("the given tax_order must contain all taxa_names(physeq). If you want to subset,
                                                                                            do it before with prune_taxa")}
        DF_CT <- DF_CT[tax_order, ]
        
        if (is.null(tax_names)) {tax_names <- paste("Taxon_", 1:nrow(DF_CT), sep = "")}
        
        if (length(tax_names) != nrow(DF_CT) ) {stop("the given tax_names must be a character vector of length = ntaxa(physeq)")}
        
        rownames(DF_CT) <- make.unique(tax_names)
        DF_CT$Taxa <- factor(rownames(DF_CT), levels = rownames(DF_CT), ordered = TRUE)
        
        DF_CT <- tidyr::gather(DF_CT, key = Sample , value = Count, -Taxa)
        
        LookUpDF <- data.frame(Sample = sample_names(physeq), Group = sample_data(physeq)[[group_var]])
        DF_CT$Group <- LookUpDF$Group[match(DF_CT$Sample, LookUpDF$Sample)]
        
        group_fac <- factor(sample_data(physeq)[[group_var]])
        
        fac_levels <- levels(group_fac)
        
        # - get the level combis -
        fac_levels_num <- setNames(seq_along(fac_levels), fac_levels) 
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
        # --
        
        
        plot_list <- vector("list", length = length(i_s))
        
        for (k in seq_along(i_s)){
                i <- i_s[k]
                j <- j_s[k]
                
                group_fac_current <- droplevels(group_fac[as.numeric(group_fac) %in% c(i, j)])
                DF_CT_current <- filter(DF_CT, Group %in% group_fac_current)
                DF_CT_current$Group <- factor(DF_CT_current$Group, levels = c(fac_levels[i], fac_levels[j]), ordered = TRUE)
                
                Tr <- ggplot(DF_CT_current, aes(x = Taxa, y = Count, col = Group))
                Tr <- Tr +
                        geom_boxplot(outlier.color = NA) +
                        geom_point(size = 1, alpha = 0.6, position = position_jitterdodge()) +
                        scale_color_manual("", values = color_levels) +
                        xlab("") +
                        ylab("abundance") +
                        theme_bw() +
                        theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0),
                              legend.position = "none")
                
                if (ttestp == "yes"){
                        Tr <- Tr + ggpubr::stat_compare_means(label = "p.signif", method = "t.test") # t.test significance
                }
                
                
                Tr1 <- ggplot(DF_CT_current, aes(x = Group, y = Count, col = Group))
                Tr1 <- Tr1 +
                        geom_boxplot(outlier.color = NA) +
                        geom_point(size = 1, alpha = 0.6, position = position_jitterdodge()) +
                        facet_wrap(~ Taxa, ncol = facet_cols, scales = "free_y") +
                        scale_color_manual("", values = color_levels) +
                        xlab("") +
                        ylab("abundance") +
                        theme_bw() +
                        theme(legend.position = "none")
                
                if (ttestp == "yes"){
                        Tr1 <- Tr1 + ggpubr::stat_compare_means(label = "p.signif", method = "t.test", label.x = 1.5) # t.test significance
                }
                
                
                Tr2 <- ggplot(DF_CT_current, aes(x = Taxa, y = Count, col = Group))
                Tr2 <- Tr2 +
                        geom_violin(fill = NA) +
                        geom_point(size = 1, alpha = 0.6, position = position_jitterdodge()) +
                        scale_color_manual("", values = color_levels) +
                        xlab("") +
                        ylab("abundance") +
                        theme_bw() +
                        theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0),
                              legend.position = "none")
                
                if (ttestp == "yes"){
                        Tr2 <- Tr2 + ggpubr::stat_compare_means(label = "p.signif", method = "t.test") # t.test significance
                }
                
                Tr3 <- ggplot(DF_CT_current, aes(x = Group, y = Count, col = Group))
                Tr3 <- Tr3 +
                        geom_violin(fill = NA) +
                        geom_point(size = 1, alpha = 0.6, position = position_jitterdodge()) +
                        facet_wrap(~ Taxa, ncol = facet_cols, scales = "free_y") +
                        scale_color_manual("", values = color_levels) +
                        xlab("") +
                        ylab("abundance") +
                        theme_bw() +
                        theme(legend.position = "none")
                
                if (ttestp == "yes"){
                        Tr3 <- Tr3 + ggpubr::stat_compare_means(label = "p.signif", method = "t.test", label.x = 1.5) # t.test significance
                }
                
                
                
                # - add a pseudocount for log10 plots - # necessary for good boxplots
                DF_CT_current$Count[DF_CT_current$Count == 0] <- min(DF_CT_current$Count[DF_CT_current$Count > 0])
                # --
                
                Tr4 <- ggplot(DF_CT_current, aes(x = Taxa, y = Count, col = Group))
                Tr4 <- Tr4 +
                        geom_boxplot(outlier.color = NA) +
                        geom_point(size = 1, alpha = 0.6, position = position_jitterdodge()) +
                        scale_color_manual("", values = color_levels) +
                        xlab("") +
                        ylab("abundance") +
                        theme_bw() +
                        theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0),
                              legend.position = "none") +
                        scale_y_log10()
                
                if (ttestp == "yes"){
                        Tr4 <- Tr4 + ggpubr::stat_compare_means(label = "p.signif", method = "t.test") # t.test significance
                }
                
                
                Tr5 <- ggplot(DF_CT_current, aes(x = Group, y = Count, col = Group))
                Tr5 <- Tr5 +
                        geom_boxplot(outlier.color = NA) +
                        geom_point(size = 1, alpha = 0.6, position = position_jitterdodge()) +
                        facet_wrap(~ Taxa, ncol = facet_cols, scales = "free_y") +
                        scale_color_manual("", values = color_levels) +
                        xlab("") +
                        ylab("abundance") +
                        theme_bw() +
                        scale_y_log10() +
                        theme(legend.position = "none")
                
                if (ttestp == "yes"){
                        Tr5 <- Tr5 + ggpubr::stat_compare_means(label = "p.signif", method = "t.test", label.x = 1.5) # t.test significance
                }
                
                
                
                Tr6 <- ggplot(DF_CT_current, aes(x = Taxa, y = Count, col = Group))
                Tr6 <- Tr6 +
                        geom_violin(fill = NA) +
                        geom_point(size = 1, alpha = 0.6, position = position_jitterdodge()) +
                        scale_color_manual("", values = color_levels) +
                        xlab("") +
                        ylab("abundance") +
                        theme_bw() +
                        theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0),
                              legend.position = "none") +
                        scale_y_log10()
                
                if (ttestp == "yes"){
                        Tr6 <- Tr6 + ggpubr::stat_compare_means(label = "p.signif", method = "t.test") # t.test significance
                }
                
                Tr7 <- ggplot(DF_CT_current, aes(x = Group, y = Count, col = Group))
                Tr7 <- Tr7 +
                        geom_violin(fill = NA) +
                        geom_point(size = 1, alpha = 0.6, position = position_jitterdodge()) +
                        facet_wrap(~ Taxa, ncol = facet_cols, scales = "free_y") +
                        scale_color_manual("", values = color_levels) +
                        xlab("") +
                        ylab("abundance") +
                        theme_bw() +
                        scale_y_log10() +
                        theme(legend.position = "none")
                
                if (ttestp == "yes"){
                        Tr7 <- Tr7 + ggpubr::stat_compare_means(label = "p.signif", method = "t.test", label.x = 1.5) # t.test significance
                }
                
                plot_list[[k]] <- list(Tr, Tr1, Tr2, Tr3, Tr4, Tr5, Tr6, Tr7)
                names(plot_list)[k] <- paste(fac_levels[i], "_vs_", fac_levels[j], sep = "")
        }
        
        plot_list
}

#######################################
### FUNCTION: test_diffs_in_prevalence
#######################################
# Function performs fisher exact test on prevalence (absence presence) between the levels
# in a grouping factor
# INPUT:
# physeq: phyloseq
# group_var: name of the column in sample_data(physeq) that defines the groups
# p.adj.method, used in p.adjust
# minCount: present are taxa in species with more counts than minCount
# OUTPUT:
# list of data.frames, one data frame for each combi of levels in your grouping factor
# The data frames are ordered by p_value, and the tax_table has been cbound:)

test_diffs_in_prevalence <- function(physeq, group_var, p.adj.method = "fdr", minCount = 0L, symnum.args = list(cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))) {
        
        
        ### The Fisher test of pooled successes
        # MatrixFisher <- matrix(c(200, 80, 200, 320), nrow = 2)
        # rownames(MatrixFisher) <- c("Strain A", "Strain B")
        # colnames(MatrixFisher) <- c("Success", "Fail")
        # fisher.test(MatrixFisher, alternative = "two")
        # PValueFisher <- fisher.test(MatrixFisher, alternative = "two")$p.value
        
        if (taxa_are_rows(physeq)) { physeq <- t(physeq)}
        
        CT <- as(otu_table(physeq), "matrix")
        CT <- CT > minCount
        
        group_fac <- factor(sample_data(physeq)[[group_var]])
        
        fac_levels <- levels(group_fac)
        
        prev_list <- lapply(fac_levels, function(level){
                data.frame(Present = colSums(CT[group_fac == level, ]),
                         Absent = colSums(!CT[group_fac == level, ]))
        })
        
        names(prev_list) <- fac_levels
        
        # - get the level combis -
        fac_levels_num <- setNames(seq_along(fac_levels), fac_levels) 
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
        # --
        
        p_val_list <- vector("list", length = length(i_s))
        for (k in seq_along(i_s)){
                i <- i_s[k]
                j <- j_s[k]
                prev_PC_1 <- round(100*(prev_list[[i]][, 1]/(prev_list[[i]][1,1] + prev_list[[i]][1,2])), 1)
                prev_PC_2 <- round(100*(prev_list[[j]][, 1]/(prev_list[[j]][1,1] + prev_list[[j]][1,2])), 1)
                direction <- rep("down", length(prev_PC_1))
                direction[prev_PC_1 > prev_PC_2] <- "up"
                rowwise_compare_matrix <- cbind(prev_list[[i]], prev_list[[j]])
                p_vals <- sapply(1:nrow(rowwise_compare_matrix), function(e){
                        mat_fisher <- matrix(c(rowwise_compare_matrix[e, 1],
                                               rowwise_compare_matrix[e, 3],
                                               rowwise_compare_matrix[e, 2],
                                               rowwise_compare_matrix[e, 4]), ncol = 2)
                        fisher.test(mat_fisher)$p.value
                })
                p_vals_adj <- p.adjust(p_vals, p.adj.method)
                symnum.args$x <- p_vals
                significance <- do.call(stats::symnum, symnum.args) %>% as.character()
                symnum.args$x <- p_vals_adj
                significance_adj <- do.call(stats::symnum, symnum.args) %>% as.character()
                p_val_df <- data.frame(p_vals, p_vals_adj, significance, significance_adj, prev_PC_1, prev_PC_2, direction)
                # colnames(p_val_df)[1:2] <- paste(c("p_val_", "p_val_adj_"), fac_levels[i], "_vs_", fac_levels[j], sep = "")
                # colnames(p_val_df)[5:6] <- paste(c("prev_PC_"), c(fac_levels[i], fac_levels[j]), sep = "")
                colnames(p_val_df)[1:6] <- c("p_val", "p_val_adj", "signi.", "signi_adj", "prev_PC_grp1", "prev_PC_grp2")
                p_val_df <- cbind(as.data.frame(p_val_df), tax_table(physeq))
                p_val_df$Taxon <- colnames(CT)
                p_val_df <- arrange(p_val_df, p_val)
                p_val_df <- select(p_val_df, Taxon, 1:(ncol(p_val_df)-1))
                p_val_list[[k]] <- p_val_df
                names(p_val_list)[[k]] <- paste(fac_levels[i], "_vs_", fac_levels[j], sep = "")
        }
        
        p_val_list
        
}


#######################################
### DESeq2Apply_physeq
#######################################
## Inputs
# physeq: phyloseq object
# group_var: name of column that defines group fac in sample_data
# SFs: often you might want to give the SizeFactors already because you wanted to calculate them on non-filtered data,
# when SFs are not NULL, type is ignored
# type: type in estimateSizeFactors, ignored when Size factors given
## OUTPUT:
# list of two lists: the first: List of DESeq2 results plus of fisher.exact test, for each level combination in group factor one data_frame = list entry
# in the second list is just the size factor adjusted physeq object


DESeq2Apply_physeq <- function(physeq, group_var, SFs = NULL, type = "ratio", p.adjust.method = "fdr", symnum.args = list(cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))){
        
        if (taxa_are_rows(physeq)) {physeq <- t(physeq)}
        
        group_fac <- factor(sample_data(physeq)[[group_var]])
        fac_levels <- levels(group_fac)
        
        # NB: DESeq2 can not deal with ordered factors, sees them somehow as just one level, therefore
        sample_data(physeq)[[group_var]] <- factor(sample_data(physeq)[[group_var]], levels = levels(sample_data(physeq)[[group_var]]), ordered = FALSE)
        
        DES = phyloseq::phyloseq_to_deseq2(physeq, formula(paste("~", group_var)))
        
        
        if (is.null(SFs)){
                if(type == "ratio"){
                        GM <- apply(otu_table(physeq), 2, gm_own, zeros.count = FALSE)
                }
                
                dds <- estimateSizeFactors(DES, type = type, geoMeans = GM)
                # NB: geoMeans is ignored when type = "iterate"
                # NB2: "iterate" takes much longer than "ratio", and when using GM via gm_own with "ratio" the size factors 
                # correlate with above 99% with size factors from "iterate"
                # SFs2 <- sizeFactors(dds)
                
        } else {
                dds <- DES
                sizeFactors(dds) = SFs
                # identical(sizeFactors(dds), SFs) 
        }
        
        
        dds <- estimateDispersions(dds, quiet = TRUE) 
        dds <- nbinomWaldTest(dds)
        
        # to get the size factor adjusted physeq object
        physeq_out <- physeq
        otu_table(physeq_out) <- otu_table(t(counts(dds, normalized = TRUE)), taxa_are_rows = FALSE)
        
        # - get the level combis, NB: dds contains result infos on all combinations! -
        fac_levels_num <- setNames(seq_along(fac_levels), fac_levels) 
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
        # --
        
        result_list <- vector("list", length = length(i_s))
        
        for (k in seq_along(i_s)) {
                i <- i_s[k]
                j <- j_s[k]
                
                res <- as.data.frame(results(dds, contrast = c(group_var, fac_levels[i], fac_levels[j])))
                
                res$p_val_adj <- p.adjust(res$pvalue, method = p.adjust.method) # NB: in case of "fdr" same as default DESeq2
                
                CT <- counts(dds, normalized = TRUE)
                n1 <- sum(group_fac == fac_levels[i])
                n2 <- sum(group_fac == fac_levels[j])
                res$Median_grp1 <- apply(CT[, group_fac == fac_levels[i]], 1, median)
                res$Median_grp2 <- apply(CT[, group_fac == fac_levels[j]], 1, median)
                res$Mean_grp1 <- apply(CT[, group_fac == fac_levels[i]], 1, mean)
                res$Mean_grp2 <- apply(CT[, group_fac == fac_levels[j]], 1, mean)
                # res$baseMeanSelf <- apply(CT, 1, mean) # exactly the same as baseMean!
                res$Zeros_grp1 <- apply(CT[, group_fac == fac_levels[i]], 1, function(cnts){sum(cnts == 0)})
                res$Zeros_grp2 <- apply(CT[, group_fac == fac_levels[j]], 1, function(cnts){sum(cnts == 0)})
                res$Present_grp1 <- apply(CT[, group_fac == fac_levels[i]], 1, function(cnts){sum(cnts != 0)})
                res$Present_grp2 <- apply(CT[, group_fac == fac_levels[j]], 1, function(cnts){sum(cnts != 0)})
                res$prev_PC_grp1 <- round(100*(res$Present_grp1/n1),1)
                res$prev_PC_grp2 <- round(100*(res$Present_grp2/n2), 1)
                #res$Sparsity_grp1 <- 100*(res$Zeros_grp1/n1)
                #res$Sparsity_grp2 <- 100*(res$Zeros_grp2/n2)
                res$n1 <- n1
                res$n2 <- n2
                
                # - add fisher exact test of sparsity/prevalence again -
                Fisher <- t(sapply(1:nrow(res), FUN = function(i){
                        fisherMat <- matrix(c(res$Present_grp1[i], res$Zeros_grp1[i], res$Present_grp2[i],
                                              res$Zeros_grp2[i]), ncol = 2, dimnames = list(c("Present", "Absent"), c("grp1", "grp2")))
                        Test <- fisher.test(fisherMat)
                        cbind(Test$p.value, Test$estimate)
                }))
                
                p_val_adj_Fisher <- p.adjust(Fisher[,1], method = p.adjust.method)
                
                res$p_val_Fisher <- Fisher[,1]
                res$p_val_Fisher_adj <- p_val_adj_Fisher
                res$oddsRatioFisher <- Fisher[,2]
                # --
                
                # - add sginificance and direction -
                symnum.args$x <- res$pvalue
                res$signi. <- do.call(stats::symnum, symnum.args) %>% as.character()
                symnum.args$x <- res$p_val_adj
                res$signi_adj <- do.call(stats::symnum, symnum.args) %>% as.character()
                symnum.args$x <- res$p_val_Fisher
                res$signif_fisher <- do.call(stats::symnum, symnum.args) %>% as.character()
                symnum.args$x <- res$p_val_Fisher_adj
                res$signif_fisher_adj <- do.call(stats::symnum, symnum.args) %>% as.character()
                
                res$direction <- rep("down", nrow(res))
                res$direction[res$Median_grp1 > res$Median_grp2] <- "up"
                res$direction_fisher <- rep("down", nrow(res))
                res$direction_fisher[res$prev_PC_grp1 > res$prev_PC_grp2] <- "up"
                # --
                
                res$Taxon <- rownames(res)
                
                res <- dplyr::select(res, Taxon, teststat = stat, p_val = pvalue, p_val_adj,
                                     signi., signi_adj, direction, p_val_Fisher,
                                     p_val_Fisher_adj, signif_fisher, signif_fisher_adj,
                                     direction_fisher, Median_grp1, Median_grp2, Mean_grp1,
                                     Mean_grp2, Present_grp1, Present_grp2, Zeros_grp1, Zeros_grp2, prev_PC_grp1, prev_PC_grp2, 
                                     n1, n2, baseMean, log2FoldChange, lfcSE, oddsRatioFisher)
                # NB: I dropped here padj from DESeq since same as p_val_adj in case of p.adjust.method = "fdr"
                res <- cbind(res, tax_table(physeq))
                res <- dplyr::arrange(res, desc(abs(teststat)))
                # res <- res[order(res$p_val),]
                result_list[[k]] <- res
                names(result_list)[[k]] <- paste(fac_levels[i], "_vs_", fac_levels[j], sep = "")
        }
        
        out <- list(result_list, physeq_out)
        
}


#######################################
### wilcoxTestApply_physeq
#######################################
# Accepts several groups in group_var and results in a list of DF for the pairwise comparisons
# directly related to wilcoxTestApply from SimulationFunctions.R just not looping through a list of phyloseq objects
# but doing the job on a single phyloseq object
# does wilcoxon test and adds further info
# in addition it performs a fisher exact test on the sparsity proportions
# NB: Tested: 
# with excludeZeros you can decide on whether 0 counts should be excluded for the wilcox.test, the fisher sparsity test is
# of course not affected. 
# NB: in case in one of the two groups all counts are 0 and excludeZeros = T, then NA is given for all wilcoxon 
# statistics!
# The teststatistic is based on the standardized teststatistic, equation provided by multtest::mt.minP (compare with mtApply)
# (see equation for standStat2 in the code, results in exactly the same as standStat when uncomment)
## Input
# physeq object
# group_var: refers to a factor in sample_data(phyloseq) that defines the groups
# excludeZeros: decides on whether 0s should be considered when comparing the groups in a wilcox.test
# p.adjust.method, used as method in p.adjust
## Output
# list of dataframes with the results for each pairwise group comparison in group_var associated group_fac

wilcoxTestApply_physeq <- function(physeq, group_var, excludeZeros = FALSE, p.adjust.method = "fdr", symnum.args = list(cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))) {
        
        
        if (taxa_are_rows(physeq)) { physeq <- t(physeq)}
        
        CT <- as(otu_table(physeq), "matrix")
        
        group_fac <- factor(sample_data(physeq)[[group_var]])
        
        fac_levels <- levels(group_fac)
        
        # - get the level combis, NB: dds contains result infos on all combinations! -
        fac_levels_num <- setNames(seq_along(fac_levels), fac_levels) 
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
        # --
        
        result_list <- vector("list", length = length(i_s))
        
        for (k in seq_along(i_s)){
                i <- i_s[k]
                j <- j_s[k]
                
                res_mat <- apply(CT, 2, function(taxon_counts){
                        x <- taxon_counts[as.numeric(group_fac) == i]
                        Zeros_grp1 <- sum(x == 0)
                        # Sparsity_grp1 <- 100*(Zeros_grp1/length(x))
                        Present_grp1 <- length(x)-Zeros_grp1
                        prev_PC_grp1 <- 100*(Present_grp1/length(x))
                        if(excludeZeros){
                                x <- x[x != 0]
                        }
                        Median_grp1 <- median(x, na.rm = T) # NA in case all 0
                        Mean_grp1 <- mean(x, na.rm = T) # NaN in case all 0
                        if (is.na(Mean_grp1)){ Mean_grp1 = NA }
                        y <- taxon_counts[as.numeric(group_fac) == j]
                        Zeros_grp2 <- sum(y == 0)
                        Present_grp2 <- length(y)-Zeros_grp2
                        # Sparsity_grp2 <- 100*(Zeros_grp2/length(y))
                        prev_PC_grp2 <- 100*(Present_grp2/length(y))
                        if(excludeZeros){
                                #if(all(y == 0)){y[1] <- ceiling(mean(taxon_counts))+1}
                                y <- y[y != 0]
                        }
                        Median_grp2 <- median(y, na.rm = T)
                        Mean_grp2 <- mean(y, na.rm = T)
                        if (is.na(Mean_grp2)){ Mean_grp2 = NA }
                        
                        if (length(x) != 0 && length(y) != 0){
                                wilcTest <- wilcox.test(x = x, y = y, alternative = "two", paired = F, exact = F)
                                pValue <- wilcTest$p.value
                                W <- wilcTest$statistic
                                # calculate standardized rank sum Wilcoxon statistics as in multtest::mt.minP
                                Ranks <- rank(c(x, y))
                                n1 <- length(x)
                                n2 <- length(y)
                                # Wx <- sum(Ranks[1:n1])-(n1*(n1+1)/2) # would be the same as W
                                # how about the other W?
                                # Wy <- sum(Ranks[(n1+1):n2]) - (n2*(n2+1)/2)
                                standStat <- -1*((sum(Ranks[1:n1]) - n1*(n1+n2+1)/2)/sqrt(n1*n2*(n1+n2+1)/12))
                                
                                # # if you want to check that multtest::mt.minP would give the same statistic
                                # mati <- matrix(c(x,y), nrow = 1)
                                # grFac <- c(rep(fac_levels[i], n1), rep(fac_levels[j], n2))
                                # grFac <- factor(grFac, levels = c(fac_levels[i], fac_levels[j]))
                                # standStat2 <- multtest::mt.minP(mati, grFac, test = "wilcoxon")$teststat
                                # # identical(standStat, standStat2) # TRUE
                                # uncomment all with standStat2 to test all the way
                                
                        } else {
                                pValue = NA
                                W <- NA
                                standStat = NA
                                n1 <- length(x)
                                n2 <- length(y)
                                # standStat2 = NA
                        }
                        
                        
                        # -- add fisher exact test of presence differences (should be none in simulation) --
                        fisherMat <- matrix(c(Present_grp1, Zeros_grp1, Present_grp2, Zeros_grp2), 
                                            ncol = 2, dimnames = list(c("Present", "Absent"), c("grp1", "grp2")) )
                        Test <- fisher.test(fisherMat)
                        
                        c(teststat = standStat, p_val = pValue, Median_grp1 = Median_grp1, Median_grp2 = Median_grp2, 
                          Mean_grp1 = Mean_grp1, Mean_grp2 = Mean_grp2, n1 = n1, n2 = n2, Present_grp1 = Present_grp1, 
                          Present_grp2 = Present_grp2, Zeros_grp1 = Zeros_grp1, Zeros_grp2 = Zeros_grp2, 
                          prev_PC_grp1 = prev_PC_grp1, prev_PC_grp2 = prev_PC_grp2, W, 
                          p_val_Fisher = Test$p.value, Test$estimate) #, teststat2 = standStat2 # waste of typing effort:)
                })
                
                res_mat <- t(res_mat)
                colnames(res_mat) <- c("teststat", "p_val", "Median_grp1", "Median_grp2", "Mean_grp1", "Mean_grp2", "n1",
                                       "n2", "Present_grp1", "Present_grp2", "Zeros_grp1", "Zeros_grp2", "prev_PC_grp1", "prev_PC_grp2", "W",
                                       "p_val_Fisher", "oddsRatioFisher") #, , "teststat2"
                
                
                DF <- data.frame(Taxon = rownames(res_mat), res_mat)
                DF$p_val_adj <- p.adjust(DF$p_val, method = p.adjust.method)
                DF$p_val_Fisher_adj <- p.adjust(DF$p_val_Fisher, method = p.adjust.method)
                
                # - add sginificance and direction -
                symnum.args$x <- DF$p_val
                DF$signi. <- do.call(stats::symnum, symnum.args) %>% as.character()
                symnum.args$x <- DF$p_val_adj
                DF$signi_adj <- do.call(stats::symnum, symnum.args) %>% as.character()
                symnum.args$x <- DF$p_val_Fisher
                DF$signif_fisher <- do.call(stats::symnum, symnum.args) %>% as.character()
                symnum.args$x <- DF$p_val_Fisher_adj
                DF$signif_fisher_adj <- do.call(stats::symnum, symnum.args) %>% as.character()
                
                DF$direction <- rep("down", nrow(DF))
                DF$direction[DF$Median_grp1 > DF$Median_grp2] <- "up"
                DF$direction_fisher <- rep("down", nrow(DF))
                DF$direction_fisher[DF$prev_PC_grp1 > DF$prev_PC_grp2] <- "up"
                # --
                
                
                DF <- dplyr::select(DF, 1:3, p_val_adj, signi., signi_adj, direction, p_val_Fisher, p_val_Fisher_adj, signif_fisher,
                                    signif_fisher_adj, direction_fisher, 4:7, Present_grp1, Present_grp2, prev_PC_grp1, prev_PC_grp2)
                # teststat/standStat2 version:
                # DF <- dplyr::select(DF, 1:2, 19, 3, 20:22, 17, 23:25, 4:7, 10:15, 8:9, 16, 18)
                
                DF <- cbind(DF, tax_table(physeq))
                # DF <- dplyr::arrange(DF, desc(abs(teststat)))
                DF <- dplyr::arrange(DF, p_val)
                
                result_list[[k]] <- DF
                names(result_list)[[k]] <- paste(fac_levels[i], "_vs_", fac_levels[j], sep = "")
        }
        
        result_list
}


#######################################
### prepare_diff_abundance_results_for_plotting
#######################################
# function exclusively to save some typing, and for forcing me to make the outcomes of the diff_abundance test functions more uniform
# TbT is just because the TbT method has currently no p_val_adj yet! Has to be changed

prepare_diff_abundance_results_for_plotting <- function(res_list, physeq, taxonomic_level, TbT = "no") {
        
        if(TbT == "yes"){
                head_values <- rep(10, times = length(res_list))
                names(head_values) <- names(res_list)
        } else {
                suppressWarnings(head_values <- sapply(res_list, function(df){
                        max(which(df[, "p_val_adj"] < 0.05))
                }))
        }
        
        original_head_values <- data.frame(Comparison = c(names(head_values), "Total"), NoSignificant = c(head_values, ntaxa(physeq)),
                                           PC_Significant = 100*c(head_values, ntaxa(physeq))/ntaxa(physeq))
        
        # Show at least 10 taxa even if less are significant and show max 25 even if far more are significant
        head_values[head_values < 10] <- 10
        head_values[head_values > 25] <- 25
        
        head_values[head_values > ntaxa(physeq)] <- ntaxa(physeq) # for "Phylum" analysis it is possible that less than 10 taxa are in physeq
        
        res_table_list <- lapply(1:length(res_list), function(i){
                df <- head(res_list[[i]], head_values[i])
                df <- dplyr::select(df, which(colnames(df) == taxonomic_level), 2:ncol(df)) #simply put taxonomic level in front and remove taxon column (= column 1) 
                df
        })
        names(res_table_list) <- names(res_list)
        
        
        row_names_for_heat_maps <- lapply(res_table_list, function(df){
                the_names <- as.character(df[, taxonomic_level])
                the_names <- sapply(strsplit(the_names, split = "/"), `[`, 1)
                the_names <- paste(the_names, df$signi., df$signi_adj, sep = "_")
        }) # in case of ambiguous species assignment keep only first one
        
        tax_orders <- lapply(1:length(res_list), function(i){
                df <- res_list[[i]]
                as.character(df$Taxon[1:head_values[i]])
        })
        
        pruned_physeqs_to_test <- lapply(1:length(tax_orders), function(i){
                prune_taxa(tax_orders[[i]], physeq)
        })
        
        list(original_head_values = original_head_values, head_values = head_values, res_table_list = res_table_list, row_names_for_heat_maps = row_names_for_heat_maps,
             tax_orders = tax_orders, pruned_physeqs_to_test = pruned_physeqs_to_test)
        
}



####################################
## calculate_TbTmatrixes:
###################################

# The method is based on direct comparisons between Taxa (one to one ratio comparisons), thus avoiding compositionality effects.
# For each "host taxon" a matrix is created, where the host taxons counts are directly compared to the
# counts of each other taxon.
# Specifically, in steps, example for count table of 100 taxons and 20 samples
# Step 1: calculate taxa by taxa ratios (one matrix per host taxon) to get rid of compositionality. 
# The ratios show for example count taxon 1/ count taxon 2 over all samples, illustrating the non-compositionality affected
# ratios of taxon 1 to taxon 2 in the different samples
# Step 2: divide the ratios by geometric mean over all samples, to make ratios of taxon1/taxon2 comparable to taxon1/taxon3 ratios and so on
# Step 3: log the ratios, so the sum of ratios/values over all samples = 0. 
# (see also implementation steps below to understand the values in the final TbT matrixes even better)
# NB: Zero Counts are fully ingored here: Reasoning on this treatment of zero counts: 
# Values/Ratios where one of the two taxa is missing (count = 0) should not be used to determine
# differential abundance (instead sparsity should be tested separately)
# so if the host taxon is 0 in a sample, then the entire sample will basically be ignored for that taxon.
# Values where the host taxon or the other taxon is 0 should be set to NA, so that all 0 values come
# from cases log(x/y) where x = y but both x and y != 0
# NB: further explanations within the code

## Input: 
# - physeq = a phyloseq object
# - group_var, name of the group_fac in sample_data(physeq)
## Output: 
# - list of TbTmatrixes, one list for each combi of levels in group_fac, list items are named by level_vs_level

# Calculation/implementation Steps (assume a count table of 100 taxons and 20 samples):
# Step 1: For each host taxon (e.g. taxon 1) calculate the within-sample taxon 1/taxon ratios and log these ratios,
# i.e. log(count taxon 1/count taxon 1:100) which is simplyd log(count taxon 1) - log(count taxon 1:100).
# Step 2: ratios/values where on of the counts were 0 are set to NA, and the log(geometric mean) over all samples
# is subtracted. Done

calculate_TbTmatrixes = function(physeq, group_var){
        
        if (taxa_are_rows(physeq)) {physeq <- t(physeq)}
        
        
        group_fac <- factor(sample_data(physeq)[[group_var]])
        
        fac_levels <- levels(group_fac)
        
        # - get the level combis, NB: dds contains result infos on all combinations! -
        fac_levels_num <- setNames(seq_along(fac_levels), fac_levels) 
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
        # --
        
        TbTmatrixes_list <- vector("list", length = length(i_s))
        
        for (k in seq_along(i_s)) {
                
                i <- i_s[k]
                j <- j_s[k]
                
                
                CT <- t(as(otu_table(physeq), 'matrix')) # now taxa are rows and samples are columns
                group_fac_current <- droplevels(group_fac[as.numeric(group_fac) %in% c(i, j)])
                CT <- CT[, group_fac %in% group_fac_current]
                
                
                # == Take in-sample taxa by taxa ratios and take the log ==
                # NB: log(x/y) = log(x) - log(y)
                TbTmatrixes <- lapply(1:nrow(CT), function(i){apply(CT, 2, function(samp_cnts){log(samp_cnts[i]) - log(samp_cnts)})})
                # produces for each taxon (= host taxon) a TbTMatrix
                # NB: there are -Inf, Inf, and NaN values in the matrixes, specifically
                # 0/x = log(0) - log(x) = -Inf, x/0 = log(x) - log(0) = Inf; 0/0 = log(0) - log(0) = NaN!
                # Inf, -Inf, and NaN will be ignored for calculating the rowMeans (geometric means) in the
                # next step
                # ====
                
                # == dividing each row of a TbTMatrix by its geometric mean, and taking log ==
                # NB: because the ratios are already "logged" it is just subtracting the rowMeans, see NB2
                # NB2: remember (see keep Note: #statistics #work geometric mean): geometric_mean(x) with x = x1, ... xn = (x1*...*xn)^1/n = exp(mean(log(x))), i.e. exp((1/n)*(log(x1)+ ... + log(xn))). 
                # Therefore, log(geometric_mean(x)) = (1/n)*(log(x1)+ ... + log(xn)) (= rowMeans when row = log(x1), ... log(xn))
                # Therefore: the rowMeans of the "logged" ratios are log(geometric_mean), and thus:
                # log(Ratio/geometric_mean) = log(Ratio) - rowMean! Remember log(Ratio) is currently in the TbTMatrixes
                # NB3: all Inf, -Inf, and NaN are set to NA, thus all 0 in the matrixes are from ratios x/y where x = y and x AND y != 0.
                TbTmatrixes <- lapply(TbTmatrixes, function(mat) {
                        mat[!is.finite(mat)] <- NA # puts all Inf, -Inf, and also NaN to NA!
                        mat-rowMeans(mat, na.rm = TRUE)
                        #newM[is.na(newM)] <- 0
                })
                # ====
                names(TbTmatrixes) <- rownames(TbTmatrixes[[1]])
                TbTmatrixes_list[[k]] <- TbTmatrixes
                names(TbTmatrixes_list)[[k]] <- paste(fac_levels[i], "_vs_", fac_levels[j], sep = "")
        }
        
        TbTmatrixes_list
        
}



####################################
## evaluate_TbTmatrixes: 
###################################
# NB: you need to give same physeq as in calculate_TbTmatrixes, of course two functions could be combined, would
# save some double calculations
# the function calculates the TbT_groupSums as "teststatistic" in the TbTmatrixes 
# it adds count infos from physeq, NB: Median_grp1 and so on are calculated only from non-zero
# values because TbT method ignores zeros
## Input: 
# - TbTmatrixes_list: The list with the list of TbTmatrixes for each level combination in group_var defined group factor
# - physeq: same physeq object as used in calculate_TbTmatrixes
# - group_var: defining the group factor in physeq, same again as in calculate_TbTmatrixes
# - p.adjust.method: here only for p.adjust on fisher p_values
## Output: 
# - list of result DF for each combination of levels in group_var defined group fac

evaluate_TbTmatrixes <- function(TbTmatrixes_list, physeq, group_var, p.adjust.method = "fdr", symnum.args = list(cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))) {
        
        if(!identical(length(TbTmatrixes_list[[1]]), ntaxa(physeq))){stop("TbTmatrixes can not fit to physeq")}
        
        CT <- t(as(otu_table(physeq), "matrix"))
        CT[CT == 0] <- NA # because the method ignores 0 counts, I want to ignore them when calculating Mean_grp1 and so on
        
        
        group_fac <- factor(sample_data(physeq)[[group_var]])
        
        fac_levels <- levels(group_fac)
        
        # - get the level combis, NB: dds contains result infos on all combinations! -
        fac_levels_num <- setNames(seq_along(fac_levels), fac_levels) 
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
        # --
        
        result_list <- vector("list", length = length(i_s))
        
        for (k in seq_along(i_s)) {
                
                i <- i_s[k]
                j <- j_s[k]
                
                TbTmatrixes <- TbTmatrixes_list[[k]]
                group_fac_current <- droplevels(group_fac[as.numeric(group_fac) %in% c(i, j)])
                
                sumGrp1 <- sapply(TbTmatrixes, function(mat){sum(mat[, group_fac_current == fac_levels[i]], na.rm = T)})
                sumGrp2 <- sapply(TbTmatrixes, function(mat){sum(mat[, group_fac_current == fac_levels[j]], na.rm = T)})
                sumAll <- sapply(TbTmatrixes, function(mat){sum(mat, na.rm = T)})
                DF <- data.frame(Taxon = names(sumAll), sumGrp1 = sumGrp1, sumGrp2 = sumGrp2, sumAll = sumAll, TbT_groupSum = pmax(sumGrp1, sumGrp2))
                
                CT_current <- CT[, group_fac %in% group_fac_current]
                
                DF$Median_grp1 <- apply(CT_current[, group_fac_current == fac_levels[i]], 1, median, na.rm = T)
                DF$Median_grp2 <- apply(CT_current[, group_fac_current == fac_levels[j]], 1, median, na.rm = T)
                DF$Mean_grp1 <- apply(CT_current[, group_fac_current == fac_levels[i]], 1, mean, na.rm = T)
                DF$Mean_grp2 <- apply(CT_current[, group_fac_current == fac_levels[j]], 1, mean, na.rm = T)
                DF$Present_grp1 <- apply(CT_current[, group_fac_current == fac_levels[i]], 1, function(cnts){sum(!is.na(cnts))})
                DF$Present_grp2 <- apply(CT_current[, group_fac_current == fac_levels[j]], 1, function(cnts){sum(!is.na(cnts))})
                DF$Zeros_grp1 <- apply(CT_current[, group_fac_current == fac_levels[i]], 1, function(cnts){sum(is.na(cnts))})
                DF$Zeros_grp2 <- apply(CT_current[, group_fac_current == fac_levels[j]], 1, function(cnts){sum(is.na(cnts))})
                n1 <- sum(group_fac_current == fac_levels[i])
                n2 <- sum(group_fac_current == fac_levels[j])
                # DF$Sparsity_grp1 <- 100*(DF$Zeros_grp1/n1)
                # DF$Sparsity_grp2 <- 100*(DF$Zeros_grp2/n2)
                DF$prev_PC_grp1 <- 100*(DF$Present_grp1/n1)
                DF$prev_PC_grp2 <- 100*(DF$Present_grp2/n2)
                DF$n1 <- n1
                DF$n2 <- n2
                # -- add fisher exact test of presence differences (should be none in simulation) --
                Fisher <- t(sapply(1:nrow(DF), FUN = function(i){
                        fisherMat <- matrix(c(DF$Present_grp1[i], DF$Zeros_grp1[i], DF$Present_grp2[i],
                                              DF$Zeros_grp2[i]), ncol = 2, dimnames = list(c("Present", "Absent"), c("grp1", "grp2")))
                        Test <- fisher.test(fisherMat)
                        cbind(Test$p.value, Test$estimate)
                }))
                
                DF$p_val_Fisher <- Fisher[,1]
                DF$p_val_Fisher_adj <- p.adjust(DF$p_val_Fisher, method = p.adjust.method)
                DF$oddsRatioFisher <- Fisher[,2]
                
                # - add sginificance and direction -
                symnum.args$x <- DF$p_val_Fisher
                DF$signif_fisher <- do.call(stats::symnum, symnum.args) %>% as.character()
                symnum.args$x <- DF$p_val_Fisher_adj
                DF$signif_fisher_adj <- do.call(stats::symnum, symnum.args) %>% as.character()
                
                DF$direction <- rep("down", nrow(DF))
                DF$direction[DF$Median_grp1 > DF$Median_grp2] <- "up"
                DF$direction_fisher <- rep("down", nrow(DF))
                DF$direction_fisher[DF$prev_PC_grp1 > DF$prev_PC_grp2] <- "up"
                # --
                
                DF <- dplyr::select(DF, Taxon, TbT_groupSum, direction, p_val_Fisher,
                                    p_val_Fisher_adj, signif_fisher, signif_fisher_adj,
                                    direction_fisher, Median_grp1, Median_grp2, Mean_grp1,
                                    Mean_grp2, Present_grp1, Present_grp2, Zeros_grp1, Zeros_grp2, prev_PC_grp1, prev_PC_grp2, 
                                    n1, n2, sumGrp1, sumGrp2, sumAll, oddsRatioFisher)
                
                
                DF <- cbind(DF, tax_table(physeq))
                DF <- dplyr::arrange(DF, desc(abs(TbT_groupSum)))
                result_list[[k]] <- DF
                names(result_list)[[k]] <- paste(fac_levels[i], "_vs_", fac_levels[j], sep = "")
        }
        
        result_list
}


####################################
## calculate_raw_TbTmatrixes:
###################################
# see calculate_TbTmatrixes, this one just calculates the "raw" ratio matrixes i.e. without log and gm

## Input: 
# - physeq = a phyloseq object
# - group_var, name of the group_fac in sample_data(physeq)
## Output: 
# - list of TbTmatrixes, one list for each combi of levels in group_fac, list items are named by level_vs_level


calculate_raw_TbTmatrixes = function(physeq, group_var){
        
        if (taxa_are_rows(physeq)) {physeq <- t(physeq)}
        
        
        group_fac <- factor(sample_data(physeq)[[group_var]])
        
        fac_levels <- levels(group_fac)
        
        # - get the level combis, NB: dds contains result infos on all combinations! -
        fac_levels_num <- setNames(seq_along(fac_levels), fac_levels) 
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
        # --
        
        TbTmatrixes_list <- vector("list", length = length(i_s))
        
        for (k in seq_along(i_s)) {
                
                i <- i_s[k]
                j <- j_s[k]
                
                
                CT <- t(as(otu_table(physeq), 'matrix')) # now taxa are rows and samples are columns
                group_fac_current <- droplevels(group_fac[as.numeric(group_fac) %in% c(i, j)])
                CT <- CT[, group_fac %in% group_fac_current]
                
                
                
                TbTmatrixes <- lapply(1:nrow(CT), function(i){apply(CT, 2, function(samp_cnts){samp_cnts[i]/samp_cnts})})
                # produces for each taxon (= host taxon) a TbTMatrix
                # NB: there are Inf, and NaN values in the matrixes, specifically
                # 0/x = 0, x/0 = Inf; 0/0 = NaN!
                
                names(TbTmatrixes) <- rownames(TbTmatrixes[[1]])
                TbTmatrixes_list[[k]] <- TbTmatrixes
                names(TbTmatrixes_list)[[k]] <- paste(fac_levels[i], "_vs_", fac_levels[j], sep = "")
        }
        
        TbTmatrixes_list
        
}



####################################
## plot_taxa_ratios_levelPairs 
###################################
## Input: 
# - TbTmatrixes_list: The list with the lists of TbTmatrixes for each level combi in group factor
# NB: here it should be raw_TbTmatrixes with 0 = 0/x, Inf = x/0, NaN = 0/0, all these values will be ignored in the ratio plots by being set to NA!!
# NB: NAs are also ignored in t.test or wilcoxon test
# - physeq: used for TbTmatrixes_list generation
# - group_var: the factor used to separate the samples
# - tax_names: the names of the taxa in physeq you want to use, e.g. taxa_names(physeq) if you like them, if NULL T1 to Tn will be used
# - taxa_nom: the nominator taxon of the abundance ratios, NB: only 1 allowed, must be included in tax_names
# - taxa_den: the denominator taxa of the abundance ratios, several allowed, all must be included in tax_names, if NULL all are used,
# i.e you get plots facet_wrapped around the taxa_den taxa
# NB: only taxa are plotted for which statistical test is possible!
# - test: either "t.test" or "wilcoxon"
# - p_adjust_method

## Output: 
# - list of pVals data frame (the pValues from t.test or wilcox.test after p.adjust) plus Violin plot plus Boxplot,
# so for each level combi in group_var a list of length 2.


plot_taxa_ratios_levelPairs <- function(TbTmatrixes_list, physeq, group_var, tax_names = NULL,
                                           taxa_nom = "Firmicutes", taxa_den = NULL, color_levels, test = "t.test", p_adjust_method = "fdr",
                                           symnum.args = list(cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))) {
        
        if(!identical(length(TbTmatrixes_list[[1]]), ntaxa(physeq))){stop("TbTmatrixes can not fit to physeq")}
        
        group_fac <- factor(sample_data(physeq)[[group_var]])
        
        fac_levels <- levels(group_fac)
        
        # - get the level combis -
        fac_levels_num <- setNames(seq_along(fac_levels), fac_levels) 
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
        # --
        
        result_list <- vector("list", length = length(i_s))
        
        LookUpDF <- data.frame(Sample = sample_names(physeq), Group = sample_data(physeq)[[group_var]])
        # order so the samples of a group are shown beside of each other (not really necessary for a violin plot)
        LookUpDF <- LookUpDF[order(match(LookUpDF$Group, levels(LookUpDF$Group))), ]
        
        
        for (k in seq_along(i_s)) {
                
                i <- i_s[k]
                j <- j_s[k]
                
                TbTmatrixes <- TbTmatrixes_list[[k]]
                
                if (is.null(tax_names)){
                        tax_names <- paste("T", 1:length(TbTmatrixes), sep = "_")
                } else {
                        if(!identical(length(TbTmatrixes), length(tax_names))){stop("tax_names do not fit in length to TbTmatrixes")}
                }
                
                names(TbTmatrixes) <- make.unique(tax_names)
                names(TbTmatrixes) <- tax_names #tax_names should be unique
                
                TbTmatrix <- TbTmatrixes[[taxa_nom]] # the only relevant matrix for this plot, it is where tax_nom counts have been divided by the counts of all other taxa
                
                if (is.null(TbTmatrix)) {stop("taxa_nom not found")}
                
                rownames(TbTmatrix) <- tax_names
                
                TbT_DF <- as.data.frame(TbTmatrix)
                TbT_DF$Taxon <- rownames(TbT_DF)
                # - use taxa_den (denominator) to restrict the taxa to which taxa_nom is compared to -
                if (is.null(taxa_den)) {taxa_den <- tax_names}
                TbT_DF <- TbT_DF[TbT_DF$Taxon %in% taxa_den, ]
                # --
                
                # - change to long DF -
                TbT_DF_l <- gather(TbT_DF, key = Sample, value = Ratio, -Taxon)
                # --
                # - add the group_var level information using the original level order in group_var -
                TbT_DF_l$Group <- as.character(LookUpDF$Group[match(TbT_DF_l$Sample, LookUpDF$Sample)])
                group_fac_current <- droplevels(group_fac[as.numeric(group_fac) %in% c(i, j)]) # keeps the current levels in original order
                TbT_DF_l$Group <- factor(TbT_DF_l$Group, levels = levels(group_fac_current), ordered = T)
                # --
                
                # - change all ratios where either nominator taxon or denominator taxon had count = 0 to NA -
                # remember: 0/0 = NaN (not finite), 0/x = 0, x/0 = Inf
                TbT_DF_l$Ratio[!is.finite(TbT_DF_l$Ratio) | TbT_DF_l$Ratio == 0] <- NA
                # --
                
                # - find taxa that would throw an error in statistical test and remove those taxa from DF -
                # first find the taxa that would throw an error in t.test or wilcox.test
                var_plus_length_check <- group_by(TbT_DF_l, Taxon, Group) %>% summarise(Variance = var(Ratio, na.rm = T), NotNA = sum(!is.na(Ratio)))
                if (test == "t.test"){
                        var_plus_length_check <- filter(var_plus_length_check, !(Variance > 0) | NotNA < 2)
                } else if (test == "wilcox.test") {
                        var_plus_length_check <- filter(var_plus_length_check, NotNA < 1)
                }
                
                if (nrow(var_plus_length_check) != 0) { # should always remove at least the case where taxa_denominator = taxa_nom, since those ratios are all 1, hence no variance
                        TbT_DF_l <- filter(TbT_DF_l, !(Taxon %in% unique(var_plus_length_check$Taxon)))
                }
                # --
                
                # - use ggpubr::compare_menas to calculate all pValues of the Ratios for the different taxa_den between current group_levels -
                pVals <- ggpubr::compare_means(formula = Ratio ~ Group, data = TbT_DF_l, group.by = "Taxon", method = test, p.adjust.method = p_adjust_method, symnum.args = symnum.args)
                
                pVals <- dplyr::arrange(pVals, p)
                
                # NB: the pVals change when you log the ratios (scale_y_log10())
                TbT_DF_l$RatioLog10 <- log10(TbT_DF_l$Ratio)
                pValsLog <- ggpubr::compare_means(formula = RatioLog10 ~ Group, data = TbT_DF_l, group.by = "Taxon", method = test, p.adjust.method = p_adjust_method, symnum.args = symnum.args)
                pValsLog <- dplyr::arrange(pValsLog, p)
                # --
                
                # NB: I plot now for log and non log independently, even though ggpubr is so smart to change the p-values when you just use
                # Tr + scale_y_log10(). I plot independently because log might change the order!
                TbT_DF_l_log <- TbT_DF_l
                
                # - order taxa based on pVals result -
                TbT_DF_l$Taxon <- factor(TbT_DF_l$Taxon, levels = pVals$Taxon, ordered = TRUE)
                TbT_DF_l_log$Taxon <- factor(TbT_DF_l_log$Taxon, levels = pValsLog$Taxon, ordered = TRUE)
                # --
                
                
                # - plot: NB: also for non removed taxa some samples might have NA ratios that will be removed -
                Tr <- ggplot(TbT_DF_l, aes(x = Group, y = Ratio, col = Group))
                Tr <- Tr +
                        geom_violin() +
                        geom_point(size = 1, alpha = 0.6, position = position_jitterdodge(dodge.width = 1)) +
                        # scale_color_manual(values = c(color_lookup$color[i], color_lookup$color[j])) +
                        scale_color_manual(values = color_levels) +
                        facet_wrap(~ Taxon, scales = "free_y") +
                        xlab("") +
                        ylab(paste("abundance ratio of", taxa_nom, "to stated taxon")) +
                        theme_bw() +
                        theme(legend.position = "none")
                
                
                Tr <- Tr + ggpubr::stat_compare_means(label = "p.signif", method = test, label.x = 1.5)
                
                
                Tr1 <- ggplot(TbT_DF_l, aes(x = Group, y = Ratio, col = Group))
                Tr1 <- Tr1 +
                        geom_boxplot(outlier.color = NA) +
                        geom_point(size = 1, alpha = 0.6, position = position_jitterdodge(dodge.width = 1)) +
                        # scale_color_manual(values = c(color_lookup$color[i], color_lookup$color[j])) +
                        scale_color_manual(values = color_levels) +
                        facet_wrap(~ Taxon, scales = "free_y") +
                        xlab("") +
                        ylab(paste("abundance ratio of", taxa_nom, "to stated taxon")) +
                        theme_bw() +
                        theme(legend.position = "none")
                
                
                Tr1 <- Tr1 + ggpubr::stat_compare_means(label = "p.signif", method = test, label.x = 1.5) # p.format
                
                
                Tr2 <- ggplot(TbT_DF_l_log, aes(x = Group, y = Ratio, col = Group))
                Tr2 <- Tr2 +
                        geom_violin() +
                        geom_point(size = 1, alpha = 0.6, position = position_jitterdodge(dodge.width = 1)) +
                        # scale_color_manual(values = c(color_lookup$color[i], color_lookup$color[j])) +
                        scale_color_manual(values = color_levels) +
                        facet_wrap(~ Taxon, scales = "free_y") +
                        xlab("") +
                        ylab(paste("abundance ratio of", taxa_nom, "to stated taxon")) +
                        theme_bw() +
                        theme(legend.position = "none")
                
                
                Tr2 <- Tr2 + ggpubr::stat_compare_means(label = "p.signif", method = test, label.x = 1.5)
                Tr2 <- Tr2 + scale_y_log10()
                
                
                Tr3 <- ggplot(TbT_DF_l_log, aes(x = Group, y = Ratio, col = Group))
                Tr3 <- Tr3 +
                        geom_boxplot(outlier.color = NA) +
                        geom_point(size = 1, alpha = 0.6, position = position_jitterdodge(dodge.width = 1)) +
                        # scale_color_manual(values = c(color_lookup$color[i], color_lookup$color[j])) +
                        scale_color_manual(values = color_levels) +
                        facet_wrap(~ Taxon, scales = "free_y") +
                        xlab("") +
                        ylab(paste("abundance ratio of", taxa_nom, "to stated taxon")) +
                        theme_bw() +
                        theme(legend.position = "none")
                
                
                Tr3 <- Tr3 + ggpubr::stat_compare_means(label = "p.signif", method = test, label.x = 1.5)
                Tr3 <- Tr3 + scale_y_log10()
                
                
                # --
                
                result_list[[k]] <- list(pVals = pVals, Tr = Tr, Tr1 = Tr1, pValsLog = pValsLog, Tr2 = Tr2, Tr3 = Tr3)
                # rm(Tr, TbT_DF_l, TbT_DF, TbT_DF_l_test, pVals)
        }
        
        names(result_list) <- names(TbTmatrixes_list)
        result_list
}


####################################
## plot_taxa_ratios_AllLevels 
###################################
# see plot_taxa_ratios_levelPairs: Here you directly calculate the count by count ratio matrix only for the tax_nom, you still facet by taxa_den
# (denominator) but you keep all levels in all plots. Makes extensive use of ggpubr, NB: ggpubr is so smart to adjust p-values when you use
# scale_y_log10

## Output: 
# - list(pVals = pVals, Tr = Tr, Tr1 = Tr1, pValsLog = pValsLog, Tr2 = Tr2, Tr3 = Tr3)


plot_taxa_ratios_AllLevels <- function(physeq, group_var, tax_names = NULL,
                                        taxa_nom = "Firmicutes", taxa_den = NULL, color_levels, test = "t.test", p_adjust_method = "fdr",
                                        symnum.args = list(cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))) {
        
        
        # - check that given tax_names fit to physeq and change taxa_names of physeq -
        if (is.null(tax_names)){
                tax_names <- paste("T", 1:length(TbTmatrixes), sep = "_")
        } 
        
        if(!identical(ntaxa(physeq), length(tax_names))){stop("tax_names do not fit in length to physeq")}
        
        tax_names <- make.unique(tax_names)
        taxa_names(physeq) <- tax_names
        # --
        
        
        if (taxa_are_rows(physeq)) {physeq <- t(physeq)}
        
        # - calculate the matrix taxa_nom/(all other taxa) -         
        CT <- t(as(otu_table(physeq), 'matrix')) # now taxa are rows and samples are columns
        
        i <- which(rownames(CT) == taxa_nom)
        if (length(i) != 1) {stop("taxa_nom not found in tax_names or tax_names not unique!")}
        
        
        TbTmatrix <- apply(CT, 2, function(samp_cnts){samp_cnts[i]/samp_cnts})
        # produces for each taxon (= host taxon) a TbTMatrix
        # NB: there are possibly Inf, and NaN values in the matrix, specifically
        # 0/x = 0, x/0 = Inf; 0/0 = NaN!
        # --
        
        TbT_DF <- as.data.frame(TbTmatrix)
        TbT_DF$Taxon <- rownames(TbT_DF)
        # - use taxa_den (denominator) to restrict the taxa to which taxa_nom is compared to -
        if (is.null(taxa_den)) {taxa_den <- tax_names}
        TbT_DF <- TbT_DF[TbT_DF$Taxon %in% taxa_den, ]
        # --
        
        # - change to long DF -
        TbT_DF_l <- gather(TbT_DF, key = Sample, value = Ratio, -Taxon)
        # --
        
        # - add the group_var level information  -
        LookUpDF <- data.frame(Sample = sample_names(physeq), Group = sample_data(physeq)[[group_var]])
        TbT_DF_l$Group <- as.character(LookUpDF$Group[match(TbT_DF_l$Sample, LookUpDF$Sample)])
        group_fac <- factor(sample_data(physeq)[[group_var]])
        TbT_DF_l$Group <- factor(TbT_DF_l$Group, levels = levels(group_fac), ordered = T)
        # --
        
        # - change all ratios where either nominator taxon or denominator taxon had count = 0 to NA -
        # remember: 0/0 = NaN (not finite), 0/x = 0, x/0 = Inf
        TbT_DF_l$Ratio[!is.finite(TbT_DF_l$Ratio) | TbT_DF_l$Ratio == 0] <- NA
        # --
        
        # - find taxa that would throw an error in statistical test and remove those taxa from DF -
        # first find the taxa that would throw an error in t.test or wilcox.test
        var_plus_length_check <- group_by(TbT_DF_l, Taxon, Group) %>% summarise(Variance = var(Ratio, na.rm = T), NotNA = sum(!is.na(Ratio)))
        if (test == "t.test"){
                var_plus_length_check <- filter(var_plus_length_check, !(Variance > 0) | NotNA < 2)
        } else if (test == "wilcox.test") {
                var_plus_length_check <- filter(var_plus_length_check, NotNA < 1)
        }
        
        if (nrow(var_plus_length_check) != 0) { # should always remove at least the case where taxa_denominator = taxa_nom, since those ratios are all 1, hence no variance
                TbT_DF_l <- filter(TbT_DF_l, !(Taxon %in% unique(var_plus_length_check$Taxon)))
        }
        # --
        
        # - use ggpubr::compare_menas to calculate all pValues of the Ratios for the different taxa_den between current group_levels -
        pVals <- ggpubr::compare_means(formula = Ratio ~ Group, data = TbT_DF_l, group.by = "Taxon", method = test, p.adjust.method = p_adjust_method, symnum.args = symnum.args)
        
        pVals <- dplyr::arrange(pVals, p)
        
        # NB: the pVals change when you log the ratios (scale_y_log10())
        TbT_DF_l$RatioLog10 <- log10(TbT_DF_l$Ratio)
        pValsLog <- ggpubr::compare_means(formula = RatioLog10 ~ Group, data = TbT_DF_l, group.by = "Taxon", method = test, p.adjust.method = p_adjust_method, symnum.args = symnum.args)
        pValsLog <- dplyr::arrange(pValsLog, p)
        # --
        
        # NB: I plot now for log and non log independently, even though ggpubr is so smart to change the p-values when you just use
        # Tr + scale_y_log10(). I plot independently because log might change the order!
        TbT_DF_l_log <- TbT_DF_l
        
        # - order taxa based on pVals result -
        TbT_DF_l$Taxon <- factor(TbT_DF_l$Taxon, levels = unique(pVals$Taxon), ordered = TRUE)
        TbT_DF_l_log$Taxon <- factor(TbT_DF_l_log$Taxon, levels = unique(pValsLog$Taxon), ordered = TRUE)
        # --
        
        # - since you might have more than two levels in each plot you need to set the comparisons argument in stat_compare_means -
        comparisonList <- get_unique_facLevel_combinations(levels(group_fac))
        # --
        
        # - plot: NB: also for non removed taxa some samples might have NA ratios that will be removed -
        
        Tr <- ggplot(TbT_DF_l, aes(x = Group, y = Ratio, col = Group))
        Tr <- Tr +
                geom_violin() +
                geom_point(size = 1, alpha = 0.6, position = position_jitterdodge(dodge.width = 1)) +
                # scale_color_manual(values = c(color_lookup$color[i], color_lookup$color[j])) +
                scale_color_manual(values = color_levels) +
                facet_wrap(~ Taxon, scales = "free_y") +
                xlab("") +
                ylab(paste("abundance ratio of", taxa_nom, "to stated taxon")) +
                theme_bw() +
                theme(legend.position = "none")
        
        Tr <- Tr + ggpubr::stat_compare_means(comparisons = comparisonList, label = "p.signif", method = "t.test", hide.ns = TRUE)
        
        Tr1 <- ggplot(TbT_DF_l, aes(x = Group, y = Ratio, col = Group))
        Tr1 <- Tr1 +
                geom_boxplot(outlier.color = NA) +
                geom_point(size = 1, alpha = 0.6, position = position_jitterdodge(dodge.width = 1)) +
                # scale_color_manual(values = c(color_lookup$color[i], color_lookup$color[j])) +
                scale_color_manual(values = color_levels) +
                facet_wrap(~ Taxon, scales = "free_y") +
                xlab("") +
                ylab(paste("abundance ratio of", taxa_nom, "to stated taxon")) +
                theme_bw() +
                theme(legend.position = "none")
        
        
        Tr1 <- Tr1 + ggpubr::stat_compare_means(comparisons = comparisonList, label = "p.signif", method = test, hide.ns = TRUE) # p.format
        
        
        Tr2 <- ggplot(TbT_DF_l_log, aes(x = Group, y = Ratio, col = Group))
        Tr2 <- Tr2 +
                geom_violin() +
                geom_point(size = 1, alpha = 0.6, position = position_jitterdodge(dodge.width = 1)) +
                # scale_color_manual(values = c(color_lookup$color[i], color_lookup$color[j])) +
                scale_color_manual(values = color_levels) +
                facet_wrap(~ Taxon, scales = "free_y") +
                xlab("") +
                ylab(paste("abundance ratio of", taxa_nom, "to stated taxon")) +
                theme_bw() +
                theme(legend.position = "none")
        
        
        Tr2 <- Tr2 + ggpubr::stat_compare_means(comparisons = comparisonList, label = "p.signif", method = test, hide.ns = TRUE)
        Tr2 <- Tr2 + scale_y_log10()
        
        
        Tr3 <- ggplot(TbT_DF_l_log, aes(x = Group, y = Ratio, col = Group))
        Tr3 <- Tr3 +
                geom_boxplot(outlier.color = NA) +
                geom_point(size = 1, alpha = 0.6, position = position_jitterdodge(dodge.width = 1)) +
                # scale_color_manual(values = c(color_lookup$color[i], color_lookup$color[j])) +
                scale_color_manual(values = color_levels) +
                facet_wrap(~ Taxon, scales = "free_y") +
                xlab("") +
                ylab(paste("abundance ratio of", taxa_nom, "to stated taxon")) +
                theme_bw() +
                theme(legend.position = "none")
        
        Tr3 <- Tr3 + ggpubr::stat_compare_means(comparisons = comparisonList, label = "p.signif", method = test, hide.ns = TRUE)
        Tr3 <- Tr3 + scale_y_log10()
        
        
        list(pVals = pVals, Tr = Tr, Tr1 = Tr1, pValsLog = pValsLog, Tr2 = Tr2, Tr3 = Tr3)
        
}


####################################
## plot_sampleSums_FirmtoBac 
###################################
# to check if there is a trend that sequencing size affects F to B ratios

## Output: 
# - list(pVals = pVals, Tr = Tr, Tr1 = Tr1, pValsLog = pValsLog, Tr2 = Tr2, Tr3 = Tr3)


plot_sampleSums_FirmtoBac <- function(physeq, group_var, tax_names = NULL,
                                       taxa_nom = "Firmicutes", taxa_den = NULL, color_levels, test = "t.test", p_adjust_method = "fdr",
                                       symnum.args = list(cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))) {
        
        
        # - check that given tax_names fit to physeq and change taxa_names of physeq -
        if (is.null(tax_names)){
                tax_names <- paste("T", 1:length(TbTmatrixes), sep = "_")
        } 
        
        if(!identical(ntaxa(physeq), length(tax_names))){stop("tax_names do not fit in length to physeq")}
        
        tax_names <- make.unique(tax_names)
        taxa_names(physeq) <- tax_names
        # --
        
        
        if (taxa_are_rows(physeq)) {physeq <- t(physeq)}
        
        # - calculate the matrix taxa_nom/(all other taxa) -         
        CT <- t(as(otu_table(physeq), 'matrix')) # now taxa are rows and samples are columns
        
        i <- which(rownames(CT) == taxa_nom)
        if (length(i) != 1) {stop("taxa_nom not found in tax_names or tax_names not unique!")}
        
        
        TbTmatrix <- apply(CT, 2, function(samp_cnts){samp_cnts[i]/samp_cnts})
        # produces for each taxon (= host taxon) a TbTMatrix
        # NB: there are possibly Inf, and NaN values in the matrix, specifically
        # 0/x = 0, x/0 = Inf; 0/0 = NaN!
        # --
        
        TbT_DF <- as.data.frame(TbTmatrix)
        TbT_DF$Taxon <- rownames(TbT_DF)
        # - use taxa_den (denominator) to restrict the taxa to which taxa_nom is compared to -
        if (is.null(taxa_den)) {taxa_den <- tax_names}
        TbT_DF <- TbT_DF[TbT_DF$Taxon %in% taxa_den, ]
        # --
        
        # - change to long DF -
        TbT_DF_l <- gather(TbT_DF, key = Sample, value = Ratio, -Taxon)
        # --
        
        # - add the group_var level information and the total_counts -
        LookUpDF <- data.frame(Sample = sample_names(physeq), Group = sample_data(physeq)[[group_var]], total_count = sample_sums(physeq))
        TbT_DF_l$Group <- as.character(LookUpDF$Group[match(TbT_DF_l$Sample, LookUpDF$Sample)])
        group_fac <- factor(sample_data(physeq)[[group_var]])
        TbT_DF_l$Group <- factor(TbT_DF_l$Group, levels = levels(group_fac), ordered = T)
        TbT_DF_l$total_count <- LookUpDF$total_count[match(TbT_DF_l$Sample, LookUpDF$Sample)]
        # --
        
        # - change all ratios where either nominator taxon or denominator taxon had count = 0 to NA -
        # remember: 0/0 = NaN (not finite), 0/x = 0, x/0 = Inf
        TbT_DF_l$Ratio[!is.finite(TbT_DF_l$Ratio) | TbT_DF_l$Ratio == 0] <- NA
        # --
        
        
        # - plot: NB: also for non removed taxa some samples might have NA ratios that will be removed -
        
        Tr <- ggplot(TbT_DF_l, aes(x = total_count, y = Ratio))
        Tr <- Tr +
                geom_point(size = 1, alpha = 0.6, aes(col = Group)) +
                # scale_color_manual(values = c(color_lookup$color[i], color_lookup$color[j])) +
                scale_color_manual("", values = color_levels) +
                facet_wrap(~ Taxon, scales = "free_y") +
                xlab("sample_sums()") +
                ylab(paste("abundance ratio of", taxa_nom, "to stated taxon")) +
                theme_bw()
        
        
        
        Tr
        
}





#######################################
#### plot_bar_own
#######################################
# based on plot_bar from phyloseq, difference, orders and colors the Samples based on group_var, and orders the fill based on abundance
# I guess inputs can be guessed on

plot_bar_own <- function(physeq, x = "Sample", y = "Abundance", group_var, color_levels, fill = NULL,
                         color_sample_names = TRUE, facet_grid = NULL){
        
        if(taxa_are_rows(physeq)) { physeq <- t(physeq) }
        
        if (!is.factor(sample_data(physeq)[[group_var]])) {sample_data(physeq)[[group_var]] <- as.factor(sample_data(physeq)[[group_var]])}
        
        if (is.null(fill)) { fill = "Phylum"}
        
        mdf <- phyloseq::psmelt(physeq)
        
        # order samples according to levels
        LookUpDF <- data.frame(Sample = sample_names(physeq), Group = sample_data(physeq)[[group_var]])
        LookUpDF <- LookUpDF[order(match(LookUpDF$Group, levels(LookUpDF$Group))), ]
        mdf$Sample <- factor(mdf$Sample, levels = LookUpDF$Sample, ordered = TRUE)
        # mdf$Sample <- factor(mdf$Sample, levels = c("A-15A", "A-5A", "A-2A", "A-1A", "B-15A", "B-5A", "B-2A", "B-1A"), ordered = TRUE)
        
        # order fill levels according to abundance over all samples
        mdf[, fill] <- as.character(mdf[, fill])
        mdf[is.na(mdf[, fill]), fill] <- "NA"
        sums <- group_by_(mdf, fill) %>% summarise(sum_abundance = sum(Abundance)) %>% arrange(sum_abundance)
        mdf[, fill] <- factor(mdf[, fill], levels = as.character(sums[[1]]), ordered = TRUE)
        

        # - define names of x axis using color_levels (which must be a named character vector) -        
        colxaxis <- color_levels[LookUpDF$Group]
        # --
        
        
        
        if (length(levels(mdf[, fill])) <= 15) {
                fill_colors <- QuantColors15[1:length(levels(mdf[, fill]))]
                names(fill_colors) <- rev(levels(mdf[, fill]))
        } else {
                fill_colors <- rev(viridis(length(levels(mdf[, fill]))))
                names(fill_colors) <- rev(levels(mdf[, fill]))
        }
        
        
        Tr <- ggplot(mdf, aes_string(x = x, y = y, fill = fill))
        Tr <- Tr + 
                geom_bar(stat = "identity", position = "stack") +
                theme_bw() +
                scale_fill_manual(values = fill_colors) +
                xlab("") +
                theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0, colour = colxaxis))
        
        if (!is.null(facet_grid)) {
                Tr <- Tr + facet_grid(facet_grid)
        }
        
        
        Tr
}


#######################################
#### plot_bar_own_compare
#######################################
# generates abundance barplots (see plot_bar_own) to compare ps to ps_tca, i.e. to see how SFs adjustment affects the abundances 

plot_bar_own_compare <- function(physeq, physeq2, x = "Sample", y = "Abundance", group_var, color_levels, fill = NULL,
                         color_sample_names = TRUE){
        
        # - prepare mdf of ps physeq -
        if(taxa_are_rows(physeq)) { physeq <- t(physeq) }
        
        if (!is.factor(sample_data(physeq)[[group_var]])) {sample_data(physeq)[[group_var]] <- as.factor(sample_data(physeq)[[group_var]])}
        
        if (is.null(fill)) { fill = "Phylum"}
        
        mdf <- phyloseq::psmelt(physeq)
        
        # order samples according to levels
        LookUpDF <- data.frame(Sample = sample_names(physeq), Group = sample_data(physeq)[[group_var]])
        LookUpDF <- LookUpDF[order(match(LookUpDF$Group, levels(LookUpDF$Group))), ]
        mdf$Sample <- factor(mdf$Sample, levels = LookUpDF$Sample, ordered = TRUE)
        # mdf$Sample <- factor(mdf$Sample, levels = c("A-15A", "A-5A", "A-2A", "A-1A", "B-15A", "B-5A", "B-2A", "B-1A"), ordered = TRUE)
        
        # order fill levels according to abundance over all samples
        mdf[, fill] <- as.character(mdf[, fill])
        mdf[is.na(mdf[, fill]), fill] <- "NA"
        sums <- group_by_(mdf, fill) %>% summarise(sum_abundance = sum(Abundance)) %>% arrange(sum_abundance)
        mdf[, fill] <- factor(mdf[, fill], levels = as.character(sums[[1]]), ordered = TRUE)
        
        # --
        
        # - prepare mdf of ps_tca physeq -
        if(taxa_are_rows(physeq2)) { physeq2 <- t(physeq2) }
        
        if (!is.factor(sample_data(physeq2)[[group_var]])) {sample_data(physeq2)[[group_var]] <- as.factor(sample_data(physeq2)[[group_var]])}
        
        if (is.null(fill)) { fill = "Phylum"}
        
        mdf2 <- phyloseq::psmelt(physeq2)
        
        # order samples according to levels
        LookUpDF <- data.frame(Sample = sample_names(physeq2), Group = sample_data(physeq2)[[group_var]])
        LookUpDF <- LookUpDF[order(match(LookUpDF$Group, levels(LookUpDF$Group))), ]
        mdf2$Sample <- factor(mdf2$Sample, levels = LookUpDF$Sample, ordered = TRUE)
        # mdf2$Sample <- factor(mdf2$Sample, levels = c("A-15A", "A-5A", "A-2A", "A-1A", "B-15A", "B-5A", "B-2A", "B-1A"), ordered = TRUE)
        
        # order fill levels according to abundance over all samples
        mdf2[, fill] <- as.character(mdf2[, fill])
        mdf2[is.na(mdf2[, fill]), fill] <- "NA"
        sums <- group_by_(mdf2, fill) %>% summarise(sum_abundance = sum(Abundance)) %>% arrange(sum_abundance)
        mdf2[, fill] <- factor(mdf2[, fill], levels = as.character(sums[[1]]), ordered = TRUE)
        # --
        
        mdf$Typer <- "before"
        mdf2$Typer <- "after SF adjustment"
        
        mdf <- rbind(mdf, mdf2)
        mdf$Typer <- factor(mdf$Typer, levels = c("before", "after SF adjustment"), ordered = TRUE)
        
        # - define names of x axis using color_levels (which must be a named character vector) -        
        colxaxis <- color_levels[LookUpDF$Group]
        # --
        
        
        
        if (length(levels(mdf[, fill])) <= 15) {
                fill_colors <- QuantColors15[1:length(levels(mdf[, fill]))]
                names(fill_colors) <- rev(levels(mdf[, fill]))
        } else {
                fill_colors <- rev(viridis(length(levels(mdf[, fill]))))
                names(fill_colors) <- rev(levels(mdf[, fill]))
        }
        
        
        Tr <- ggplot(mdf, aes_string(x = x, y = y, fill = fill))
        Tr <- Tr + 
                geom_bar(stat = "identity", position = "stack") +
                theme_bw() +
                scale_fill_manual(values = fill_colors) +
                xlab("") +
                facet_grid(~ Typer) +
                theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0, colour = colxaxis))
        
        Tr
}



####################################
## create_raw_TbT_TilePlots: 
###################################
# wilcoxon.test is used
## Input: 
# - TbTmatrixes_list: The list with the lists of TbTmatrixes for each level combi in group factor
# NB: here it should be raw_TbTmatrixes with 0 = 0/x, Inf = x/0, NaN = 0/0, all these values will be ignored!!
# - physeq: used for TbTmatrixes_list generation
# - p_adjust: method for p.adjust, NB: if not bonferroni or none, the tile plots are not necessarily symmetric anymore
## Output: 
# - list of pValMatrixes plus TileTr for each level combination. 
# NB: I negated the p-values if host taxon was more abundant in grp2 compared to other taxon!!

create_raw_TbT_TilePlots <- function(TbTmatrixes_list, physeq, group_var, tax_names = NULL,
                                     test = "wilcoxon", p_adjust = "none") {
        
        if(!identical(length(TbTmatrixes_list[[1]]), ntaxa(physeq))){stop("TbTmatrixes can not fit to physeq")}
        
        if(test != "wilcoxon" & test != "t.test"){stop("test unknown, must be wilcoxon or t.test")}
        
        group_fac <- factor(sample_data(physeq)[[group_var]])
        
        fac_levels <- levels(group_fac)
        
        # - get the level combis -
        fac_levels_num <- setNames(seq_along(fac_levels), fac_levels) 
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
        # --
        
        result_list <- vector("list", length = length(i_s))
        
        for (k in seq_along(i_s)) {
                
                i <- i_s[k]
                j <- j_s[k]
                
                TbTmatrixes <- TbTmatrixes_list[[k]]
                
                if (is.null(tax_names)){
                        tax_names <- paste("T", 1:length(TbTmatrixes), sep = "_")
                } else {
                        if(!identical(length(TbTmatrixes), length(tax_names))){stop("tax_names do not fit in length to TbTmatrixes")}
                }
                
                tax_names <- make.unique(tax_names)
                names(TbTmatrixes) <- tax_names
                
                TbTmatrixes <- lapply(TbTmatrixes, function(mat){
                        rownames(mat) <- tax_names
                        mat
                })
                group_fac_current <- droplevels(group_fac[as.numeric(group_fac) %in% c(i, j)])
                
                # ntaxa * ntaxa wilcoxon tests take time!
                pValMatrix <- sapply(TbTmatrixes, function(mat){
                        apply(mat, 1, function(taxon_ratios){
                                x <- taxon_ratios[group_fac_current == fac_levels[i]]
                                x <- x[is.finite(x) & x != 0] # removes all ratios in which one of the two taxa was not present!
                                y <- taxon_ratios[group_fac_current == fac_levels[j]]
                                y <- y[is.finite(y) & y != 0] # removes 0/0 = NaN, 0/x = 0, x/0 = Inf
                                if (test == "wilcoxon"){
                                        if (length(x) > 0 && length(y) > 0){
                                                pValue <- wilcox.test(x = x, y = y, alternative = "two", paired = F, exact = F)$p.value
                                                # NB: wilcox.test ignores in default setting (maybe see na.action) NA, NaN, Inf, -Inf
                                                # For plot: change sign of pValue to negative if taxon is more abundant in group 1. 
                                                Ranks <- rank(c(x[!is.na(x)], y[!is.na(y)]))
                                                n1 <- length(x[!is.na(x)])
                                                n2 <- length(y[!is.na(y)])
                                                Wx <- sum(Ranks[1:n1])-(n1*(n1+1)/2)
                                                Wy <- sum(Ranks[(n1+1):(n1+n2)])-(n2*(n2+1)/2)
                                                if(Wx > Wy){pValue <- -1*pValue}
                                                pValue
                                                
                                        } else {
                                                pValue = 1
                                        }
                                        
                                } else if (test == "t.test") {
                                        if (length(x) > 1 && length(y) > 1 && var(x) > 0 && var(y) > 0){
                                                pValue <- t.test(x = x, y = y, alternative = "two")$p.value
                                                if (mean(x, na.rm = T) > mean(y, na.rm = T)){pValue <- -1*pValue}
                                                pValue
                                        } else {
                                                pValue <- 1
                                        }
                                        
                                }
                                
                        })
                })
                
                # make sure diagonal is all NA (can be exceptions especially for t.test)
                diag(pValMatrix) <- NA
                
                signs <- pValMatrix < 0
                signs[is.na(signs)] <- FALSE
                
                pValMatrix <- abs(pValMatrix)
                for (e in 1:nrow(pValMatrix)){
                        pValMatrix[e, ] <- p.adjust(pValMatrix[e, ], method = p_adjust)
                } # equal to t(apply(pValMatrix, 1, p.adjust, method = p_adjust))
                
                pValMatrix[signs] <- pValMatrix[signs]*(-1)
                
                # -- add a tile plot of the pValMatrix --
                
                DF <- as.data.frame(pValMatrix)
                DF[is.na(DF)] <- 2 # just to avoid missing values in plot and have a clear non-pValue value to mark self comparisons as black
                DF$HostTaxon <- rownames(pValMatrix)
                DF <- tidyr::gather(DF, key = Taxon , value = pValue, - HostTaxon)
                DF$Taxon <- factor(DF$Taxon, levels = rownames(pValMatrix), ordered = TRUE)
                DF$HostTaxon <- factor(DF$HostTaxon, levels = rev(rownames(pValMatrix)), ordered = TRUE)
                
                # # add color to the taxa names so you see up and down
                # colyaxis <- vector(mode = "character", length = nrow(DF))
                # colyaxis[] <- "black"
                # colyaxis[grepl("TP-U", levels(DF$HostTaxon))] <- "#E69F00"
                # colyaxis[grepl("TP-D", levels(DF$HostTaxon))] <- "#009E73"
                # colxaxis <- vector(mode = "character", length = nrow(DF))
                # colxaxis[] <- "black"
                # colxaxis[grepl("TP-U", levels(DF$Taxon))] <- "#E69F00"
                # colxaxis[grepl("TP-D", levels(DF$Taxon))] <- "#009E73"
                DF$Fill <- "not significant"
                DF$Fill[DF$pValue < 0.05 & DF$pValue > 0] <- "up (p < 0.05)"
                DF$Fill[DF$pValue > -0.05 & DF$pValue < 0] <- "down (p < 0.05)"
                DF$Fill[DF$pValue == 2] <- " "
                DF$Fill <- factor(DF$Fill, levels = c("up (p < 0.05)", "not significant", "down (p < 0.05)", " "), ordered = T)
                TileTr <- ggplot(DF, aes(x = Taxon, y = HostTaxon, fill = Fill))
                TileTr <- TileTr + 
                        geom_raster() + 
                        ggtitle(names(TbTmatrixes_list)[k]) +
                        scale_fill_manual("", values = c("not significant" = "gray98", "up (p < 0.05)" = "#AA4499", "down (p < 0.05)" = "#88CCEE", " " = "black")) +
                        scale_x_discrete(position = "top") +
                        labs(x=NULL, y=NULL) +
                        theme_tufte(base_family="Helvetica") +
                        theme(plot.title=element_text(hjust=0)) +
                        theme(axis.ticks=element_blank()) +
                        theme(axis.text=element_text(size=7)) +
                        theme(legend.title=element_blank()) +
                        theme(legend.text=element_text(size=6)) +
                        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0)) # colour = colxaxis
                result_list[[k]] <- list(pValMatrix = pValMatrix, TileTr = TileTr)
        }
        
        names(result_list) <- names(TbTmatrixes_list)
        result_list
}

