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
        
        # - get the level combis -
        fac_levels_num <- setNames(seq_along(fac_levels), fac_levels) 
        i_sAndj_s <- get_unique_facLevel_combis(fac_levels)
        i_s <- i_sAndj_s[["i_s"]]
        j_s <- i_sAndj_s[["j_s"]]
        # --
        
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