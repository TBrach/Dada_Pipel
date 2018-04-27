draw_taxa_heatmap <- function (physeq, taxa_data, group, compare = NULL, block = NULL, 
          max_abundance_for_color = NULL, log_transform = FALSE, border_color = "grey60", 
          filename = NA, show_rownames = TRUE, gradient_steps = c(0.01, 
                                                                  0.1, 0.5, 1), cellwidth = NA, cellheight = NA, gaps_col = NULL, 
          custom_palette = NULL) 
{
        if (!is.null(taxa_data)) {
                shown_taxa <- rownames(taxa_data)
                pruned_ps <- prune_taxa(shown_taxa, physeq)
        }
        else {
                shown_taxa <- taxa_names(physeq)
                pruned_ps <- physeq
        }
        tax_table <- tax_table(pruned_ps)
        tax_labels <- apply(tax_table, 1, get_pretty_taxon_name)
        if (length(tax_labels) != length(unique(tax_labels))) {
                rownames_bak <- names(tax_labels)
                tax_labels <- paste0(tax_labels, " [", rownames(tax_table), 
                                     "]")
                names(tax_labels) <- rownames_bak
                rm(rownames_bak)
        }
        my_heat_map <- make_heat_map_from_physeq(pruned_ps, group = group, 
                                                 compare = compare, block = block, max_abundance_for_color = max_abundance_for_color, 
                                                 log_transform = log_transform, border_color = border_color, 
                                                 tax_labels = tax_labels, taxa_data = taxa_data, tax_order = shown_taxa, 
                                                 show_rownames = show_rownames, cellwidth = cellwidth, 
                                                 cellheight = cellheight, gaps_col = gaps_col, gradient_steps = gradient_steps, 
                                                 filename = filename, custom_palette = custom_palette)
        return(my_heat_map)
}
<bytecode: 0x124e06e98>
        <environment: namespace:microbiomeX>