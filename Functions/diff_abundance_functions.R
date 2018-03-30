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
                colnames(p_val_df)[1:2] <- paste(c("p_val_", "p_val_adj_"), fac_levels[i], "_vs_", fac_levels[j], sep = "")
                colnames(p_val_df)[5:6] <- paste(c("prev_PC_"), c(fac_levels[i], fac_levels[j]), sep = "")
                p_val_df <- cbind(as.data.frame(p_val_df), tax_table(physeq))
                p_val_df$Taxon <- colnames(CT)
                p_val_df <- arrange(p_val_df, p_val_df[,1])
                p_val_df <- select(p_val_df, Taxon, 1:(ncol(p_val_df)-1))
                p_val_list[[k]] <- p_val_df
                names(p_val_list)[[k]] <- paste(fac_levels[i], "_vs_", fac_levels[j], sep = "")
        }
        
        p_val_list
        
}