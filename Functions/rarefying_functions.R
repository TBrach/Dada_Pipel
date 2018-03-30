cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#######################################
### rarefy_sample
#######################################
# simple rarefy function based on sample(), see phyloseq:::rarefaction_subsample

rarefy_sample <- function (sample_cnts, size) {
        # attention, only works if nrow(ProbMat) is not NULL
        simsample <- integer(length(sample_cnts)) 
        suppressWarnings(draws <- sample(1:length(sample_cnts), size = size, replace = TRUE, prob = sample_cnts))
        drawtable <- table(draws)
        simsample[as(names(drawtable), "integer")] <- drawtable
        return(simsample)
}


#######################################
### raref_curve_richness
#######################################
# read rarefaction_curve_own, this fast version uses the vegan::rarefy function, so it only generates rarefaction
# for richness. NB: vegan::rarefy used here does averaging (read help)
# and is super fast. I still do not use the SE provided by rarefy, so that could be added in here

# # NB: in case you wanted to add sample labels to a Tr_richness_col
# df <- Tr_richness_col$data
# df <- filter(df, step == max(df$step))
# df$step <- df$step + 400
# Tr_richness_col <- Tr_richness_col + 
#         geom_label(data = df, aes(label = Sample)) +
#         scale_x_continuous(limits = c(-100, max(df$step) + 2500))

raref_curve_richness <- function(physeq, group_var = NULL, max_total = NULL, step_size = 200, col_levels, seed = 123) {
        
        
        if (taxa_are_rows(physeq)) {
                seqtab <- t(as(otu_table(physeq), "matrix"))
        } else {
                seqtab <- as(otu_table(physeq), "matrix")
        }
        
        
        if (is.null(max_total)) {
                max_total <- quantile(sample_sums(physeq), probs = .25)
        }
        
        steps <- seq(from = 0, to = max_total, by = step_size)
        NoSamples <- nsamples(physeq)
        NoSteps <- length(steps)
        richness_matrix <- matrix(data = NA, nrow = NoSamples, ncol = NoSteps)
        totalamplicons <- sample_sums(physeq)
        
        set.seed(seed)
        
        # ptm <- Sys.time()
        for (i in 1:length(steps)) {
                richness_matrix[,i] <- vegan::rarefy(seqtab, sample = steps[i], se = FALSE)
        }
        # Sys.time() - ptm
        
        # set richness for samples at which steps > totalamplicons to NA
        for (i in 1:NoSamples) {
                NaIndex <- which(totalamplicons[i] < steps)[1]
                if (!is.na(NaIndex)){
                        richness_matrix[i, NaIndex:ncol(richness_matrix)] <- NA
                }
        }
        
        rownames(richness_matrix) <- rownames(seqtab)
        colnames(richness_matrix) <- paste("step_", steps, sep = "")
        richness_df <- as.data.frame(richness_matrix)
        
        
        plot_div_df3 <- function (div_df, type = "richness") {
                
                div_df$Sample <- rownames(div_df)
                div_df$Total <- totalamplicons
                div_df <- tidyr::gather(div_df, key = step, value = div, -Sample, -Total)
                div_df <- tidyr::separate(div_df, col = step, into = c("term", "step"), sep = "_")
                div_df$step <- as.numeric(div_df$step)
                
                Tr <- ggplot(div_df, aes(x = step, y = div, group = Sample, col = Total))
                Tr <- Tr +
                        geom_line() +
                        scale_color_gradient2("total counts bef.", low = cbPalette[6], mid = cbPalette[1], high = cbPalette[2], midpoint = median(div_df$Total)) +
                        xlab("total counts in sample") +
                        ylab(type) +
                        theme_bw()
        }
        
        Tr_richness_grad <- plot_div_df3(richness_df, type = "richness")
        
        
        
        if (!is.null(group_var)){
                Group <- sample_data(physeq)[[group_var]] 
                if (!is.null(Group) && !is.factor(Group)) {
                        Group <- as.factor(Group)
                }
        } else { Group <- NULL}
        
        
        if (!is.null(Group)) {
                
                richness_df$Group <- Group
                
                pairwise.tt_richness <- lapply(richness_df[, -ncol(richness_df)], function(step){ 
                        ptt <- pairwise.t.test(x = step, g = richness_df$Group, alternative = "two", p.adjust.method = "none", var.equal = F, pool.sd = F)
                        ptt$p.value
                })
                
                plot_div_df2 <- function (div_df, type = "richness") {
                        
                        div_df$Sample <- rownames(div_df)
                        div_df <- tidyr::gather(div_df, key = step, value = div, -Sample, -Group)
                        div_df <- tidyr::separate(div_df, col = step, into = c("term", "step"), sep = "_")
                        div_df$step <- as.numeric(div_df$step)
                        
                        Tr <- ggplot(div_df, aes(x = step, y = div, group = Sample, col = Group))
                        Tr <- Tr +
                                geom_line() +
                                xlab("total counts in sample") +
                                scale_color_manual("", values = color_levels) +
                                ylab(type) +
                                theme_bw(12)
                }
                
                plot_div_df_group <- function (div_df, type = "alpha diversity") {
                        div_df <- tidyr::gather(div_df, key = step, value = div, - Group)
                        div_df <- tidyr::separate(div_df, col = step, into = c("term", "step"), sep = "_")
                        div_df$step <- as.numeric(div_df$step)
                        div_df <- group_by(div_df, Group, step)
                        div_df <- dplyr::summarise(div_df, Mean_div = mean(div, na.rm = T), SD_div = sd(div, na.rm = T), n = n(), SE_div = SD_div/sqrt(n))
                        
                        Tr <- ggplot(div_df, aes(x = step, y = Mean_div, col = Group))
                        Tr <- Tr + 
                                geom_line() +
                                geom_point(size = 1) +
                                geom_errorbar(aes(ymin = Mean_div-SE_div, ymax = Mean_div+SE_div)) +
                                scale_color_manual("", values = color_levels) +
                                xlab("total counts in sample") +
                                ylab(paste(type, "(group mean +- SE)")) +
                                theme_bw(12)
                }
                Tr_richness_col <- plot_div_df2(richness_df, type = "richness")
                
                Tr_richness_group <- plot_div_df_group(richness_df, type = "richness")
                
        } else {
                Tr_richness_col <- NA
                Tr_richness_group <- NA
                pairwise.tt_richness <- NA
        }
        
        outlist <- list(rarefaction_richness_df = richness_df, Tr_richness_grad = Tr_richness_grad, Tr_richness_col = Tr_richness_col,
                        Tr_richness_group = Tr_richness_group, pairwise.tt_richness = pairwise.tt_richness)
        
        
}

