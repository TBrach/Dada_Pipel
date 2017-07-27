cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#######################################
### FUNCTION: assignTaxonomyaddSpecies
#######################################

# Function that simply calls assignTaxonomy and then addSpecies from dada2
## Input
# seqtab: your sequence table from dada2
# minBoot: minBoot for assignTaxonomy
# allowMultiple: default 3, allowMultiple for addSpecies
# PathToRefs: path to the folder with bot the RefDataBase and SpeciesDB in
# RefDataBase: the reference database for assignTaxonomy
# SpeciesDB: the reference database for addSpecies
# PathToSave: the folder into which taxa and taxa.species will be saved, if NULL working directory is used
## Output:
# taxa, taxa.species are saved as Taxonomy.RData


assignTaxonomyaddSpecies <- function(seqtab, 
                                     minBoot = 50, 
                                     allowMultiple = TRUE,
                                     PathToRefs = NULL,
                                     RefDataBase = "silva_nr_v128_train_set.fa.gz",
                                     SpeciesDB = "silva_species_assignment_v128.fa.gz",
                                     PathToSave = getwd(),
                                     tryRC = FALSE){
        
        ## Packages
        try(library(dada2), biocLite("dada2"))
        
        
        ## Reference Databases
        if(is.null(PathToRefs)){
                stop("PathToRefs must be given")
        }
        
        RefDB <- file.path(PathToRefs, RefDataBase)
        
        if(!file.exists(RefDB)){
                stop("could not find the reference database for assignTaxonomy")
        }
        
        SpecDB <- file.path(PathToRefs, SpeciesDB)
        
        if(!file.exists(SpecDB)){
                stop("could not find the reference database for addSpecies")
        }
        
        ## list the input
        InputSave <- list(seqtab = seqtab, minBoot = minBoot,
                          allowMultiple = allowMultiple,
                          tryRC = tryRC,
                          PathToRefs = PathToRefs, 
                          RefDataBase = RefDataBase,
                          SpeciesDB = SpeciesDB,
                          PathToSave = PathToSave,
                          tryRC = tryRC)
        
        
        ## run assignTaxonomy and save
        ptm <- proc.time()
        taxa <- assignTaxonomy(seqtab, refFasta = RefDB, verbose = TRUE, minBoot = minBoot, tryRC = tryRC)
        TimeForassignTaxonomy <- (proc.time()-ptm)[3]
        
        InputSave[[length(InputSave) + 1]] <- TimeForassignTaxonomy
        
        save(taxa, InputSave, file = file.path(PathToSave, "Taxonomy.RData"))
        
        message("assignTaxonomy has been run")
        
        ## run addSpecies and save
        ptm <- proc.time()
        taxa.species <- addSpecies(taxa, refFasta = SpecDB, verbose = TRUE, allowMultiple = allowMultiple)
        TimeForaddSpecies <- (proc.time()-ptm)[3]
        
        InputSave[[length(InputSave) + 1]] <- TimeForaddSpecies
        
        save(taxa, taxa.species, InputSave, file = file.path(PathToSave, "Taxonomy.RData"))


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
### FUNCTION: check_assignment_vs_prevalence
#######################################
# checks if there is a trend for better assignment for more prealent SVs

check_assignment_vs_prevalence <- function(taxa, seqtab, prevalences = seq(0, 90, by = 10)){
        
        prev <- colSums(seqtab != 0)
        remaining_SVs <- lapply(prevalences, function(preva) {
                prev > (preva/100)*nrow(seqtab)
        })
        No_SVs <- sapply(remaining_SVs, sum)
        filtered.taxas <- lapply(remaining_SVs, function(indexes){
                taxa[indexes,]
        })
        
        assignment_distributions <- lapply(filtered.taxas, get_assignemnt_distribution)
        
        PC_assigned <- sapply(assignment_distributions, function(distri){distri[["PC_assigned"]]})
        
        rownames(PC_assigned) <- colnames(taxa)
        colnames(PC_assigned) <- paste("Prev_", prevalences, "_", No_SVs, sep = "")
        
        PC_assigned <- as.data.frame(PC_assigned)
        PC_assigned$Level <- rownames(PC_assigned)
        df_long <- gather(PC_assigned, key = Prev, value = PC, -Level)
        df_long <- separate(df_long, col = Prev, into = c("Type", "Prevalence", "No_SVs"), sep = "_")
        df_long$Level <- factor(df_long$Level, levels = colnames(taxa), ordered = TRUE)
        df_long$Prevalence <- factor(df_long$Prevalence, levels = as.character(prevalences), ordered = TRUE)
        
        
        tr <- ggplot(df_long, aes(x = Prevalence, y = PC, col = Level))
        tr <- tr + 
                geom_point(size = 2) +
                scale_colour_manual("", values = cbPalette[2:8]) +
                scale_x_discrete(breaks = prevalences, labels = paste(prevalences, " (", No_SVs, ")", sep = "")) +
                ylab("percentage assigned") +
                xlab("prevalence (No SVs)") +
                theme_bw(12) +
                theme(axis.text.x = element_text(angle=90, vjust=0.5))
        
        list(PC_assigned, tr)
        
}


#######################################
### FUNCTION: check_assignment_vs_abundance
#######################################
# checks if there is a trend for better assignment for more abundant SVs

check_assignment_vs_abundance <- function(taxa, seqtab, abundanceQuantiles = seq(0, 90, by = 10)){
        
        abundances <- colSums(seqtab)
        abQuantiles <- quantile(abundances, probs = abundanceQuantiles/100)
        remaining_SVs <- lapply(abQuantiles, function(quant) {
                abundances >= quant
        })
        No_SVs <- sapply(remaining_SVs, sum)
        filtered.taxas <- lapply(remaining_SVs, function(indexes){
                taxa[indexes,]
        })
        
        assignment_distributions <- lapply(filtered.taxas, get_assignemnt_distribution)
        
        PC_assigned <- sapply(assignment_distributions, function(distri){distri[["PC_assigned"]]})
        
        rownames(PC_assigned) <- colnames(taxa)
        colnames(PC_assigned) <- paste("Ab_", abundanceQuantiles, "_", round(abQuantiles), "_", No_SVs, sep = "")
        
        PC_assigned <- as.data.frame(PC_assigned)
        PC_assigned$Level <- rownames(PC_assigned)
        df_long <- gather(PC_assigned, key = Ab, value = PC, -Level)
        df_long <- separate(df_long, col = Ab, into = c("Type", "Quant", "abQuant", "No_SVs"), sep = "_")
        df_long$Level <- factor(df_long$Level, levels = colnames(taxa), ordered = TRUE)
        df_long$Quant <- factor(df_long$Quant, levels = as.character(abundanceQuantiles), ordered = TRUE)
        
        tr <- ggplot(df_long, aes(x = Quant, y = PC, col = Level))
        tr <- tr + 
                geom_point(size = 2) +
                scale_colour_manual("", values = cbPalette[2:8]) +
                scale_x_discrete(breaks = abundanceQuantiles, labels = paste(abundanceQuantiles, " (", round(abQuantiles), ", ", No_SVs, ")", sep = "")) +
                ylab("percentage assigned") +
                xlab("abundance quantile (abundance, No SVs)") +
                theme_bw(12) +
                theme(axis.text.x = element_text(angle=90, vjust=0.5))
        
        list(PC_assigned, tr)
        
}


#######################################
### FUNCTION: plotTaxLevelvsAbundPrev
#######################################

# Function to determine how many sequences could be assigned to the different taxonomic levels compared to the abundance
# and prevalence of the SV
## Input
# taxa: the output of the assignTaxonomy (dada2) and maybe addSpecies command. nrow = number of sequences and ncol number taxonomic levels.
# seqtab: The abundance table output from Dada2_wrap
## Output
# list of plots ab vs assignment and list of plots prevalence vs assignment

plotTaxLevelvsAbundPrev <- function(taxa, seqtab, color = "#E69F00"){
        
        taxa_vs_ab_pref <- as.data.frame(cbind(Abundance = colSums(seqtab), Prevalence = colSums(seqtab != 0), !is.na(taxa)))
        
        TrList_ab <- list()
        for (i in 3:ncol(taxa_vs_ab_pref)){
                Tr <- ggplot(taxa_vs_ab_pref, aes_string(x = "Abundance", y = colnames(taxa_vs_ab_pref)[i]))
                Tr <- Tr + geom_jitter(col = color, size = 1.5, alpha = 0.5) +
                        scale_y_continuous(limits = c(-0.5, 1.5), breaks = c(0,1), labels = c('No', "Yes")) +
                        xlab("abundance of SV") +
                        ylab(paste("assigned to ", colnames(taxa_vs_ab_pref)[i], sep = "")) +
                        theme_bw() 
                TrList_ab[[i-2]] <- Tr
        }
        names(TrList_ab) <- colnames(taxa_vs_ab_pref)[3:ncol(taxa_vs_ab_pref)]
        
        TrList_prev <- list()
        for (i in 3:ncol(taxa_vs_ab_pref)){
                Tr <- ggplot(taxa_vs_ab_pref, aes_string(x = "Prevalence", y = colnames(taxa_vs_ab_pref)[i]))
                Tr <- Tr + geom_jitter(col = color, size = 1.5, alpha = 0.5) +
                        scale_y_continuous(limits = c(-0.55, 1.55), breaks = c(0,1), labels = c('No', "Yes")) +
                        xlab("SV is present in No of samples") +
                        ylab(paste("assigned to ", colnames(taxa_vs_ab_pref)[i], sep = "")) +
                        theme_bw() 
                TrList_prev[[i-2]] <- Tr
        }
        names(TrList_prev) <- colnames(taxa_vs_ab_pref)[3:ncol(taxa_vs_ab_pref)]
        
        list(TrList_ab = TrList_ab, TrList_prev = TrList_prev)
        
}


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
### rarefaction_curve_own
#######################################
# creates total amplicon steps from 0 to max_total. Then rarefies otu_table(physeq) to these steps and calculates each time
# richness and shannon. It plots these rarefaction curves. If group_var (factor in sample_data(physeq)) is given, it also
# plots grouped (Mean plus SE) and calculates the pairwise t.test p.values between the groups.
# type can be vegan or sample, rarefaction is then either based on vegan::rrarefy or rarefy_sample (sample())
# Outputs a named list with all results.
# NB: depending on step_size and size of otu_table(physeq) expect it to take 2-4 minutes!
# NB2: runs only one rarefaction per step, so depends on seed!

rarefaction_curve_own <- function(physeq, group_var = NULL, max_total = NULL, step_size = 200, type = "vegan", seed = 123) {
        
        
        if (taxa_are_rows(physeq)) {
                seqtab <- t(as(otu_table(physeq), "matrix"))
        } else {
                seqtab <- as(otu_table(physeq), "matrix")
        }
        
        Group <- sample_data(physeq)[[group_var]]
        
        if (!is.null(Group) && !is.factor(Group)) {
                Group <- as.factor(Group)
        }
        
        if (is.null(max_total)) {
                max_total <- quantile(rowSums(seqtab), probs = .25)
        }
        
        steps <- seq(from = 0, to = max_total, by = step_size)
        NoSamples <- nsamples(ps)
        NoSteps <- length(steps)
        richness_matrix <- matrix(data = NA, nrow = NoSamples, ncol = NoSteps)
        shannon_matrix <- matrix(data = NA, nrow = NoSamples, ncol = NoSteps)
        totalamplicons <- rowSums(seqtab)
        
        set.seed(seed)
        # ptm <- Sys.time()
        if (type == "vegan") {
                for (i in 1:length(steps)) {
                        rare_CT <- rrarefy(seqtab, sample = steps[i])
                        richness_matrix[,i] <- rowSums(rare_CT != 0)
                        shannon_matrix[,i] <- vegan::diversity(rare_CT, index = "shannon")
                }
        } else if (type == "sample") {
                #ptm <- Sys.time()
                for (i in 1:length(steps)) {
                        rare_CT <- t(apply(seqtab, 1, function(cnts){rarefy_sample(cnts, size = steps[i])}))
                        richness_matrix[,i] <- rowSums(rare_CT != 0)
                        shannon_matrix[,i] <- vegan::diversity(rare_CT, index = "shannon")
                }
                # Sys.time() - ptm
        } else {
                stop("type neither vegan nor sample")
        }
        
        # Sys.time() - ptm
        
        # set alpha-diversities for samples at which steps > totalamplicons to NA
        for (i in 1:NoSamples) {
                NaIndex <- which(totalamplicons[i] < steps)[1]
                if (!is.na(NaIndex)){
                        richness_matrix[i, NaIndex:ncol(richness_matrix)] <- NA
                        shannon_matrix[i, NaIndex:ncol(shannon_matrix)] <- NA
                }
        }
        
        rownames(richness_matrix) <- rownames(seqtab)
        colnames(richness_matrix) <- paste("step_", steps, sep = "")
        rownames(shannon_matrix) <- rownames(seqtab)
        colnames(shannon_matrix) <- paste("step_", steps, sep = "")
        
        richness_df <- as.data.frame(richness_matrix)
        shannon_df <- as.data.frame(shannon_matrix)
        
        plot_div_df <- function (div_df, type = "alpha diversity") {
                
                div_df$Sample <- rownames(div_df)
                div_df <- tidyr::gather(div_df, key = step, value = div, -Sample)
                div_df <- tidyr::separate(div_df, col = step, into = c("term", "step"), sep = "_")
                div_df$step <- as.numeric(div_df$step)
                
                Tr <- ggplot(div_df, aes(x = step, y = div, group = Sample))
                Tr <- Tr +
                        geom_line() +
                        xlab("total amplicons") +
                        ylab(type) +
                        theme_bw(12)
        }
        
        Tr_richness <- plot_div_df(richness_df, type = "richness")
        Tr_shannon <- plot_div_df(shannon_df, type = "shannon")
        
        
        if (!is.null(Group)) {
                
                richness_df$Group <- Group
                shannon_df$Group <- Group
                
                pairwise.tt_richness <- lapply(richness_df[, -ncol(richness_df)], function(step){ 
                        ptt <- pairwise.t.test(x = step, g = richness_df$Group, alternative = "two", p.adjust.method = "none", var.equal = F, pool.sd = F)
                        ptt$p.value
                })
                
                pairwise.tt_shannon <- lapply(shannon_df[, -ncol(shannon_df)], function(step){ 
                        ptt <- pairwise.t.test(x = step, g = shannon_df$Group, alternative = "two", p.adjust.method = "none", var.equal = F, pool.sd = F)
                        ptt$p.value
                })
                
                
                plot_div_df_group <- function (div_df, type = "alpha diversity") {
                        div_df <- tidyr::gather(div_df, key = step, value = div, -Group)
                        div_df <- tidyr::separate(div_df, col = step, into = c("term", "step"), sep = "_")
                        div_df$step <- as.numeric(div_df$step)
                        div_df <- group_by(div_df, Group, step)
                        div_df <- dplyr::summarise(div_df, Mean_div = mean(div, na.rm = T), SD_div = sd(div, na.rm = T), n = n(), SE_div = SD_div/sqrt(n))
                        
                        Tr <- ggplot(div_df, aes(x = step, y = Mean_div, col = Group))
                        Tr <- Tr + 
                                geom_line() +
                                geom_point(size = 1) +
                                geom_errorbar(aes(ymin = Mean_div-SE_div, ymax = Mean_div+SE_div)) +
                                scale_color_manual("", values = cbPalette[2:8]) +
                                xlab("total amplicons") +
                                ylab(paste(type, "(group mean +- SE)")) +
                                theme_bw(12)
                }
                
                Tr_richness_group <- plot_div_df_group(richness_df, type = "richness")
                Tr_shannon_group <- plot_div_df_group(shannon_df, type = "shannon")
                
        }
        
        outlist <- list(richness_df = richness_df, shannon_df = shannon_df,
                        Tr_richness = Tr_richness, Tr_shannon = Tr_shannon)
        
        if (!is.null(Group)) {
                outlist[[5]] <- pairwise.tt_richness
                outlist[[6]] <- pairwise.tt_shannon
                outlist[[7]] <- Tr_richness_group
                outlist[[8]] <- Tr_shannon_group
        }
        
        return(outlist)
        
}


#######################################
### rarefaction_curve_own_fast
#######################################
# read rarefaction_curve_own, this fast version uses the vegan::rarefy function, so it only generates rarefaction
# for richness. NB: rarefaction_curve_own only does one rarefaction for each step, so in principle you had to run it
# several times and average (since I do not it gives more bumpy curves). vegan::rarefy used here does averaging (read help)
# and is super fast. I still do not use the SE provided by rarefy, so that could be added in here

# # NB: in case you wanted to add sample labels to a Tr_richness_col
# df <- Tr_richness_col$data
# df <- filter(df, step == max(df$step))
# df$step <- df$step + 400
# Tr_richness_col <- Tr_richness_col + 
#         geom_label(data = df, aes(label = Sample)) +
#         scale_x_continuous(limits = c(-100, max(df$step) + 2500))

rarefaction_curve_own_fast <- function(physeq, group_var = NULL, max_total = NULL, step_size = 200, seed = 123) {
        
        if (taxa_are_rows(physeq)) {
                seqtab <- t(as(otu_table(physeq), "matrix"))
        } else {
                seqtab <- as(otu_table(physeq), "matrix")
        }
        
        Group <- sample_data(physeq)[[group_var]]
        
        if (!is.null(Group) && !is.factor(Group)) {
                Group <- as.factor(Group)
        }
        
        if (is.null(max_total)) {
                max_total <- quantile(sample_sums(physeq), probs = .25)
        }
        
        steps <- seq(from = 0, to = max_total, by = step_size)
        NoSamples <- nsamples(ps)
        NoSteps <- length(steps)
        richness_matrix <- matrix(data = NA, nrow = NoSamples, ncol = NoSteps)
        totalamplicons <- rowSums(seqtab)
        
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
        
        plot_div_df <- function (div_df, type = "alpha diversity") {
                
                div_df$Sample <- rownames(div_df)
                div_df <- tidyr::gather(div_df, key = step, value = div, -Sample)
                div_df <- tidyr::separate(div_df, col = step, into = c("term", "step"), sep = "_")
                div_df$step <- as.numeric(div_df$step)
                
                Tr <- ggplot(div_df, aes(x = step, y = div, group = Sample))
                Tr <- Tr +
                        geom_line() +
                        xlab("total amplicons") +
                        ylab(type) +
                        theme_bw(12)
        }
        
        Tr_richness <- plot_div_df(richness_df, type = "richness")
        
        
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
                                xlab("total amplicons") +
                                scale_color_manual("", values = cbPalette[2:8]) +
                                ylab(type) +
                                theme_bw(12)
                }
                
                plot_div_df3 <- function (div_df, type = "richness") {
                        
                        div_df$Sample <- rownames(div_df)
                        div_df$Total <- totalamplicons
                        div_df <- tidyr::gather(div_df, key = step, value = div, -Sample, -Group, -Total)
                        div_df <- tidyr::separate(div_df, col = step, into = c("term", "step"), sep = "_")
                        div_df$step <- as.numeric(div_df$step)
                        
                        Tr <- ggplot(div_df, aes(x = step, y = div, group = Sample, col = Total))
                        Tr <- Tr +
                                geom_line() +
                                scale_color_gradient2("total amp.", low = cbPalette[6], mid = cbPalette[1], high = cbPalette[2], midpoint = 
                                                              median(div_df$Total)) +
                                xlab("total amplicons") +
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
                                scale_color_manual("", values = cbPalette[2:8]) +
                                xlab("total amplicons") +
                                ylab(paste(type, "(group mean +- SE)")) +
                                theme_bw(12)
                }
                Tr_richness_col <- plot_div_df2(richness_df, type = "richness")
                
                Tr_richness_grad <- plot_div_df3(richness_df, type = "richness")
                
                Tr_richness_group <- plot_div_df_group(richness_df, type = "richness")
                
        }
        
        outlist <- list(rarefaction_richness_df = richness_df, Tr_richness = Tr_richness)
        
        if (!is.null(Group)) {
                outlist[[3]] <- pairwise.tt_richness
                outlist[[4]] <- Tr_richness_col
                outlist[[5]] <- Tr_richness_grad
                outlist[[6]] <- Tr_richness_group
                names(outlist)[3:6] <- c("pairwise.tt_richness", "Tr_richness_col", "Tr_richness_grad", "Tr_richness_group")
        }
        
        return(outlist)
        
}


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
### calculate_alphadiversity
#######################################
# NB: I use the estimate_richness function here (that is also used in plot_richness from phyloseq).
# You could easily calculate Shannon, Chao1, Observed self, see alphaDiversityMeasures.Rmd
# estimate_richness itself uses functions from the vegan package

calculate_alphadiversity <- function(physeq, measures = c("Observed", "Shannon")) {
        DF_alpha <- suppressWarnings(phyloseq::estimate_richness(physeq, measures = measures))
        rownames(DF_alpha) <- sample_names(physeq)
        if ("Observed" %in% colnames(DF_alpha)){
                colnames(DF_alpha)[colnames(DF_alpha) == "Observed"] <- "Richness" 
                measures[measures == "Observed"] <- "Richness"
        }
        DF_alpha$Total <- sample_sums(physeq)
        ncol_df_alpha <- ncol(DF_alpha)
        # add residuals of lm fits
        fitlist <- list()
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
                phynew <- phyloseq(otu_table(Mat, taxa_are_rows = TRUE), sample_data(physeq), tax_table(physeq))
        } else {
                if (!identical(rownames(otu_table(physeq)), names(SFs))) {stop("names SFs do not fit to physeq")}
                Mat <- sweep(otu_table(physeq),1,SFs, "/") 
                phynew <- phyloseq(otu_table(Mat, taxa_are_rows = FALSE), sample_data(physeq), tax_table(physeq))
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
                
                phynew
                
        }
}


#######################################
### FUNCTION: KeepAmplicons
#######################################

# Function to find the amplicons that are present in at least "Percentage" of samples
## Input
# taxa: the output of the assignTaxonomy (dada2) and maybe addSpecies command. nrow = number of sequences and ncol number taxonomic levels.
# seqtab: The abundance table output from Dada2_wrap
# percentage: only amplicons present in at least this percentage of samples will be kept
## Output
# taxseq: list with the shortened taxa and seqtab

KeepAmplicons <- function(taxa, seqtab, Percentage = 10){
        
        if(!identical(colnames(seqtab), rownames(taxa))){
                stop("taxa and seqtab do not fit togetehr")
        }
        
        
        PerCSampleValue <- ceiling((Percentage/100)*dim(seqtab)[1])
        KeepAmplis <- colnames(seqtab)[colSums(seqtab != 0) >= PerCSampleValue]
        seqtab.keep <- seqtab[,KeepAmplis]
        taxa.keep <- taxa[KeepAmplis,]
        taxseq <- list(taxa.keep, seqtab.keep)
        return(taxseq)
}






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

# Explanation of original gm function:
# Note zero.propagate just makes sure you get a 0 if there is any zero in your vector
# if there is no zero exp(mean(log(x), na.rm = TRUE)) == exp(sum(log(x[x > 0]), na.rm=TRUE) / length(x))
# I think this zero.propagate is silly, I therefore decided to make my own function, gm_own below

gm = function(x, na.rm=TRUE, zero.propagate = FALSE){
        if(any(x < 0, na.rm = TRUE)){
                return(NaN)
        }
        if(zero.propagate){
                if(any(x == 0, na.rm = TRUE)){
                        return(0)
                }
                exp(mean(log(x), na.rm = na.rm))
        } else {
                exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
        }
}



#######################################
### filttaxa_by_prevalence
#######################################
# filters taxa by prevalence and provides plots that illustrate how the sample_sums
# were changed by the filtering

## Input:
# physeq
# prevalence: in %, only take stay that are present in more than prevalence % of the samples unless MaxCountCheck is used
# MaxCountCheck, MaxCount: if max count is used a taxa can remain even if not fulfilling the prevalence
# criterion if its max count in one of the samples is bigger than MaxCount 
# the idea is that some samples might depend very much on a rare taxa

## Output:
# list with the filtered physeq, the sparsity measures, and 3 different plots (3-5)


filttaxa_by_prevalence <- function(physeq, prevalence = 30, MaxCountCheck = FALSE, MaxCount = 0.2*(min(sample_sums(physeq)))) {
        
        if (MaxCountCheck) {
                physeq_f <- filter_taxa(physeq, function(x){(sum(x != 0) > (prevalence/100)*length(x)) || (max(x) > MaxCount)}, prune = TRUE)
        } else {
                physeq_f <- filter_taxa(physeq, function(x){(sum(x != 0) > (prevalence/100)*length(x))}, prune = TRUE)
        }
        Sparsity <- 100*(sum(otu_table(physeq) == 0)/(ntaxa(physeq)*nsamples(physeq)))
        Sparsity_f <- 100*(sum(otu_table(physeq_f) == 0)/(ntaxa(physeq_f)*nsamples(physeq_f)))
        PCCountsRemoved <- 100*(1-(sample_sums(physeq_f)/sample_sums(physeq)))
        DF <- data.frame(Sample = names(PCCountsRemoved), PCRemoved = PCCountsRemoved)
        DF$Sample <- as.factor(DF$Sample)
        DF <- dplyr::arrange(DF, PCRemoved)
        # relevel so samples are shown from min NoReads to max NoReads
        LevelsWant <- as.character(DF$Sample) 
        for (i in 1:length(LevelsWant)) {
                DF$Sample <- relevel(DF$Sample, ref = LevelsWant[i])
        }
        
        DF2 <- data.frame(before = sample_sums(physeq), after = sample_sums(physeq_f))
        DF2$Sample <- rownames(DF2)
        DF2$Sample <- as.factor(DF2$Sample)
        DF2 <- dplyr::arrange(DF2, before)
        LevelsWant <- as.character(DF2$Sample) 
        for (i in 1:length(LevelsWant)) {
                DF2$Sample <- relevel(DF2$Sample, ref = LevelsWant[i])
        }
        DF2 <- tidyr::gather(DF2, key = physeq, value = count, -Sample)
        
        Tr <- ggplot(DF, aes(x = Sample, y = PCRemoved))
        Tr <- Tr + geom_point(col = "#E69F00", size = 2) +
                ylab("% of amplicon reads removed") +
                xlab("") +
                theme_bw() +
                theme(panel.grid.minor = element_blank(),
                      panel.grid.major.y = element_blank(),
                      panel.grid.major.x = element_line(color = "#999999", size = .15),
                      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6))
        
        Tr2 <- ggplot(data = DF, aes(x = PCRemoved))
        Tr2 <- Tr2 + geom_histogram(binwidth = 2, col = "black", fill = "#E69F00") +
                geom_rug() +
                ylab("No Samples") + 
                xlab("% of amplicon reads removed") +
                theme_bw() + 
                theme(panel.grid.minor = element_blank(),
                      panel.grid.major.y = element_blank(),
                      panel.grid.major.x = element_line(color = "#999999", size = .15))
        
        cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
        
        Tr3 <- ggplot(data = DF2, aes(x = Sample, y = count, color = physeq))
        Tr3 <- Tr3 + geom_point(size = 2) +
                ylab("sample_sums") +
                xlab("") +
                scale_color_manual(values = cbPalette[c(4,2)]) +
                theme_bw() +
                theme(panel.grid.minor = element_blank(),
                      panel.grid.major.y = element_blank(),
                      panel.grid.major.x = element_line(color = "#999999", size = .15),
                      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
                      legend.title = element_blank())
        
        
        list(FilteredPhyseq = physeq_f, Sparsity = data.frame(Before = Sparsity, After = Sparsity_f), PlotReadsRemoved = Tr,
             HistReadsRemoved = Tr2, PlotSampleSumsBeforeAfter = Tr3)
        
}











# ----------------------------- obsolete functions ----------------------------------------------------

# #######################################
# ### adjust_LS
# #######################################
# # own implementation of DESeq2 library size adjustment with little tweaks
# 
# ## Input:
# # physeq
# # zeros.count: zeros.count of gm_own, if TRUE 0 will be considered when calculating the geometric mean
# # if FALSE not and thus the geometric means will be bigger (see gm_own)
# # percentile: the percentile to select the size factor SF of a sample based on its count ratios to the reference sample. In
# # DESeq percentile = 50, i.e. the median is used. NOTE when ignore.zero.ratios = FALSE, you might get SF = 0 in case that more 
# # than 50% of the taxa of in a sample are 0. Function will through a warning then. s
# # ignore.zero.ratios: if TRUE, the SF of each sample (i.e. the percentile) will be calculated based on the pool of non-zero ratios of the sample
# # plots: if TRUE SFs and plots will be given in an output list, otherwise just the adjusted physeq is returned
# 
# 
# adjust_LS <- function(physeq, zeros.count = FALSE, percentile = 50, ignore.zero.ratios = TRUE, plots = FALSE)  {
#         
#         # ---- Step 1: Use the geometric means for each taxa over all samples to calculate a "representative" Reference Sample ------
#         if(taxa_are_rows(physeq)){
#                 GM <- apply(otu_table(physeq), 1, gm_own, zeros.count = zeros.count)   
#         } else {
#                 GM <- apply(otu_table(physeq), 2, gm_own, zeros.count = zeros.count) 
#         }
#         
#         # ---- Step 2: Calculate Size factors --------
#         # -- 2a: Calculate Ratio Matrix by dividing counts of each sample by Reference Sample
#         
#         if(taxa_are_rows(physeq)){
#                 RatioMatrix <- sweep(otu_table(physeq), 1, GM, "/")
#         } else {
#                 RatioMatrix <- sweep(otu_table(physeq), 2, GM, "/") 
#         }
#         
#         # -- 2b: Insert: Calculate 0 Percentage for each sample
#         
#         if(taxa_are_rows(physeq)){
#                 zeroPC <- 100*(colSums(otu_table(physeq) == 0)/ntaxa(physeq))
#         } else {
#                 zeroPC <- 100*(rowSums(otu_table(physeq) == 0)/ntaxa(physeq))
#         }
#         
#         
#         # -- 2c: calculate SFs for each sample depending on given percentile and whether to ignore.zero.ratios
#         # if ignore.zero.ratios = TRUE, the percentile refers only to the non zero ratios in each sample
#         # NB: DESEQ uses ignore.zero.ratios = TRUE
#         if(ignore.zero.ratios){
#                 if (taxa_are_rows(physeq)) {
#                         SFs <- apply(RatioMatrix, 2, function(x){x <- x[x>0]; quantile(x, probs = percentile/100, na.rm = T)})
#                 } else {
#                         SFs <- apply(RatioMatrix, 1, function(x){x <- x[x>0]; quantile(x, probs = percentile/100, na.rm = T)})
#                 }  
#         } else {
#                 if (taxa_are_rows(physeq)) {
#                         SFs <- apply(RatioMatrix, 2, quantile, probs = percentile/100, na.rm = T)
#                 } else {
#                         SFs <- apply(RatioMatrix, 1, quantile, probs = percentile/100, na.rm = T) 
#                 }   
#         }
#         
#         if(min(SFs) == 0) { warning("in at least one sample the Size Factor was 0!") }
#         
#         
#         # --- 3: calculate the new counts and put into a physeq object
#         
#         if(taxa_are_rows(physeq)){
#                 if (!identical(colnames(otu_table(physeq)), names(SFs))) {stop("names SFs do not fit to physeq")}
#                 Mat <- sweep(otu_table(physeq),2,SFs, "/")
#                 phynew <- phyloseq(otu_table(Mat, taxa_are_rows = TRUE), sample_data(physeq), tax_table(physeq))
#         } else {
#                 if (!identical(rownames(otu_table(physeq)), names(SFs))) {stop("names SFs do not fit to physeq")}
#                 Mat <- sweep(otu_table(physeq),1,SFs, "/") 
#                 phynew <- phyloseq(otu_table(Mat, taxa_are_rows = FALSE), sample_data(physeq), tax_table(physeq))
#         }
#         
#         if (plots){
#                 
#                 # ---- step 4: Generate Plots that illustrate the process
#                 # -- 4a: calculate for each sample a histogram of the ratios used to determine the SF
#                 
#                 cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#                 
#                 # Define the histos function
#                 histos <- function(x, SF, zPC, SampleName) {
#                         x <- data.frame(x = x)
#                         Tr <- ggplot(x, aes(x = x))
#                         Tr <- Tr + geom_histogram(binwidth = 0.3, col = "black", fill = "#E69F00") +
#                                 geom_rug() +
#                                 geom_vline(xintercept = SF, col = "#009E73") +
#                                 ylab("Frequency") + 
#                                 xlab("Count/(Count of Ref Sample)") +
#                                 ggtitle(paste("Sample: ", SampleName, "; zeroPC: ", round(zPC, 2), "; size factor: ", round(SF,2), sep = "")) +
#                                 theme_bw() + 
#                                 theme(panel.grid.minor = element_blank(),
#                                       panel.grid.major.y = element_blank(),
#                                       panel.grid.major.x = element_line(color = "#999999", size = .15))
#                 }
#                 
#                 
#                 if(taxa_are_rows(physeq) && (!all.equal(colnames(RatioMatrix), names(SFs), names(zeroPC)))) {warning("names messed up!")}
#                 if(!taxa_are_rows(physeq) && (!all.equal(rownames(RatioMatrix), names(SFs), names(zeroPC)))) {warning("names messed up!")}
#                 
#                 if(ignore.zero.ratios){
#                         if(taxa_are_rows(physeq)){
#                                 HistoList <- lapply(1:nsamples(physeq), function(i){histos(x = RatioMatrix[,i][RatioMatrix[,i] > 0], SF = SFs[i], zPC = zeroPC[i], SampleName = names(SFs)[i])})
#                         } else {
#                                 HistoList <- lapply(1:nsamples(physeq), function(i){histos(x = RatioMatrix[i,][RatioMatrix[i,] > 0], SF = SFs[i], zPC = zeroPC[i], SampleName = names(SFs)[i])})
#                         }  
#                 } else {
#                         if(taxa_are_rows(physeq)){
#                                 HistoList <- lapply(1:nsamples(physeq), function(i){histos(x = RatioMatrix[,i], SF = SFs[i], zPC = zeroPC[i], SampleName = names(SFs)[i])})
#                         } else {
#                                 HistoList <- lapply(1:nsamples(physeq), function(i){histos(x = RatioMatrix[i,], SF = SFs[i], zPC = zeroPC[i], SampleName = names(SFs)[i])}) 
#                         }   
#                 }
#                 
#                 # -- 3b: compare calculated SFs to library sizes of the samples
#                 if(!identical(names(SFs), names(sample_sums(physeq)))){warning("Probably some Mix Up CHECK!")}
#                 comp <- data.frame(Sample = names(SFs), TotAmps = sample_sums(physeq), SFs = SFs)
#                 comp$SFsNormed <- comp$SFs/median(comp$SFs)
#                 comp$TotAmpsNormed <- comp$TotAmps/mean(comp$TotAmps)
#                 comp <- dplyr::arrange(comp, desc(TotAmpsNormed))
#                 comp$Sample <- as.factor(comp$Sample)
#                 LevelsWant <- as.character(comp$Sample)
#                 for(i in 1:length(LevelsWant)){
#                         comp$Sample <- relevel(comp$Sample, LevelsWant[i])
#                 }
#                 
#                 comp <- comp[c(1,4,5)]
#                 names(comp)[c(2,3)] <- c("SizeFactor_DESeq", "TotalAmplicons_relAb")
#                 comp <- tidyr::gather(comp, key = Corrector, value = NormedValue, -Sample)
#                 Tr <- ggplot(comp, aes(x = Sample, y = NormedValue, color = Corrector))
#                 Tr <- Tr + geom_point(size = 2) +
#                         xlab("") +
#                         ylab("correction value (normalized to median)") +
#                         ggtitle(paste("Median SF: ", round(median(SFs),3), " Median TotAmp: ", round(median(sample_sums(physeq)),3), sep = "")) +
#                         scale_color_manual(values = cbPalette[c(4,2)]) +
#                         theme_bw() +
#                         theme(panel.grid.minor = element_blank(),
#                               panel.grid.major.y = element_blank(),
#                               panel.grid.major.x = element_line(color = "#999999", size = .15),
#                               axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
#                               legend.title = element_blank())
#                 
#                 # -- 3c: save Histograms of the SFs and sample_sums
#                 histo2 <- function(x, xtitle, gtitle) {
#                         x <- data.frame(x = x)
#                         Tr <- ggplot(x, aes(x = x))
#                         Tr <- Tr + geom_histogram(binwidth = diff(range(x))/60, col = "black", fill = "#E69F00") +
#                                 geom_rug() +
#                                 geom_vline(xintercept = median(x$x), col = "#009E73", size = 1) +
#                                 ylab("No Samples") + 
#                                 xlab(xtitle) +
#                                 ggtitle(gtitle) +
#                                 theme_bw() + 
#                                 theme(panel.grid.minor = element_blank(),
#                                       panel.grid.major.y = element_blank(),
#                                       panel.grid.major.x = element_line(color = "#999999", size = .15))
#                 }
#                 
#                 Tr2 <- histo2(SFs, xtitle = "Size Factors", gtitle = "Size Factors a la DESeq")
#                 Tr3 <- histo2(sample_sums(physeq), xtitle = "Total Amplicons", gtitle = "Size Factors a la relative abundance")
#                 
#                 List <- list(Physeq = phynew, SFs = SFs, SizeFactorCompare = Tr, SFHistos = list(Tr2, Tr3), RefSample = GM, RatioMatrix = RatioMatrix, zeroPC = zeroPC, HistoList = HistoList)
#                 
#                 
#         } else {
#                 
#                 phynew
#                 
#         }
#         
#         
# }


# #######################################
# ### FUNCTION: plotTaxLevelvsprevalence
# #######################################
# 
# # Function to determine how many sequences could be assigned to the different taxonomic levels compared to the number of smaples the amplicons are present
# ## Input
# # taxa: the output of the assignTaxonomy (dada2) and maybe addSpecies command. nrow = number of sequences and ncol number taxonomic levels.
# # seqtab: The abundance table output from Dada2_wrap
# ## Output
# # TLvsprevalence: list of 2 data frames that state for each prevalence level the assigned taxonomic levels, first data frame as numbers, second as percentage 
# ## requires: 
# # dplyr!
# 
# plotTaxLevelvsprevalence <- function(taxa, seqtab){
#         
#         TLvsprevalence <- as.data.frame(cbind(colSums(seqtab != 0), apply(taxa, 2, function(x){!is.na(x)})))
#         colnames(TLvsprevalence)[1] <- "InNoSamples"
#         TLvsprevalence <- dplyr::group_by(TLvsprevalence, InNoSamples)
#         
#         if("Species" %in% colnames(taxa)){
#                 
#                 TLvsprevalence <- as.data.frame(dplyr::summarise(TLvsprevalence, NoAmplicons = n(), Kingdom = sum(Kingdom), Phylum = sum(Phylum), Class = sum(Class),
#                                                                  Order = sum(Order), Family = sum(Family), Genus = sum(Genus), Species = sum(Species)))
#                 TLPC <- dplyr::mutate(TLvsprevalence, Kingdom = Kingdom/NoAmplicons, Phylum = 100*(Phylum/NoAmplicons), Class = 100*(Class/NoAmplicons), Order =
#                                               100*(Order/NoAmplicons), Family = 100*(Family/NoAmplicons), Genus = 100*(Genus/NoAmplicons), Species = 100*(Species/NoAmplicons))
#                 
#         } else {
#                 
#                 TLvsprevalence <- as.data.frame(dplyr::summarise(TLvsprevalence, NoAmplicons = n(), Kingdom = sum(Kingdom), Phylum = sum(Phylum), Class = sum(Class),
#                                                                  Order = sum(Order), Family = sum(Family), Genus = sum(Genus)))
#                 TLPC <- dplyr::mutate(TLvsprevalence, Kingdom = Kingdom/NoAmplicons, Phylum = 100*(Phylum/NoAmplicons), Class = 100*(Class/NoAmplicons), Order =
#                                               100*(Order/NoAmplicons), Family = 100*(Family/NoAmplicons), Genus = 100*(Genus/NoAmplicons))
#                 
#         }
#         
#         TLvsprevalence <- list(Counts = TLvsprevalence, PerC = TLPC)
#         
#         return(TLvsprevalence)
#         
# }