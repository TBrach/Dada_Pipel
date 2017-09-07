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
        InputSave <- list(seqtab = seqtab, 
                          minBoot = minBoot,
                          allowMultiple = allowMultiple,
                          tryRC = tryRC,
                          PathToRefs = PathToRefs, 
                          RefDataBase = RefDataBase,
                          SpeciesDB = SpeciesDB,
                          PathToSave = PathToSave)
        
        
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
### FUNCTION: construct_phylogenetic_tree
#######################################
## Background: 
# the entire function is basically a copy from: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4955027/ (page 8)
# the function uses the packages DECIPHER (bioconductor) and phangorn (cran)
# I had once a clash of the DECIPHER package with dada2 when dada2 was loaded first, so maybe function needs to be run in a fresh session
# you should run the wrapper on the server because the last phangorn::optim.pml takes long time
## Inputs, 
# savepath = determines were the output data list will be saved
# all other inputs: refer to arguments for phangorn::optim.pml and I have honestly no clue about them at the moment
# so read code and help
## output: 
# a list with the GTR: generalized time-reversible with Gamma rate variation model fit containing the tree, and the Input data saved
# the list is saved as phylogenetic_tree.rds in savepath


construct_phylogenetic_tree <- function(seqtab.nochim, savepath,
                                        k = 4, inv = .2, 
                                        model = "GTR", rearrangement = "stochastic",
                                        trace = 0) {
        
        library(DECIPHER)
        library(phangorn); packageVersion("phangorn")
        
        RVer <- R.Version()
        RVer <- RVer$version.string
        PackageVersions <- data.frame(Package = c("R", "DECIPHER", "phangorn"),
                                      Version = c(RVer,
                                                  as.character(packageVersion("DECIPHER")),
                                                  as.character(packageVersion("phangorn"))))
        
        
        seqs <- colnames(seqtab.nochim) # simple character string
        names(seqs) <- seqs
        
        Inputs <- list(seqs = seqs,
                       savepath = savepath,
                       k = k,
                       inv = inv,
                       model = model,
                       rearrangement = rearrangement,
                       trace = trace)
        
        # ---- multiple alignment of SVs using "DECIPHER" (takes only 10 seconds) ----
        
        alignment <- DECIPHER::AlignSeqs(DNAStringSet(seqs), anchor = NA)
        
        # --------
        
        # ---- use alignment to construct phylogenetic tree (phangorn package) -----
        
        # -- First construct a neighbor-joining tree --
        
        phang.align <- phangorn::phyDat(as(alignment, "matrix"), type = "DNA") # just transformation of DNA data into phyDat format
        dm <- phangorn::dist.ml(phang.align) # computes pairwise distances
        
        treeNJ <- phangorn::NJ(dm)
        
        # NB: They warn: tip order != sequence order
        # however, I got a TRUE here
        # all.equal(treeNJ$tip.label, seqs, check.attributes = FALSE)
        
        fit <- phangorn::pml(treeNJ, data = phang.align) # warns: negative edges length changed to 0!
        
        # ----
        
        # -- use the neighbor-joining tree as a start point to fit a GTR+G+I (generalized time-reversible with Gamma rate variation) maximum 
        # likelihood tree --
        
        fitGTR <- update(fit, k = k, inv = inv)
        message("starting the time consuming fit now")
        fitGTR <- phangorn::optim.pml(fitGTR, model = model, optInv = TRUE, optGamma = TRUE,
                                      rearrangement = rearrangement, control = pml.control(trace = trace))
        message("done with the time consuming fit:)")
        
        
        out <- list(fitGTR = fitGTR,
                    Inputs = Inputs,
                    PackageVersions = PackageVersions)
        
        saveRDS(out, file = file.path(savepath, "phylog_tree.rds"))
        # ----
        
        # --------
        
        return(out)
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
        
        if (!is.null(group_var)){
                Group <- sample_data(physeq)[[group_var]] 
                if (!is.null(Group) && !is.factor(Group)) {
                        Group <- as.factor(Group)
                }
        } else { Group <- NULL}
        
        
        if (is.null(max_total)) {
                max_total <- quantile(sample_sums(physeq), probs = .25)
        }
        
        steps <- seq(from = 0, to = max_total, by = step_size)
        NoSamples <- nsamples(physeq)
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
                        xlab("total counts in sample") +
                        ylab(type) +
                        theme_bw(12)
        }
        
        Tr_richness <- plot_div_df(richness_df, type = "richness")
        
        plot_div_df3 <- function (div_df, type = "richness") {
                
                div_df$Sample <- rownames(div_df)
                div_df$Total <- totalamplicons
                div_df <- tidyr::gather(div_df, key = step, value = div, -Sample, -Total)
                div_df <- tidyr::separate(div_df, col = step, into = c("term", "step"), sep = "_")
                div_df$step <- as.numeric(div_df$step)
                
                Tr <- ggplot(div_df, aes(x = step, y = div, group = Sample, col = Total))
                Tr <- Tr +
                        geom_line() +
                        scale_color_gradient2("total counts bef.", low = cbPalette[6], mid = cbPalette[1], high = cbPalette[2], midpoint = 
                                                      median(div_df$Total)) +
                        xlab("total counts in sample") +
                        ylab(type) +
                        theme_bw(12)
        }
        
        Tr_richness_grad <- plot_div_df3(richness_df, type = "richness")
        
        
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
                                scale_color_manual("", values = cbPalette[2:8]) +
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
                                xlab("total counts in sample") +
                                ylab(paste(type, "(group mean +- SE)")) +
                                theme_bw(12)
                }
                Tr_richness_col <- plot_div_df2(richness_df, type = "richness")
                
                Tr_richness_group <- plot_div_df_group(richness_df, type = "richness")
                
        } 
        
        outlist <- list(rarefaction_richness_df = richness_df, Tr_richness = Tr_richness, Tr_richness_grad = Tr_richness_grad)
        
        if (!is.null(Group)) {
                outlist[[4]] <- pairwise.tt_richness
                outlist[[5]] <- Tr_richness_col
                outlist[[6]] <- Tr_richness_group
                names(outlist)[4:6] <- c("pairwise.tt_richness", "Tr_richness_col", "Tr_richness_group")
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
        # add residuals of lm fits to total_amplicons
        fitlist <- list()
        for (i in 1:length(measures)){
                fit <- lm(DF_alpha[,measures[i]] ~ DF_alpha[,"Total"])
                fitlist[[i]] <- fit
                DF_alpha[, ncol_df_alpha + i] <- residuals(fit)
                colnames(DF_alpha)[ncol_df_alpha + i] <- paste(measures[i], "resids", sep = "_")
        }
        
        names(fitlist) <- measures
        
        # if you have, add lm fits to FilteredReads
        fitlist_FilteredReads <- list()
        
        if (!is.null(sample_data(physeq)$FilteredReads)){
                DF_alpha$filtered_reads <- sample_data(physeq)$FilteredReads
                
                
                ncol_df_alpha <- ncol(DF_alpha)
                for (i in 1:length(measures)){
                        fit <- lm(DF_alpha[,measures[i]] ~ DF_alpha[,"filtered_reads"])
                        fitlist_FilteredReads[[i]] <- fit
                        DF_alpha[, ncol_df_alpha + i] <- residuals(fit)
                        colnames(DF_alpha)[ncol_df_alpha + i] <- paste(measures[i], "resids", "FilteredReads", sep = "_")
                }
        }
        
        names(fitlist_FilteredReads) <- measures
        
        DF_alpha <- data.frame(Sample = rownames(DF_alpha), DF_alpha, sample_data(physeq))
        
        outlist <- list(DF_alpha = DF_alpha, fitlist = fitlist, fitlist_FilteredReads = fitlist_FilteredReads)
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
                
                phynew
                
        }
}

#######################################
### FUNCTION: alpha_diversity_wrapper
#######################################
# just wraps around other functions here to save writing time
# currently only works for alpha_div_measures = c("Observed", "Shannon")
alpha_diversity_wrapper <- function(physeq, alpha_div_measures = c("Observed", "Shannon")){
        
        DF_alpha_list <- calculate_alphadiversity(physeq = physeq, measures = alpha_div_measures)
        DF_alpha <- DF_alpha_list[[1]]
        
        # just I prefer Richness over Observed
        alpha_div_measures2 <- alpha_div_measures
        if ("Observed" %in% alpha_div_measures) {
                alpha_div_measures2[alpha_div_measures2 == "Observed"] <- "Richness" 
        }
        
        pairwise.tt_rich <- pairwise.t.test(x = DF_alpha$Richness, g = DF_alpha[[group_var]], alternative = "two", p.adjust.method = "none", var.equal = F, pool.sd = F)
        # compare, to get the same results as individual t.test you need pool.sd = F
        # t.test(x = DF_alpha$Richness[DF_alpha$Group == "Old"], y = DF_alpha$Richness[DF_alpha$Group == "MiddleAged"], alternative = "two", var.equal = F)$p.value
        
        pairwise.tt_shannon <- pairwise.t.test(x = DF_alpha$Shannon, g = DF_alpha[[group_var]], alternative = "two", p.adjust.method = "none", var.equal = F, pool.sd = F)
        
        if(var(DF_alpha$Total) > 0) {
                pairwise.tt_totalcounts <- pairwise.t.test(x = DF_alpha$Total, g = DF_alpha[[group_var]], alternative = "two", p.adjust.method = "none", var.equal = F, pool.sd = F)
        } else {
                pairwise.tt_totalcounts <- NA
        }
        
        pairwise.tt_rich_resids <- pairwise.t.test(x = DF_alpha$Richness_resids, g = DF_alpha[[group_var]], alternative = "two", p.adjust.method = "none", var.equal = F, pool.sd = F)
        
        pairwise.tt_shannon_resids <- pairwise.t.test(x = DF_alpha$Shannon_resids, g = DF_alpha[[group_var]], alternative = "two", p.adjust.method = "none", var.equal = F, pool.sd = F)
        
        if(var(DF_alpha$filtered_reads) > 0) {
                pairwise.tt_filtered_reads <- pairwise.t.test(x = DF_alpha$filtered_reads, g = DF_alpha[[group_var]], alternative = "two", p.adjust.method = "none", var.equal = F, pool.sd = F)
        } else {
                pairwise.tt_filtered_reads <- NA
        }
        
        pairwise.tt_rich_resids_filt <- pairwise.t.test(x = DF_alpha$Richness_resids_FilteredReads, g = DF_alpha[[group_var]], alternative = "two", p.adjust.method = "none", var.equal = F, pool.sd = F)
        
        pairwise.tt_shannon_resids_filt <- pairwise.t.test(x = DF_alpha$Shannon_resids_FilteredReads, g = DF_alpha[[group_var]], alternative = "two", p.adjust.method = "none", var.equal = F, pool.sd = F)
        
        # add the linear models correcting for total or filtered amplicons
        fitter_rich_totalcounts <- lm(DF_alpha[["Richness"]] ~ DF_alpha[["Total"]] + DF_alpha[[group_var]])
        fitter_shannon_totalcounts <- lm(DF_alpha[["Shannon"]] ~ DF_alpha[["Total"]] + DF_alpha[[group_var]])
        fitter_rich_filteredreads <- lm(DF_alpha[["Richness"]] ~ DF_alpha[["filtered_reads"]] + DF_alpha[[group_var]])
        fitter_shannon_filteredreads <- lm(DF_alpha[["Shannon"]] ~ DF_alpha[["filtered_reads"]] + DF_alpha[[group_var]])
        fitter.list <- list(fitter_rich_totalcounts = fitter_rich_totalcounts, fitter_shannon_totalcounts = fitter_shannon_totalcounts,
                            fitter_shannon_filteredreads = fitter_shannon_filteredreads, 
                            fitter_shannon_filteredreads = fitter_shannon_filteredreads)
        
        
        TrListBP <- boxplot_alphaDiv_fromDF(DF = DF_alpha, color = group_var, group = group_var, measures = c(alpha_div_measures2, paste(alpha_div_measures2, "resids", sep = "_"), paste(alpha_div_measures2, "resids", "FilteredReads", sep = "_")))
        
        TrList_lm <- plot_alphaDivVstotalCounts_fromList(DF_List = DF_alpha_list, measures = alpha_div_measures2, color = group_var)
        
        TrList_lm_filteredReads <- plot_alphaDivVsfilteredReads_fromList(DF_List = DF_alpha_list, measures = alpha_div_measures2, color = group_var)
        
        out <- list(DF_alpha_list = DF_alpha_list, TrListBP = TrListBP, TrList_lm = TrList_lm, TrList_lm_filteredReads = TrList_lm_filteredReads, pairwise.tt_rich = pairwise.tt_rich,
                    fitter_list = fitter.list, pairwise.tt_shannon = pairwise.tt_shannon, pairwise.tt_totalcounts = pairwise.tt_totalcounts, pairwise.tt_rich_resids = pairwise.tt_rich_resids,
                    pairwise.tt_shannon_resids = pairwise.tt_shannon_resids, pairwise.tt_filtered_reads = pairwise.tt_filtered_reads,
                    pairwise.tt_rich_resids_filt = pairwise.tt_rich_resids_filt, pairwise.tt_shannon_resids_filt = pairwise.tt_shannon_resids_filt)
        
}

#######################################
### FUNCTION: arrange_p_values
#######################################
# see Generalized_Phyloseq_Analysis to understand it, it is just about putting a lot of p-value matrixes into one wide data.frame

arrange_p_values <- function(pairwise_list){
        pvalue_types <- strsplit(names(pairwise_list), split = "t_")
        pvalue_types <- sapply(pvalue_types, `[`, 2)
        pairwise_list <- lapply(pairwise_list, function(pairw) {
                pmat <- pairw$p.value
                rowCol <- expand.grid(rownames(pmat), colnames(pmat))
                labs <- rowCol[as.vector(lower.tri(pmat,diag=T)),]
                cbind(labs, p_value = pmat[lower.tri(pmat, diag = T)])
        })
        pairwise_list <- lapply(1:length(pairwise_list), function(i){cbind(pairwise_list[[i]], Type = pvalue_types[i])})
        pairwise_long <- do.call("rbind", pairwise_list)
        pairwise_wide <- spread(pairwise_long, key = Type, value = p_value)
}



#######################################
### FUNCTION: calc_distances
#######################################

calc_distances <- function(physeq, dist_methods = c("bray")) {
        
        dist_list <- vector("list", length(dist_methods))
        names(dist_list) = dist_methods
        
        for (i in dist_methods) {
                iDist <- phyloseq::distance(physeq, method=i)
                dist_list[[i]] = iDist
        }
        
        return(dist_list)
        
}

#######################################
### FUNCTION: calc_ordination_from_distances
#######################################

calc_ordination_from_distances <- function(physeq, dist_list, ordination_type = "PCoA", group_var = NULL, coord_cor = FALSE){
        
        ordination_list <- vector("list", length(dist_list))
        DFList <- vector("list", length(dist_list))
        DF_taxa_List <- vector("list", length(dist_list))
        # TrList <- vector("list", length(dist_list))
        TrList_own <- vector("list", length(dist_list))
        TrList_taxa <- vector("list", length(dist_list))
        
        axes <- 1:2 # currently only allowed to plot first and second
        
        for (i in seq_along(dist_list)) {
                
                ordination <- phyloseq::ordinate(physeq, method = ordination_type, distance = dist_list[[i]])
                ordination_list[[i]] <- ordination
                DF <- phyloseq::plot_ordination(physeq, ordination_list[[i]], color = group_var, justDF = TRUE)
                DFList[[i]] <- DF # just the first two axes cbind to sample_data in physeq
                
                x = colnames(DF)[1]
                y = colnames(DF)[2]
                Tr <- ggplot(DF, aes_string(x = x, y = y, col = group_var)) 
                Tr <- Tr + geom_point() +
                        scale_color_manual("", values = cbPalette[2:8]) +
                        theme_bw(12) +
                        ggtitle(names(dist_list)[i])
                
                # for labelling axes
                if (ordination_type == "PCoA" || ordination_type == "NMDS") {
                        if (ordination_type == "PCoA") {
                                eigvec <- phyloseq:::extract_eigenvalue.pcoa(ordination)
                        } else {
                                eigvec <- phyloseq:::extract_eigenvalue.default(ordination)
                        }
                        
                        if (length(eigvec[axes]) > 0){
                                fracvar = eigvec[axes]/sum(eigvec)
                                percvar = round(100 * fracvar, 1)
                                strivar = as(c(Tr$label$x, Tr$label$y), "character")
                                strivar = paste0(strivar, " (", percvar, " %)")
                                Tr <- Tr + xlab(strivar[1]) + ylab(strivar[2]) 
                        }
                        
                        if (!is.null(eigvec) && coord_cor) {
                                Tr <- Tr + coord_fixed(sqrt(eigvec[2] / eigvec[1]))
                        }
                        
                }
                 
                TrList_own[[i]] <- Tr
                rm(Tr)
                        
                
                # TrList[[i]] <- phyloseq::plot_ordination(physeq, ordination_list[[i]], color = group_var) + ggtitle(names(dist_list)[i])
                
                DF_taxa <- phyloseq::plot_ordination(physeq, ordination_list[[i]], type = "taxa", color = "Phylum", justDF = TRUE)
                DF_taxa_List[[i]] <- DF_taxa
                
                x = colnames(DF_taxa)[1]
                y = colnames(DF_taxa)[2]
                Tr <- ggplot(DF_taxa, aes_string(x = x, y = y, col = "Phylum")) 
                Tr <- Tr + geom_point() +
                        # scale_color_manual("", values = cbPalette[2:8]) +
                        theme_bw(12) +
                        ggtitle(names(dist_list)[i])
                
                # for labelling axes
                if (ordination_type == "PCoA" || ordination_type == "NMDS") {
                        if (ordination_type == "PCoA") {
                                eigvec <- phyloseq:::extract_eigenvalue.pcoa(ordination)
                        } else {
                                eigvec <- phyloseq:::extract_eigenvalue.default(ordination)
                        }
                        
                        if (length(eigvec[axes]) > 0){
                                fracvar = eigvec[axes]/sum(eigvec)
                                percvar = round(100 * fracvar, 1)
                                strivar = as(c(Tr$label$x, Tr$label$y), "character")
                                strivar = paste0(strivar, " (", percvar, " %)")
                                Tr <- Tr + xlab(strivar[1]) + ylab(strivar[2]) 
                        }
                        
                        if (!is.null(eigvec) && coord_cor) {
                                Tr <- Tr + coord_fixed(sqrt(eigvec[2] / eigvec[1]))
                        }
                        
                }
                
                TrList_taxa[[i]] <- Tr
                rm(Tr)
                
        }
        
        names(ordination_list) <- names(TrList_taxa) <- names(DFList) <- names(DF_taxa_List) <- names(TrList_own) <- names(dist_list)
        out <- list(ordination_list = ordination_list, DFList = DFList, DF_taxa_List = DF_taxa_List, ordination_Tr_samples = TrList_own, ordination_Tr_taxa = TrList_taxa)
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




############################################
## simulate_totalabVSrichness_rarefaction ##
############################################
# Input:
# physeq is needed to simulate the sample, i.e. the composition of the original sample
# specifically: the sample in physeq with the highest richness is used as base_sample (only the SVs > 0)
# no_low_extra: defines the Number of low abundance SVs that will be added to the base_sample, the "counts"
# of these extra SVs are random numbers between 0 and min_count in base_sample
# NB: the low_extra SVs make sure that the generated samples have more low abundance counts than the base_sample, and
# therefore plateau slower in rarefaction curves
# NB2: you could of course think about other ways to generate the original proportion
# the idea of the simulation is then: 1 million amplicons of the original proportion go to the sequencer (DNA)
# S1 is size1 amplicons from this original
# S2 is size2 amplicons from this original
# S3 is size2 amplicons rarefied from S1
# richnesses are compared and rarefaction curves are generated, all is saved in the output list

# further Inputs:
# type: "sample" or "vegan" decides whether the rarefaction steps are based on sample() or vegan::rrarefy command
# sd: is used to spread the low_extra_SVs "counts" in an rnorm command, the bigger the higher the spread, but note anyway bewtween 0 and min_count
# no_DNA_total_amplicons: "the number of amplicons generated and put in the seqeuncer", the higher the more secure that the richness of the
# DNA sample was 100%


simulate_totalabVSrichness_rarefaction <- function(physeq, size1 = 50000, size2 = 15000, nsims = 100, no_low_extra = 92, type = "sample", seed = 1576,
                                                   sd = 2, no_DNA_total_amplicons = 1e6){
        
        set.seed(seed)
        seqtab <- as(otu_table(physeq), "matrix")
        base_sample <- seqtab[which.max(rowSums(seqtab != 0)),]
        min_count <- min(base_sample[base_sample > 0])
        total_SVs <- no_low_extra + sum(base_sample > 0)
        
        DNA_richness <- vector("numeric", length = nsims)
        S1_richness <- vector("numeric", length = nsims)
        S2_richness <- vector("numeric", length = nsims)
        S3_richness <- vector("numeric", length = nsims)
        S1s <- matrix(nrow = nsims, ncol = total_SVs)
        
        for (i in 1:nsims) {
                
                factors <- abs(rnorm(n = no_low_extra, sd = 2))
                factors <- factors/max(factors)
                extra_low_SVs <- factors*min_count
                
                sample_prop <- c(extra_low_SVs, base_sample[base_sample > 0])
                sample_prop <- sample_prop/sum(sample_prop)
                
                DNA_sample <- round(sample_prop * no_DNA_total_amplicons)
                DNA_richness[i] <- sum(DNA_sample > 0)
                
                if (type == "vegan") {
                        S1 <- rrarefy(matrix(DNA_sample, nrow = 1), sample = size1)
                        S2 <- rrarefy(matrix(DNA_sample, nrow = 1), sample = size2)
                        S3 <- rrarefy(S1, sample = size2)
                        S1 <- as.vector(S1)
                        S2 <- as.vector(S2)
                        S3 <- as.vector(S3)
                        S1_richness[i] <- sum(S1 > 0)
                        S2_richness[i] <- sum(S2 > 0)
                        S3_richness[i] <- sum(S3 > 0)
                        S1s[i,] <- S1
                        
                } else if (type == "sample") {
                        
                        S1 <- rarefy_sample(DNA_sample, size = size1)
                        S2 <- rarefy_sample(DNA_sample, size = size2)
                        S3 <- rarefy_sample(S1, size = size2)
                        S1_richness[i] <- sum(S1 > 0)
                        S2_richness[i] <- sum(S2 > 0)
                        S3_richness[i] <- sum(S3 > 0)
                        S1s[i,] <- S1
                        
                } else {
                        stop("wrong type")
                }
                
        }
        
        DF <- data.frame(DNA = DNA_richness, S1 = S1_richness, S2 = S2_richness, S3 = S3_richness)
        
        DFl <- gather(DF, key = sample, value = richness)
        DFl$Sample <- factor(DFl$sample, levels <- c("DNA", "S1", "S2", "S3"), ordered = TRUE)
        
        pair_tt <- pairwise.t.test(x = DFl$richness, g = DFl$Sample, alternative = "two", p.adjust.method = "none", var.equal = F, pool.sd = F)
        
        Tr <- ggplot(DFl, aes(x = sample, y = richness, col = sample))
        
        Tr <- Tr +
                geom_boxplot() +
                geom_jitter(alpha = .7) +
                scale_color_manual("", values = cbPalette[2:8]) +
                theme_bw(12)
        
        
        # add rarefaction curves for 5 of the S1s plus the base_sampled supplied with 0s for the low_extra
        CT <- rbind(S1s[sample(nsims, 5),], c(rep(0, no_low_extra), base_sample[base_sample > 0]))
        samdf <- sample_data(ps)
        samdf <- samdf[1:nrow(CT),]
        rownames(CT) <- rownames(samdf)
        colnames(CT) <- paste("T_", 1:ncol(CT), sep = "")
        pss <- phyloseq(otu_table(CT, taxa_are_rows = FALSE), 
                        sample_data(samdf))
        
        Curves <- rarefaction_curve_own_fast(physeq = pss, group_var = NULL)
        
        # show distribution of SVs with abundance below 20 and above 0
        low_abund_SVs <- apply(CT, 1, function(x){x[x > 0 & x < 20]})
        df_plot <- data.frame(Sample = rep(names(low_abund_SVs), sapply(low_abund_SVs, length)),
                              Abund = unlist(low_abund_SVs))
        Trr <- ggplot(df_plot, aes(x = Abund, group = Sample, col = Sample))
        Trr <- Trr + geom_density() +
                scale_color_manual("", values = cbPalette[2:8], labels = c(paste("Sim", 1:5), "dada2_sam")) +
                xlab("Abundances below 20") +
                theme_bw(20)
        
        out <- list(S1s = S1s, DF_richness = DF, pair_tt_of_richness = pair_tt, Tr = Tr, CT = CT, Curves = Curves, Trr = Trr, sample_prop = sample_prop, base_sample = base_sample)
        
}



#######################################
### FUNCTION: distance_t_analyse
#######################################
# INPUT:
# dist_list: named list of dist objects (so for each distance one object)
# physeq: phyloseq object used to make the dist objects
# group_var: character identifying the grouping variable in the sample_data of the phyloseq
# OUTPUT:
# list of two named lists, one with the plots and one with the data_frames of the p_values from pairwise t.tests

distance_t_analyse <- function(dist_list, physeq, group_var) {
        
        TrList <- vector(mode = "list", length = length(dist_list))
        pValList <- vector(mode = "list", length = length(dist_list))
        
        for (i in 1:length(dist_list)){
                DistMat <- as(dist_list[[i]], "matrix")
                rowCol <- expand.grid(rownames(DistMat), colnames(DistMat))
                labs <- rowCol[as.vector(lower.tri(DistMat, diag=F)),]
                df <- cbind(labs, DistMat[lower.tri(DistMat, diag=F)])
                colnames(df) <- c("Row","Col","Distance")
                samdf <- data.frame(Sample = sample_names(physeq), Group = sample_data(physeq)[[group_var]])
                df$Row_Group <- samdf$Group[match(df$Row, samdf$Sample)]
                df$Col_Group <- samdf$Group[match(df$Col, samdf$Sample)]
                
                # add GroupvsGroup using sort (order + match) to make sure that Old vs Young and Young vs Old both become Young vs Old
                df$GroupvsGroup <- apply(df[, c("Row_Group", "Col_Group")], 1, function(x){
                        paste(x[order(match(x, levels(samdf$Group)))], collapse = " vs ")
                })
                
                # to make sure group comparisons are in the "order" of levels "group_var"
                df$GroupvsGroupOrder <- apply(df[, c("Row_Group", "Col_Group")], 1, function(x){
                        x1 <- match(x[1], levels(samdf$Group))
                        x2 <- match(x[2], levels(samdf$Group))
                        if (x1 <= x2) {
                                as.numeric(paste(x1, x2, sep = ""))
                        } else {
                                as.numeric(paste(x2, x1, sep = ""))
                        }
                })
                
                df$Type <- "between"
                df$Type[df$Row_Group == df$Col_Group] <- "within"
                # NB: for within group distances, you have samples*(samples-1)/2 distances, for between group distances
                # you have samplesgrp1 * samplesgrp2 distances.
                
                GvsGLeveldf <- unique(df[, c("GroupvsGroup", "GroupvsGroupOrder", "Type")])
                GroupvsGroupLevels <- GvsGLeveldf[order(GvsGLeveldf$GroupvsGroupOrder), ]$GroupvsGroup
                GvsGLeveldfBetween <- GvsGLeveldf[GvsGLeveldf$Type == "between", ]
                GroupvsGroupLevelsBetween <- GvsGLeveldfBetween[order(GvsGLeveldfBetween$GroupvsGroupOrder), ]$GroupvsGroup
                df$GroupvsGroup <- factor(df$GroupvsGroup, levels = GroupvsGroupLevels, ordered = TRUE)
                
                pairwise_dist_ttest <- pairwise.t.test(x = df$Distance, g = df$GroupvsGroup, alternative = "two", p.adjust.method = "none", var.equal = F, pool.sd = F)
                pairwise_dist_ttest <- pairwise_dist_ttest$p.value
                rowCol <- expand.grid(rownames(pairwise_dist_ttest), colnames(pairwise_dist_ttest))
                labs <- rowCol[as.vector(lower.tri(pairwise_dist_ttest,diag=T)),]
                df_p <- cbind(labs, pairwise_dist_ttest[lower.tri(pairwise_dist_ttest,diag=T)])
                colnames(df_p) <- c("Distances_2","Distances_1","p_value")
                df_p <- df_p[, c("Distances_1", "Distances_2", "p_value")]
                df_p$Significance <- ""
                df_p$Significance[df_p$p_value <= 0.05] <- "*"
                df_p$Significance[df_p$p_value <= 0.01] <- "**"
                df_p$Significance[df_p$p_value <= 0.001] <- "***"
                
                pValList[[i]] <- df_p
                
                
                # in case there are more than two groups in group_var I want a faceted plot, i.e. one facet for each GroupvsGroupLevelsBetween, for this I need to duplicate some data
                
                if (length(GroupvsGroupLevels) > 1) {
                        df_list <- lapply(GroupvsGroupLevelsBetween, function(level){
                                grps <- unlist(strsplit(level, " vs "))
                                df_current <- df[df$Row_Group %in% grps & df$Col_Group %in% grps, ]
                                df_current$Level <- level
                                df_current
                        })
                        df_plot <- do.call(rbind, df_list)
                        df_plot$Level <- factor(df_plot$Level, levels = GroupvsGroupLevelsBetween, ordered = TRUE)
                        
                        Tr <- ggplot(df_plot, aes(x = GroupvsGroup, y = Distance, col = GroupvsGroup))
                        Tr <- Tr + geom_boxplot(outlier.color = NA) + 
                                geom_jitter(position = position_jitter(width = 0.3, height = 0), alpha = 0.65) +
                                facet_wrap(~ Level, ncol = 2, scales = "free_x") +
                                theme_bw() +
                                xlab("") +
                                ylab(paste(names(dist_list)[i], "distance", sep = " ")) +
                        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                              legend.title = element_blank())
                        if (length(GroupvsGroupLevels) <= 7) {
                                Tr <- Tr +
                                        scale_color_manual(values = cbPalette[2:8])
                        }
                        
                } else {
                        Tr <- ggplot(df, aes(x = GroupvsGroup, y = Distance, col = GroupvsGroup))
                        Tr <- Tr + geom_boxplot(outlier.color = NA) + 
                                geom_jitter(position = position_jitter(width = 0.3, height = 0), alpha = 0.65) +
                                theme_bw() +
                                xlab("") +
                                ylab(paste(names(dist_list)[i], "distance", sep = " ")) + 
                                theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                                      legend.title = element_blank())
                        if (length(GroupvsGroupLevels) <= 7) {
                                Tr <- Tr +
                                        scale_color_manual(values = cbPalette[2:8])
                        }
                        
                }
                
                TrList[[i]] <- Tr
        }
        
        names(TrList) <- names(dist_list)
        names(pValList) <- names(dist_list)
        out <- list(DistanceBoxplots = TrList, DistancePValues = pValList)
        
}


#######################################
### FUNCTION: pairwise.perm.manova.own
#######################################
# Function is very much based on pairwise.perm.manova {RVAideMemoire}
# but it also records R2 while looping through vegan::adonis, and generates a
# result data frame in which the results are shown in the order of the group_fac levels
# INPUT:
# dist_obj: dist object
# group_fac: the factor that groups the samples
# nperm: permutations in adonis
# p.adj.method: method to adjust p.values
# OUTPUT:
# data.frame showing p.values and R2 and adjusted p.values for the different between group comparisons


pairwise.perm.manova.own <- function(dist_obj, group_fac, nperm = 999, 
                                     p.adj.method = "none") {
        
        if (!("dist" %in% class(dist_obj))){
                stop("dist_obj must be of class dist")
        }
        
        group_fac <- factor(group_fac)
        
        fac_levels <- levels(group_fac)
        fac_levels_num <- setNames(seq_along(fac_levels), fac_levels) # see pairwise.table
        
        i_s <- outer(fac_levels_num, fac_levels_num, function(ivec, jvec){
                sapply(seq_along(ivec), function(x){
                        i <- ivec[x]
                })
        })
        j_s <- outer(fac_levels_num, fac_levels_num, function(ivec, jvec){
                sapply(seq_along(ivec), function(x){
                        j <- jvec[x]
                })
        })
        i_s <- i_s[upper.tri(i_s)]
        j_s <- j_s[upper.tri(j_s)]
        
        p_vals <- vector(mode = "numeric", length = length(i_s))
        r2s <- vector(mode = "numeric", length = length(i_s))
        
        for (k in seq_along(i_s)){
                i <- i_s[k]
                j <- j_s[k]
                group_fac2 <- droplevels(group_fac[as.numeric(group_fac) %in% c(i, j)])
                dist_obj_mat <- as.matrix(dist_obj)
                rows <- which(group_fac %in% levels(group_fac2))
                dist_obj2 <- as.dist(dist_obj_mat[rows, rows])
                fit <- vegan::adonis(dist_obj2 ~ group_fac2, permutations = nperm)
                p_vals[k] <- fit$aov.tab[1, "Pr(>F)"]
                r2s[k] <- fit$aov.tab[1, "R2"]
        }
        
        
        
        result_df <- data.frame(Comparison = paste0(names(fac_levels_num[i_s]), "_vs_", names(fac_levels_num[j_s]), sep = ""),
                                addonis_p_value = p_vals, adonis_R2 = r2s, p_val_adj = p.adjust(p_vals, p.adj.method))
        
}


#######################################
### overviewPhyseq##
#################


overviewPhyseq <- function(physeq, group_var, max_rel_ab_for_color = NULL){
        
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
        
        
        # # == Tr3: Histogram of template sample SD on relative abundances excluding zeros == 
        # Temp <- data.frame(Sample = sample_names(template), Sample_SD_RA_NZ = 100*apply(OTUTTNA_RA, 1, sd, na.rm = TRUE))
        # Tr3 <- ggplot(Temp, aes(x = Sample_SD_RA_NZ))
        # Tr3 <- Tr3 + geom_histogram(binwidth = diff(range(Temp$Sample_SD_RA_NZ))/30, fill = cbPalette[2]) +
        #         geom_rug() +
        #         geom_vline(xintercept = median(Temp$Sample_SD_RA_NZ), col = "#009E73", lty = "dashed") +
        #         ggtitle(paste("Sample SD (RA, NZ), distribution of ", nsamples(template), " template samples (median green line)", sep = "")) +
        #         theme_bw() +
        #         xlab("Sample_SD_RA_NZ")
        
        
        
        # == Heatmap plot ==
        
        # here Samples should be columns and taxa rows
        # zeros should be of special color to see sparsity
        # then to have the option to give Count limits, so that every count above that value gets the max colour
        # The latter point can be done with the scale package and the option: oob = swish
        # please read here: <https://github.com/tidyverse/ggplot2/issues/866>
        # see also the na.value option not used
        
        DF_CT <- as.data.frame(t(CT))
        rownames(DF_CT) <- paste("Taxon_", 1:nrow(DF_CT), sep = "")
        DF_CT$Taxa <- rownames(DF_CT)
        DF_CT <- tidyr::gather(DF_CT, key = Sample , value = Count, -Taxa)
        DF_CT$Taxa <- factor(DF_CT$Taxa, levels = rev(unique(DF_CT$Taxa)), ordered = TRUE)
        LookUpDF <- data.frame(Sample = sample_names(physeq), Group = sample_data(physeq)[[group_var]])
        # I later want the samples coloured by group_var
        LookUpDF <- LookUpDF[order(match(LookUpDF$Group, levels(LookUpDF$Group))), ]
        DF_CT$Sample <- factor(DF_CT$Sample, levels = LookUpDF$Sample, ordered = TRUE)
        
        # add color to the taxa so you see up and down
        if (length(levels(LookUpDF$Group)) <= 7){
                color_lookup <- data.frame(level = levels(LookUpDF$Group), color = cbPalette[2:(length(levels(LookUpDF$Group)) + 1)])
                colxaxis <- as.character(color_lookup$color[match(LookUpDF$Group, color_lookup$level)])
        } else {
                colxaxis <- rep("black", nrow(LookUpDF))
        }
        
        
        DF_CT_RA <- dplyr::group_by(DF_CT, Sample) %>% dplyr::mutate(rel_ab = Count/sum(Count))
        
        if (is.null(max_rel_ab_for_color) || max_rel_ab_for_color < 0 || max_rel_ab_for_color > 1) {
                max_rel_ab_for_color <- max(DF_CT_RA$rel_ab)
        }
        
        # Test that really only 0 entries are shown in red:
        # DF_CT_RA$rel_ab[DF_CT_RA$Taxa %in% c("Taxon_1", "Taxon_2", "Taxon_3")] <- 1e-12 # will show you it only gets red as soon as it is <- 1e-14
        # Test that values above max_rel_ab_for_color are yellow
        # DF_CT_RA$rel_ab[DF_CT_RA$Taxa %in% c("Taxon_1", "Taxon_2", "Taxon_3")] <- .9
        
        
        hmtRA <- ggplot(DF_CT_RA, aes(x = Sample, y = Taxa, fill = rel_ab))
        hmtRA <- hmtRA + 
                geom_raster() + 
                scale_fill_gradientn(limits = c(0, max_rel_ab_for_color), colors = c("red", viridis(5)), values = c(0, 1e-14, 0.15, 0.3, 0.45, 1), oob = squish) +
                scale_x_discrete(position = "top") +
                #coord_equal() +
                labs(x=NULL, y=NULL) +
                theme_tufte(base_family = "Helvetica") +
                theme(axis.ticks=element_blank(),
                      axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0,
                                                 colour = colxaxis))
        
        list(OverviewDF = Overview, heatmap_ra = hmtRA)
        
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