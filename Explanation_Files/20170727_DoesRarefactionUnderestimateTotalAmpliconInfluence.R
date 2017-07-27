# type = "pois" # or "nb"
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
size1 <- 50000
size2 <- 15000
seed <- 154
nsims <- 100

# - use a base sample -
seqtab <- as(otu_table(ps), "matrix")
base_sample <- seqtab[which.max(rowSums(seqtab != 0)),]
# the sample has 308 SVs above 0
sum(base_sample > 0)
# the lowest count is 4
min_count <- min(base_sample[base_sample >0])
# --


# - generate a sample with 400 SVs -
set.seed(seed)

DNA_start_richness <- vector("numeric", length = nsims)
S1_richness <- vector("numeric", length = nsims)
S2_richness <- vector("numeric", length = nsims)
S3_richness <- vector("numeric", length = nsims)

# extra 
S1_min <- vector("numeric", length = nsims)
S1_max <- vector("numeric", length = nsims)
S1_minsec <- vector("numeric", length = nsims)
S1_sd <- vector("numeric", length = nsims)

for (i in 1:nsims) {
        # let's assume the real richness would be 400, so get 92 SVs more with "counts" below 4
        factors <- abs(rnorm(n = 92, sd = 2))
        factors <- factors/max(factors)
        extra_low_SVs <- factors*min_count
        
        # sample_prop <- c(extra_low_SVs, base_sample[base_sample > 0])
        sample_prop <- base_sample[base_sample > 0]
        sample_prop <- sample_prop/sum(sample_prop)
        
        # and after a fair successful PCR
        DNA_start_sample <- round(sample_prop*1e6)
        DNA_start_richness[i] <- sum(DNA_start_sample > 0)
        
        S1 <- rarefy_sample(DNA_start_sample, size = size1)
        S2 <- rarefy_sample(DNA_start_sample, size = size2)
        S3 <- rarefy_sample(S1, size = size2)
        S1_richness[i] <- sum(S1 > 0)
        S2_richness[i] <- sum(S2 > 0)
        S3_richness[i] <- sum(S3 > 0)
        
        S1_min[i] <- min(S1)
        S1_max[i] <- max(S1)
        S1_minsec[i] <- sort(unique(S1))[2]
        S1_sd[i] <- sd(S1)
}

DF <- data.frame(DNA = DNA_start_richness, S1 = S1_richness, S2 = S2_richness, S3 = S3_richness)

DFl <- gather(DF, key = sample, value = richness)
DFl$Sample <- factor(DFl$sample, levels <- c("DNA", "S1", "S2", "S3"), ordered = TRUE)

pairwise.t.test(x = DFl$richness, g = DFl$Sample, alternative = "two", p.adjust.method = "none", var.equal = F, pool.sd = F)

Tr <- ggplot(DFl, aes(x = sample, y = richness, col = sample))

Tr <- Tr +
        geom_boxplot() +
        geom_jitter() +
        scale_color_manual("", values = cbPalette[2:8]) +
        theme_bw(12)
        
        



# --------------------- other sample attempts ----------------------------------


# if (type == "pois") {
#         proportions_real <- rpois(n = 500, lambda = .1)
#         proportions_real <- proportions_real/sum(proportions_real)
#         sum(proportions_real == 0)
#         proportions_real[proportions_real == 0] <- 0.0005
#         # add two random large samples
#         proportions_real[sample(500, 2)] <- 0.12
#         proportions_real <- proportions_real/sum(proportions_real)
#         hist(proportions_real)
#         range(proportions_real)
# }
# 
# 
# # if (type == "nb"){
# #         proportions_real <- rnbinom(n = 500, size = 40, prob = 0.1)
# #         proportions_real <- proportions_real/sum(proportions_real)
# #         range(proportions_real)
# #         sum(proportions_real == 0)
# #         hist(proportions_real)
# # }
# 
# # --