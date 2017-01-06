data("GlobalPatterns")
GP = GlobalPatterns


# go to relative abundance and filter taxa
GPr  = transform_sample_counts(GlobalPatterns, function(x) x / sum(x) )
GPfr = filter_taxa(GPr, function(x) mean(x) > 1e-5, TRUE)

# only keep taxa of Phylum Chlamydiae
GP.chl = subset_taxa(GP, Phylum=="Chlamydiae")
GP.chl = prune_samples(sample_sums(GP.chl)>=20, GP.chl) # why no filter_samples?
# go to relative abundance
transform_sample_counts(GP.chl, function(OTU) OTU/sum(OTU))

### Combine based on taxonomic rank
## the tax_glom function is what you want when combining based on taxonomic ranks. However, function takes long time
unique(tax_table(GP)[, "Phylum"]) == "Bacteroidetes"
gpsfb = subset_taxa(GP, Phylum == "Bacteroidetes")
unique(tax_table(gpsfb)[, "Phylum"])
# all Bacteroidetes
## How about the OTU counts
identical(otu_table(gpsfb), otu_table(GP)[rownames(otu_table(gpsfb)),]) # TRUE
# so no summing here or things like that

# Now tax_glom at the "Family" level
length(unique(tax_table(gpsfb)[, "Family"])) # I expect 16 taxa after tax_glom
gpsfbg = tax_glom(gpsfb, "Family")
ntaxa(gpsfbg) # only 15, guess NA ignored?
intersect(unique(tax_table(gpsfbg)[, "Family"]), unique(tax_table(gpsfb)[, "Family"])) # 15 entries NA missing
setdiff(unique(tax_table(gpsfbg)[, "Family"]), unique(tax_table(gpsfb)[, "Family"])) # settdiff ignores the NA too
Uniquegpsfbg <- unique(tax_table(gpsfb)[, "Family"]) # NB it is still of class taxonomy table!
setdiff(rownames(unique(tax_table(gpsfbg)[, "Family"])), rownames(unique(tax_table(gpsfb)[, "Family"])))
# shows another possible issue, that the rownames are not consistent:
taxa_names(gpsfbg) %in% taxa_names(gpsfb) # but all would be in gpsfb

## check the OTU tables if it is summing up correctly
OT <- as.data.frame(otu_table(gpsfb))
TT <- as.data.frame(tax_table(gpsfb))
OT$Family <- TT$Family
unique(OT$Family) # only 15 levels, one is NA
OT <- group_by(OT, Family)
OTf <- summarise_each(OT, funs(sum))
# now compare
OTff <- as.data.frame(otu_table(gpsfbg))
# add Family 
OTff$Family <- as.data.frame(tax_table(gpsfbg))$Family
OTf <- arrange(OTf, Family)
OTff <- arrange(OTff, Family)
OTff <- select(OTff, Family, 1:26)
# and indeed
all.equal(OTf[-16, ], OTff) # TRUE!


### Further preprocessing
###### Filtering 
GP = filter_taxa(GlobalPatterns, function(x) sum(x > 3) > (0.2*length(x)), TRUE) # keeps only taxa that have in at least 20% of samples more than 3 reads
# down to 2572 taxa
## add the human variable
sample_data(GP)$human = factor( get_variable(GP, "SampleType") %in% c("Feces", "Mock", "Skin", "Tongue") )

############ Standardize abundances
total = median(sample_sums(GP))
standf = function(x, t=total) round(t * (x / sum(x)))
gps = transform_sample_counts(GP, standf)
range(sample_sums(gps)) # now all very similar


####### Einschub DESeq2 ####################
### Standardize a la DESEQ using geometric mean
#NB: check the percentage of 0 in your samples
sample_sums(otu_table(GP)==0)/ntaxa(GP)
# NB: DESeq2 ignores when taking the median of the ratios for size factor determination 0 values!

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
# # NB the prod(x/gm(x)) = 1 only applies when no 0 present, e.g.
# x <- 1:10
# prod(x/gm(x))

if(taxa_are_rows(GP)){
        GM <- apply(otu_table(GP), 1, gm)   
} else {
        GM <- apply(otu_table(GP), 2, gm) 
}

# calculate Size Factors for the samples
if(taxa_are_rows(GP)){
        SF <- apply(sweep(otu_table(GP), 1, GM, "/"), 2, function(x){median(x[x!=0])})   
} else {
        SF <- apply(sweep(otu_table(GP), 2, GM, "/"), 1, function(x){median(x[x!=0])}) 
}
## again ignores all 0 here!
sample_sums(otu_table(GP)==0)/ntaxa(GP)

# do it anyway:
GPDE <- GP
if(taxa_are_rows(GPDE)){
        otu_table(GPDE) <- otu_table(sweep(otu_table(GPDE), 2, SF, "/"), taxa_are_rows = TRUE)
} else {
        otu_table(GPDE) <- otu_table(sweep(otu_table(GPDE), 1, SF, "/"), taxa_are_rows = FALSE) 
}

### Compare directly to DESeq2
GPDES = phyloseq_to_deseq2(GP, ~ human)
OTUBefore <- counts(GPDES)
#GM <- apply(counts(GPDES), 1, gm)
SIZE <- estimateSizeFactors(GPDES, geoMeans = GM)
SFDES <- sizeFactors(SIZE)
# NOT SURE WHY SF is sometimes minimally bigger than SFDES
SF/SFDES
# COMPARE OTUs
OTU <- as.data.frame(counts(SIZE, normalized = TRUE))
OTUSelf <- as.data.frame(otu_table(GPDE))
all.equal(OTU, OTUSelf) #Tiny tiny differences
# go with DESeq2 back to phyloseq
GPDES <- phyloseq(otu_table(OTU, taxa_are_rows = TRUE), 
                            tax_table(GP),
                            sample_data(GP))
sample_sums(GPDES)
all.equal(sample_sums(GPDE), sample_sums(GPDES))
## NB: far from similar now!!!!
############################# Ende Einschub ########

## Filter the taxa using a cutoff of 3.0 for the Coefficient of Variation
## so they only want to keep taxa that vary

gpsf = filter_taxa(gps, function(x) sd(x)/mean(x) > 3.0, prune = TRUE) # down to 1201 taxa

# also keep a physeq with only Bacteroidetes
gpsfb = subset_taxa(gpsf, Phylum=="Bacteroidetes")

title = "plot_bar; Bacteroidetes-only"
plot_bar(gpsfb, x = "SampleType", y = "Abundance", title=title)
# see indasettelse below, basically, if you plotted sample vs Abundance and had corrected all bars should be similarly long

####### indsaettelse understanding plot_bar ##########
# the function is based on psmelt which is again based on the melt function from reshape2
# which is without possible conflicts simply
# rankNames = rank_names(physeq, FALSE)
# sampleVars = sample_variables(physeq, FALSE)
# reservedVarnames = c("Sample", "Abundance", "OTU")
# otutab = as(otu_table(physeq), "matrix")
# mdf = reshape2::melt(otutab) # generates long format (nxm of matrix in nrow) with the following three columns
# colnames(mdf)[1] <- "OTU"
# colnames(mdf)[2] <- "Sample"
# colnames(mdf)[3] <- "Abundance"
# mdf$OTU <- as.character(mdf$OTU)
# mdf$Sample <- as.character(mdf$Sample)
# sdf = data.frame(sample_data(physeq), stringsAsFactors = FALSE) # nrows = nsamples
# sdf$Sample <- sample_names(physeq)
# mdf <- merge(mdf, sdf, by.x = "Sample") # merge function nice
# # add the taxonomic ranks
# TT = access(physeq, "tax_table")
# keepTTcols <- colSums(is.na(TT)) < ntaxa(TT) # not keeping taxonomic levels that are all NA
# if (length(which(keepTTcols)) > 0 & ncol(TT) > 0) {
#         TT <- TT[, keepTTcols]
#         tdf = data.frame(TT, OTU = taxa_names(physeq)) # nrow = ntaxa(physeq)
#         mdf <- merge(mdf, tdf, by.x = "OTU")}
# mdf = mdf[order(mdf$Abundance, decreasing = TRUE), ]
# # return(mdf) so that was psmelt
# # mdf is a long data frame having all combis of Sample OTU Abundance in
# x <- "Sample"
# y <- "Abundance"
# fill <- NULL
# p = ggplot(mdf, aes_string(x = x, y = y, fill = fill))
# p = p + geom_bar(stat = "identity", position = "stack", color = "black")
# p = p + theme(axis.text.x = element_text(angle = -90, hjust = 0))
# if (!is.null(facet_grid)) {
#         p <- p + facet_grid(facet_grid)
# }
# if (!is.null(title)) {
#         p <- p + ggtitle(title)
# }
# return(p)
############################### end indsaettelse ###############

plot_bar(gpsfb, x = "SampleType", y = "Abundance", fill = "Family", title=title)
plot_bar(gpsfb, "Family", "Abundance", "Family", 
         title=title, facet_grid="SampleType~.")
