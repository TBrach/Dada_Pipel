data("GlobalPatterns")
GP = GlobalPatterns
GP = prune_taxa(taxa_sums(GlobalPatterns) > 0, GlobalPatterns)
humantypes = c("Feces", "Mock", "Skin", "Tongue")
sample_data(GP)$human <- get_variable(GP, "SampleType") %in% humantypes


######## merge_samples 
# sums up abundances for samples of same type
mergedGP = merge_samples(GP, "SampleType")
# NB merge_samples changed taxa_are_rows
taxa_are_rows(GP)
taxa_are_rows(mergedGP)
SD = merge_samples(sample_data(GP), "SampleType")
SD1 <- sample_data(mergedGP)
identical(SD, SD1) #TRUE
print(SD[, "SampleType"])
sample_names(GP)
sample_names(mergedGP)

#As emphasized earlier, the OTU abundances of merged samples are summed. Let’s investigate this ourselves looking at just the top10 most abundance OTUs.
OTUnames10 = names(sort(taxa_sums(GP), decreasing = TRUE)[1:10])
GP10  = prune_taxa(OTUnames10,  GP)
mGP10 = prune_taxa(OTUnames10, mergedGP)
ocean_samples = sample_names(subset(sample_data(GP), SampleType=="Ocean"))
taxa_are_rows(GP10) #TRUE > samples are columns
otu_table(GP10)[, ocean_samples]
rowSums(otu_table(GP10)[, ocean_samples])
taxa_are_rows(mGP10) # FALSE > samples are rows
otu_table(mGP10)["Ocean", ] # same numbers, indeed abundances were summed

## RICHNESs plots
plot_richness(GP, x = "human", color = "SampleType", title="unmerged")

SD <- sample_data(GP)
SD1 <- sample_data(mergedGP)
## they even recommend checking the sample data
identical(sample_data(mergedGP)$SampleType, sample_names(mergedGP)) #FALSE, first is 1,2,3,4
sample_data(mergedGP)$SampleType = sample_names(mergedGP)
sample_data(mergedGP)$human = sample_names(mergedGP) %in% humantypes
plot_richness(mergedGP, "human", "SampleType", title="merged")

### Summary: Not sure I want to use merge_samples, but can be used to add up abundances within samples

######### merge_taxa
### also note: agglomeration functions tip_glom or tax_glom that merge similar OTUs based on a phylogenetic or taxonomic threshold

# could not get the closedls

## you can selectively choose which taxa to merge!
x1 = merge_taxa(closedps, taxa_names(closedps)[3:27], 2)

## but how to merge based on taxonomic rank
## SEE tax_glom and tip_glom here, (in Preprocessing tutorial)

######### merge_phyloseq
# see tutorial, is to add for example sample data, but you could always construct a completely new object


