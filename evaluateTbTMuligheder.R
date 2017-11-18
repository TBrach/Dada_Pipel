# saved ps_filt_extra of the phylum level comparison

saveRDS(object = ps_filt_extra, "ps_filt_extra.rds")

# so load object here
ps_filt_extra <- readRDS(file = "ps_filt_extra.rds")
group_var = "Group"
functionpath <- "/Users/jvb740/MarieCurie_Work/BackgroundKnowledge/16S_Learning/Dada_Pipel/Functions/"
source(file.path(functionpath, "Dada_TaxonomyFunctions.R"))

raw_TbTmatrixes_list <- calculate_raw_TbTmatrixes(ps_filt_extra, group_var = group_var)
# REMEMBER: # 0/x = 0, x/0 = Inf; 0/0 = NaN!

# get only the TbTmatrixes for Young_vs_Old
TbTmatrixes <- raw_TbTmatrixes_list[[2]]
# define the group_fac for Young_vs_Old
group_fac <- sample_data(ps_filt_extra)[[group_var]]
group_fac <- droplevels(group_fac[as.numeric(group_fac) %in% c(1, 3)])

#now work on rank_evaluate_TbTmatrixes in Dada_TaxonomyFunctions
i = 1
mat <- TbTmatrixes[[i]]
# REMEMBER: # 0/x = 0, x/0 = Inf; 0/0 = NaN!
# calculate ranks, and transform NaN (i.e. where both host taxon and denominator taxon were absent to NA)

mat_rank <- t(apply(mat, 1, rank, na.last = "keep"))
# NB: in case no NA then rowSums(mat_rank) is equal to ncol(mat)*(ncol(mat) + 1)/2
# the smallest ratio is rank 1, the highest ratio is Max, so the higher the number the higher was the host taxon in comparison to this taxon
# If NA are present, then note that the ranks are smaller, the sum is only (nsamples that are not NA * nsamples that are not NA + 1) / 2

# it is therefore that I divide by the geometric mean and log
# NB: log(x/gm) = log(x) - log(gm)
mat_rank_divgm_log <- t(apply(mat_rank, 1, function(taxa_ranks){
        log(taxa_ranks) - log(gm_own_tbt(taxa_ranks))
}))

# that should have kept all NAs, and rowSums should be 0
all.equal(max(rowSums(mat_rank_divgm_log), na.rm = T), 0) # should be True

# How to continue? Think about a sample in which the host taxon is not present, The ratio was 0 for each taxon that was present, and NA for all others
# if less than 50% of taxa were present in that sample, than the median would result in NA, if not it would result in a negative number
# worse: if the host taxon was present, but the majority of other taxa wasn't, then when taking the median, you will get a high value caused by the simple fact
# that most of the taxa was not present. 



# next steps:
# - divide by geometric mean and log so rowSums are 0 in each row
# - set all NA to 0
# - get colSums to get evaluation matrix1
# - run wilcoxon test