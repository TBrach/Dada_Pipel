# see <https://joey711.github.io/phyloseq/distance.html>

# The distance function takes a phyloseq-class object and method option, 
# and returns a dist-class distance object suitable for certain ordination methods
# and other distance-based analyses. There are currently 44 explicitly supported method options in the phyloseq package,
## NB there is also an option for user-provided arbitrary methods via an interface to vegan::designdist.

# complet list of options/arguments for the method parameter
distanceMethodList

# Only sample-wise distances are currently supported (the type argument), but eventually OTU-wise (e.g. species) distances will be supported as well.

?phyloseq::distance

## example physeq
data(GlobalPatterns)
GP <- GlobalPatterns
ntaxa(GP) # 19216

## Filtering the GP to reduce processing time
# keep only human samples
sample_data(GP)$human = factor( get_variable(GP, "SampleType") %in% c("Feces", "Mock", "Skin", "Tongue") )

GPh <- subset_samples(GP, subset = human == TRUE )

## keep only taxa present in at least 30% of samples
GPhf <- filter_taxa(GPh, function(x){sum(x != 0) > .3*length(x)}, prune = TRUE) # still 3670 taxa

## do their normalisation
total = median(sample_sums(GPhf))
standf = function(x, t=total) round(t * (x / sum(x)))
GPhfs = transform_sample_counts(GPhf, standf)
range(sample_sums(GPhfs))

## change to relative abundance here
GPhfs = transform_sample_counts(GPhfs, function(x){x/sum(x)})


## keep only taxa with good variation
GPhfsf = filter_taxa(GPhfs, function(x) sd(x)/mean(x) > 3.0, prune = TRUE) # down to 388 taxa
sample_sums(GPhfsf) # again very different but I think fine since all have now same taxa

plot_bar(GPhfsf, x = "Sample", y = "Abundance", fill = "SampleType")
plot_bar(GPhfsf, x = "SampleType", y = "Abundance", fill = "Phylum")

theme_set(theme_bw())

##### I now do their enterotypes example but for this Global Pattern set 
# NB: I ran into trouble figuring out that many of these distances appear to need relative abundances!
# therefore I added the relative abundance above

dist_methods <- unlist(distanceMethodList)
# I reduce
dist_methods <- dist_methods[dist_methods %in% c("unifrac", "wunifrac", "jsd", "canberra", "bray", "jaccard")]

# NB: UniFrac requires a tree! GP has a tree
## I save plots using MDS/PCoA (test if difference) and NMDS!
plist <- vector("list", length(dist_methods))
names(plist) = dist_methods
plistPCoA <- vector("list", length(dist_methods))
names(plistPCoA) = dist_methods # NB: CAN BE LEFT OUT SAME AS MDS!!
plistNMDS <- vector("list", length(dist_methods))
names(plistNMDS) = dist_methods
for( i in dist_methods ){
        # Calculate distance matrix
        iDist <- distance(GPhfsf, method=i)
        
        ## NB her: as(iDist, "matrix") # is a nsamplesxnsamples matrix
        
        # Calculate ordination
        iMDS  <- ordinate(GPhfsf, "MDS", distance=iDist) # 
        iMDSPCoA  <- ordinate(GPhfsf, "PCoA", distance=iDist) #
        iMDSNMDS  <- ordinate(GPhfsf, "NMDS", distance=iDist)
        ## Make plot
        # Don't carry over previous plot (if error, p will be blank)
        p <- NULL
        # Create plot, store as temp variable, p
        p <- plot_ordination(GPhfsf, iMDS, color="SampleType")
        # Add title to each plot
        p <- p + ggtitle(paste("MDS using distance method ", i, sep=""))
        # Save the graphic to file.
        plist[[i]] = p
        ## Make plot
        # Don't carry over previous plot (if error, p will be blank)
        p1 <- NULL
        # Create plot, store as temp variable, p
        p1 <- plot_ordination(GPhfsf, iMDSPCoA, color="SampleType")
        # Add title to each plot
        p1 <- p1 + ggtitle(paste("PCoA using distance method ", i, sep=""))
        # Save the graphic to file.
        plistPCoA[[i]] = p1
        ## Make plot
        # Don't carry over previous plot (if error, p will be blank)
        p2 <- NULL
        # Create plot, store as temp variable, p
        p2 <- plot_ordination(GPhfsf, iMDSNMDS, color="SampleType")
        # Add title to each plot
        p2 <- p2 + ggtitle(paste("NMDS using distance method ", i, sep=""))
        # Save the graphic to file.
        plistNMDS[[i]] = p2
}

# NB arranging the plots could be done with grid.arrange This would keep the distance info
do.call("grid.arrange", c(plist, nrow=3))
# equals: grid.arrange(plist[[1]], plist[[2]], plist[[3]], plist[[4]], plist[[5]], plist[[6]],nrow = 3)
## NB
do.call("grid.arrange", c(plistPCoA, nrow=3)) # are indeed fully the same as plist


# his plyr version is nide to use facet wrap but the axis info is lost
library(plyr)

df = ldply(plist, function(x) x$data)
# basically you could loop through plist and rbind plist[[1]]$Data than add the .id from the names
names(df)[1] <- "distance"
p = ggplot(df, aes(x = Axis.1, y = Axis.2, color=SampleType))
p = p + geom_point(size=3, alpha=0.5)
p = p + facet_wrap(~distance, scales="free")
p

# NB: The since PCoA was fully the same as MDS, comparing just PCoA vs NMDS
do.call("grid.arrange", c(plistPCoA, plistNMDS, nrow=2))

# try a plyr version
df = ldply(plistPCoA, function(x) x$data)
df1 = ldply(plistNMDS, function(x) x$data)
df$Method <- "PCoA"
df1$Method <- "NMDS"
names(df1)[2:3] <- names(df)[2:3]
names(df)[1] <- "distance"
names(df1)[1] <- "distance"
dfC <- rbind(df, df1)

p = ggplot(dfC, aes(x = Axis.1, y = Axis.2, color=SampleType))
p = p + geom_point(size=3, alpha=0.5)
p = p + facet_wrap(Method~distance, scales="free")
p

