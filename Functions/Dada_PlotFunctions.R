#######################################
#### IndivQualPlot
#######################################

IndivQualPlot <- function(QStatsList, SampleName, xlim_low = -5, xlim_high = 255) {
        
        cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
        
        df <- QStatsList[[SampleName]]
        ggplot(data = df, aes(x = Cycle, y = Mean_QS))  + 
                geom_hline(yintercept = 25, color = 'darkred', linetype = "dashed", size = .35) +
                geom_line(color = "#009E73") + 
                geom_line(aes(y = Median_QS), color = "#E69F00") +
                #geom_line(aes(y = q25_QS), color = "#FC8D62", linetype = "dashed") + 
                #geom_line(aes(y = q75_QS), color = "#FC8D62", linetype = "dashed") +
                ylab("Quality Score") + 
                xlab("Cycle") +
                coord_cartesian(xlim=c(xlim_low, xlim_high),
                                ylim = c(8,42), expand = FALSE) +
                ## I added this
                ggtitle(paste("Sample:", SampleName, " No Reads: ", df$NoReads[1])) +
                theme_bw() + 
                theme(panel.grid.minor = element_blank(),
                      panel.grid.major.x = element_blank(),
                      panel.grid.major.y = element_line(color = "#999999", size = .15))
        
}

#######################################
#### QS_Median_OverviewPlot
#######################################

QS_Median_OverviewPlot <- function(QStatsList, SampleNames, Prefix = "FW", xlim_low = -5, xlim_high = 255) {
        # Input:
        # QStatsList: list of QStats data frames, such as FW_QualityStats
        # SampleNames: List of Sample names that must be the names of the QStatsList
        
        if (!all(SampleNames %in% names(QStatsList))) {
                stop("not all SampleNames were found in the QStatsList")
        }
        
        ReadLengths <- sapply(1:length(SampleNames), function(x) dim(QStatsList[[SampleNames[x]]])[1])
        x <- range(ReadLengths)/mean(ReadLengths)
        if (!isTRUE(all.equal(x[1], x[2], tolerance = .Machine$double.eps ^ 0.5))) {
                
                stop("Samples are of different read lengths") 
        }
        
        cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
        
        QList <- QStatsList[SampleNames]
        # add sample names to the data frames
        for (i in seq_along(QList)) {
                QList[[i]]$Sample <- names(QList[i])
        }
        # combine all df in the list to one big df
        df <- do.call("rbind", QList)
        
        ggplot(data = df, aes(x = Cycle, y = Median_QS, colour = NoReads, group = Sample))  + 
                geom_hline(yintercept = 25, color = 'darkred', linetype = "dashed", size = .35) +
                geom_line() +
                #geom_point() +
                ylab("Quality Score") + 
                xlab("Cycle") +
                coord_cartesian(xlim=c(xlim_low, xlim_high),
                                ylim = c(8,42), expand = FALSE) +
                scale_color_gradient2(name = "No reads", limits = c(10000, 60000), midpoint = 35000, low = "#009E73", high = "#E69F00", mid = "#999999") +
                #scale_color_gradient2(name = "No reads", limits = c(min(df$NoReads), max(df$NoReads)), midpoint = min(df$NoReads) + (range(df$NoReads)[2]- range(df$NoReads)[1])/2, low = "#009E73", high = "#E69F00", mid = "#999999")
                ggtitle(paste(Prefix, "reads, Median_QS, NoSamples:", length(SampleNames))) +
                theme_bw() + 
                theme(panel.grid.minor = element_blank(),
                      panel.grid.major.x = element_blank(),
                      panel.grid.major.y = element_line(color = "#999999", size = .15))
}

#######################################
#### QS_Mean_OverviewPlot
#######################################

QS_Mean_OverviewPlot <- function(QStatsList, SampleNames, Prefix = "FW", xlim_low = -5, xlim_high = 255) {
        # Input:
        # QStatsList: list of QStats data frames, such as FW_QualityStats
        # SampleNames: List of Sample names that must be the names of the QStatsList
        
        if (!all(SampleNames %in% names(QStatsList))) {
                stop("not all SampleNames were found in the QStatsList")
        }
        
        ReadLengths <- sapply(1:length(SampleNames), function(x) dim(QStatsList[[SampleNames[x]]])[1])
        x <- range(ReadLengths)/mean(ReadLengths)
        if (!isTRUE(all.equal(x[1], x[2], tolerance = .Machine$double.eps ^ 0.5))) {
                
                stop("Samples are of different read lengths") 
        }
        
        cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
        
        QList <- QStatsList[SampleNames]
        
        # add sample names to the data frames
        for (i in seq_along(QList)) {
                QList[[i]]$Sample <- names(QList[i])
        }
        
        df <- do.call("rbind", QList)
        ggplot(data = df, aes(x = Cycle, y = Mean_QS, colour = NoReads, group = Sample))  + 
                geom_hline(yintercept = 25, color = 'darkred', linetype = "dashed", size = .35) +
                geom_line() + 
                #geom_point() +
                ylab("Quality Score") + 
                xlab("Cycle") +
                coord_cartesian(xlim=c(xlim_low, xlim_high),
                                ylim = c(8,42), expand = FALSE)  +
                scale_color_gradient2(name = "No reads", limits = c(10000, 60000), midpoint = 35000, low = "#009E73", high = "#E69F00", mid = "#999999") +
                #scale_color_gradient2(name = "No reads", limits = c(min(df$NoReads), max(df$NoReads)), midpoint = min(df$NoReads) + (range(df$NoReads)[2]- range(df$NoReads)[1])/2, low = "#009E73", high = "#E69F00", mid = "#999999")
                ggtitle(paste(Prefix, " reads, Mean_QS, NoSamples", length(SampleNames))) +
                theme_bw() + 
                theme(panel.grid.minor = element_blank(),
                      panel.grid.major.x = element_blank(),
                      panel.grid.major.y = element_line(color = "#999999", size = .15))
        
}


#######################################
#### NoReads_DotPlot
#######################################

NoReads_DotPlot <- function(QStatsList, SampleNames) {
        # Input:
        # QStatsList: list of QStats data frames, such as FW_QualityStats
        # SampleNames: List of Sample names that must be the names of the QStatsList
        
        if (!all(SampleNames %in% names(QStatsList))) {
                stop("not all SampleNames were found in the QStatsList")
        }
        
        cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
        
        QList <- QStatsList[SampleNames]
        # add sample names to the data frames
        for (i in seq_along(QList)) {
                QList[[i]]$Sample <- names(QList[i])
        }
        # combine all df in the list to one big df
        df <- do.call("rbind", QList)
        df <- df[!duplicated(df$Sample), c("Sample", "NoReads")]
        df <- dplyr::arrange(df, desc(NoReads))
        df$Sample <- factor(df$Sample)
        
        # relevel so samples are shown from min NoReads to max NoReads
        LevelsWant <- as.character(df$Sample) 
        for (i in 1:length(LevelsWant)) {
                df$Sample <- relevel(df$Sample, ref = LevelsWant[i])
        }
        
        ggplot(data = df, aes(x = Sample, y = NoReads))  + 
                geom_hline(yintercept = 10000, color = 'darkred', linetype = "dashed", size = .35) +
                geom_point(color = "#E69F00", size = 1) +
                ylab("No Reads") + 
                xlab("") +
                theme_bw() + 
                theme(panel.grid.minor = element_blank(),
                      panel.grid.major.y = element_blank(),
                      panel.grid.major.x = element_line(color = "#999999", size = .15),
                      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6))
}


#######################################
#### NoReads_Histogram
#######################################

NoReads_Histogram <- function(QStatsList, SampleNames) {
        # Input:
        # QStatsList: list of QStats data frames, such as FW_QualityStats
        # SampleNames: List of Sample names that must be the names of the QStatsList
        
        if (!all(SampleNames %in% names(QStatsList))) {
                stop("not all SampleNames were found in the QStatsList")
        }
        
        cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
        
        QList <- QStatsList[SampleNames]
        # add sample names to the data frames
        for (i in seq_along(QList)) {
                QList[[i]]$Sample <- names(QList[i])
        }
        # combine all df in the list to one big df
        df <- do.call("rbind", QList)
        df <- df[!duplicated(df$Sample), c("Sample", "NoReads")]
        
        ggplot(data = df, aes(x = NoReads))  + 
                geom_histogram(binwidth = 1000, col = "black", fill = "#E69F00")+
                geom_vline(xintercept = 10000, color = 'darkred', linetype = "dashed", size = .35) +
                geom_rug() +
                ylab("No Samples") + 
                xlab("No Reads") +
                theme_bw() + 
                theme(panel.grid.minor = element_blank(),
                      panel.grid.major.y = element_blank(),
                      panel.grid.major.x = element_line(color = "#999999", size = .15))
        
}

#######################################
#### NoReads_DotPlot_Type
#######################################

NoReads_DotPlot_Type <- function(QStatsList = NULL, derepFs = NULL, mergers = NULL, mergers.nochim = NULL, SampleNames, sort = TRUE) {
        
        # Input:
        # QStatsList: list of QStats data frames, such as FW_QualityStats
        # SampleNames: List of Sample names that must be names of the QStatsList, derepFs, mergers, and mergers.nochim
        # sort: if the samples should be sorted in the plot based on NoReads
        
        if (all(c(is.null(QStatsList), is.null(derepFs), is.null(mergers), is.null(mergers.nochim)))) {
                stop("QStatsList, derepFs, mergers, mergers.nochim can not all be NULL")
        }
        
        if (!all(SampleNames %in% c(names(QStatsList), names(derepFs), names(mergers), names(mergers.nochim)))) {
                stop("not all SampleNames were found in the given data")
        }
        
        cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
        
        ## Construct the plotting data frame
        
        df.all <- NULL
        
        if (!is.null(QStatsList) & any(SampleNames %in% names(QStatsList))) {
                
                QStatsList <- QStatsList[SampleNames[which(SampleNames %in% names(QStatsList))]]
                
                for (i in seq_along(QStatsList)) {
                        QStatsList[[i]]$Sample <- names(QStatsList[i])
                }
                
                df.all <- do.call("rbind",QStatsList)
                df.all <- df.all[!duplicated(df.all$Sample), c("Sample", "NoReads")]
                df.all$Type <- "all"
                
        }
        
        df.filtered <- NULL
        
        if (!is.null(derepFs) & any(SampleNames %in% names(derepFs))) {
                
                derepFs <- derepFs[SampleNames[which(SampleNames %in% names(derepFs))]]
                
                df.filtered <- data.frame(Sample = names(derepFs), NoReads = 0, Type = "filtered")
                for (i in seq_along(derepFs)) {
                        df.filtered$NoReads[i] <- sum(derepFs[[i]]$uniques) 
                }
                
        }
        
        df.merged <- NULL
        
        if (!is.null(mergers) & any(SampleNames %in% names(mergers))) {
                
                mergers <- mergers[SampleNames[which(SampleNames %in% names(mergers))]]
                
                df.merged <- data.frame(Sample = names(mergers), NoReads = 0, Type = "merged")
                for(i in seq_along(mergers)) {
                        df.merged$NoReads[i] <- sum(mergers[[i]]$abundance)
                }
                
        } 
        
        df.nochim <- NULL
        
        if (!is.null(mergers.nochim) & any(SampleNames %in% names(mergers.nochim))) {
                
                mergers.nochim <- mergers.nochim[SampleNames[which(SampleNames %in% names(mergers.nochim))]]
                
                df.nochim <- data.frame(Sample = names(mergers.nochim), NoReads = 0, Type = "nochim")
                for(i in seq_along(mergers.nochim)) {
                        df.nochim$NoReads[i] <- sum(mergers.nochim[[i]]$abundance)
                }
                
        }
        
        df.plot <- rbind(df.all, df.filtered, df.merged, df.nochim)
        
        if (sort) {
                # sort always by highest level: all, filtered, merged, nochim
                if (!is.null(df.all)) {
                        df.all <- dplyr::arrange(df.all, desc(NoReads))
                        LevelsWant <- as.character(df.all$Sample)
                } else if (!is.null(df.filtered)) {
                        df.filtered <- dplyr::arrange(df.filtered, desc(NoReads))
                        LevelsWant <- as.character(df.filtered$Sample)
                } else if (!is.null(df.merged)) {
                        df.merged <- dplyr::arrange(df.merged, desc(NoReads))
                        LevelsWant <- as.character(df.merged$Sample)
                } else if (!is.null(df.nochim)) {
                        df.nochim <- dplyr::arrange(df.nochim, desc(NoReads))
                        LevelsWant <- as.character(df.nochim$Sample)
                }
                
                
                df.plot$Sample <- factor(df.plot$Sample)
                for (i in seq_along(LevelsWant)) {
                        df.plot$Sample <- relevel(df.plot$Sample, ref = LevelsWant[i])
                }
                
        }
        
        
        ggplot(data = df.plot, aes(x = Sample, y = NoReads, color = Type))  + 
                geom_hline(yintercept = 10000, color = 'darkred', linetype = "dashed", size = .35) +
                geom_point(size = 2) +
                scale_color_manual(values = c(cbPalette[2:4],cbPalette[7])) +
                ylab("No Reads") + 
                xlab("") +
                theme_bw() + 
                theme(panel.grid.minor = element_blank(),
                      panel.grid.major.y = element_blank(),
                      panel.grid.major.x = element_line(color = "#999999", size = .15),
                      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
                      legend.title = element_blank())
        
}



#######################################
#### NoReads_StepsSimple
#######################################

NoReads_StepsSimple <- function(ReadSummary, SampleNames, sort = TRUE) {
        
        
        if (!all(SampleNames %in% ReadSummary$Sample)) {
                stop("not all SampleNames were found in the given ReadSummary")
        }
        
        ReadSummary <- ReadSummary[ReadSummary$Sample %in% SampleNames,]
        
        cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
        
        ## Construct the plotting data frame
        
        ReadSummary <- subset(ReadSummary, select = c("Sample", "NoReads", "FilteredReads", "MergedReads", "MergedReadsWOBimera"))
        
        
        if (sort) {
                
                ReadSummary <- dplyr::arrange(ReadSummary, desc(NoReads))
                LevelsWant <- as.character(ReadSummary$Sample)
                for (i in seq_along(LevelsWant)) {
                        ReadSummary$Sample <- relevel(ReadSummary$Sample, ref = LevelsWant[i])
                }
                
        }
        
        colnames(ReadSummary)[2:5] <- c("all", "filtered", "merged", "nochim")
        ReadSummary <- tidyr::gather(ReadSummary, key = Type, value = NoReads, -Sample)
        
        Tr <- ggplot(data = ReadSummary, aes(x = Sample, y = NoReads, color = Type))  
        Tr <- Tr + 
                geom_hline(yintercept = 10000, color = 'darkred', linetype = "dashed", size = .35) +
                geom_point(size = 2) +
                scale_color_manual(values = c(cbPalette[2:4],cbPalette[7])) +
                ylab("No Reads") + 
                xlab("") +
                theme_bw() + 
                theme(panel.grid.minor = element_blank(),
                      panel.grid.major.y = element_blank(),
                      panel.grid.major.x = element_line(color = "#999999", size = .15),
                      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
                      legend.title = element_blank())
        Tr
        
}



#######################################
#### TotalandUniqueAmplicons
#######################################
# Input:
# FinalNumers: data frame with Sample, UniqueAmplicons and NoAmplicons columns
# sort: default TRUE: sort samples after number total amplicons
# Output:
# the Tr structure
# NB: requires tidyr


TotalandUniqueAmplicons <- function(FinalNumbers, seqtab, sort =TRUE) {
        
        ## show the final amplicon numbers in a plot
        FinalNumbers$TotalAmplicons <- FinalNumbers$NoAmplicons/100 
        
        FinalNumbersL <- tidyr::gather(FinalNumbers, key = Type, value = Number, -Sample, -NoAmplicons)
        
        if(sort) {
                
                FinalNumbersL <- dplyr::arrange(FinalNumbersL, desc(NoAmplicons))
                LevelsWant <- as.character(FinalNumbersL$Sample)
                for (i in seq_along(LevelsWant)) {
                        FinalNumbersL$Sample <- relevel(FinalNumbersL$Sample, ref = LevelsWant[i])
                }
                
        }
        
        Tr <- ggplot(FinalNumbersL, aes(x = Sample, y = Number, col = Type)) +
                geom_point() +
                scale_color_manual(values = c("#E69F00", "#009E73"), labels = c("TotalAmplicons/100", "UniqueAmplicons"))+
                ylab("No Amplicons") + 
                xlab("") +
                #scale_y_continuous(breaks = c(100, 200, 300, 400), labels = c("100\n(10000)", "200\n(20000)", "300\n(30000)", "400\n(40000)")) +
                ggtitle(paste("No unique amplicons in all samples",dim(seqtab)[2]))+
                theme_bw() + 
                theme(panel.grid.minor = element_blank(),
                      panel.grid.major.y = element_blank(),
                      panel.grid.major.x = element_line(color = "#999999", size = .15),
                      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
                      legend.title = element_blank())
        
        return(Tr)
        
}


#######################################
#### plotSVdistributions
#######################################
## REQUIRES dplyr
## Input:
# Seqtab: seqtab from dada2 wrapper
# prevalence: a percentage of samples (bw 0 and 100)
## Output: 
# TrList: List of four Trelis objects, plus distribution data.frame


plotSVdistributions <- function(seqtab, prevalence = 10) {
        
        FinalNumbersSeq <- data.frame(Sequence = colnames(seqtab), InNumberSamples = colSums(seqtab != 0), TotalAmplicons = colSums(seqtab))
        # a look at it shows you that seqtab is ordered after total counts
        FinalNumbersSeq <- group_by(FinalNumbersSeq, InNumberSamples)
        AmpliconDistribution <- dplyr::summarise(FinalNumbersSeq, UniqueAmplicons = n(), TotalAmplicons = sum(TotalAmplicons))
        AmpliconDistribution$CumSumUnique <- rev(cumsum(rev(AmpliconDistribution$UniqueAmplicons)))
        AmpliconDistribution$CumPerCUnique <- rev(cumsum(rev(AmpliconDistribution$UniqueAmplicons/ncol(seqtab))))
        AmpliconDistribution$CumSumTotal <- rev(cumsum(rev(AmpliconDistribution$TotalAmplicons)))
        AmpliconDistribution$CumPerCTotal <- rev(cumsum(rev(AmpliconDistribution$TotalAmplicons/sum(colSums(seqtab)))))
        
        PCValue <- ceiling((prevalence/100)*dim(seqtab)[1]) # tells you in how many samples a SV must be present to meet the prevalence 
        
        # Diff <- AmpliconDistribution$InNumberSamples - PCValue
        # index <- which.max(Diff[Diff<0]) + which.min(Diff[Diff>=0])
        index <- which(AmpliconDistribution$InNumberSamples >= PCValue)[1]
        PCKeptAtPCValue <- AmpliconDistribution$CumPerCTotal[index]
        SVskeptAtPCValue <- AmpliconDistribution$CumPerCUnique[index]
        
        # The number of samples the SVs are present in
        Tr <- ggplot(AmpliconDistribution, aes(x = InNumberSamples, y = UniqueAmplicons))
        Tr <- Tr + geom_point(col = "#E69F00", size = 3) +
                xlab("prevalence") +
                ylab("number of SVs") +
                theme_bw() +
                theme(panel.grid.minor = element_blank(),
                      panel.grid.major.y = element_blank(),
                      panel.grid.major.x = element_line(color = "#999999", size = .15))
        
        
        In1Index <- which(AmpliconDistribution$InNumberSamples == 1)
        if (length(In1Index) != 0) {
                Tr <- Tr + ggtitle(paste(AmpliconDistribution$UniqueAmplicons[In1Index], " of ", AmpliconDistribution$CumSumUnique[1], " SVs (", round(100*AmpliconDistribution$UniqueAmplicons[In1Index]/AmpliconDistribution$CumSumUnique[1], 1), " %)", " were only found in 1 sample", sep = ""))
        } 
        
        
        # Cumulative Percentage of SVs
        Tr1 <- ggplot(AmpliconDistribution, aes(x = InNumberSamples, y = CumPerCUnique))
        Tr1 <- Tr1 + geom_point(col = "#E69F00", size = 3) +
                xlab("prevalence") +
                ylab("cumulative percentage of SVs") +
                theme_bw() +
                theme(panel.grid.minor = element_blank(),
                      panel.grid.major.y = element_blank(),
                      panel.grid.major.x = element_line(color = "#999999", size = .15))
        
        Tr1 <- Tr1 + 
                geom_hline(yintercept = SVskeptAtPCValue, lty =  "dashed") +
                geom_vline(xintercept = (prevalence/100)*dim(seqtab)[1], lty = 'dashed') +
                ggtitle(paste("prevalence ", prevalence, " % = ", round((prevalence/100)*dim(seqtab)[1],1), "; ", AmpliconDistribution$CumSumUnique[index], " of ", AmpliconDistribution$CumSumUnique[1], 
                              " SVs (", round(100*AmpliconDistribution$CumSumUnique[index]/AmpliconDistribution$CumSumUnique[1], 1), " %) have higher prevalence", sep = ""))
        
        
        # The number of samples the total amplicons are present in
        Tr2 <- ggplot(AmpliconDistribution, aes(x = InNumberSamples, y = TotalAmplicons))
        Tr2 <- Tr2 + geom_point(col = "#E69F00", size = 3) +
                xlab("prevalence") +
                ylab("amplicons") +
                theme_bw() +
                theme(panel.grid.minor = element_blank(),
                      panel.grid.major.y = element_blank(),
                      panel.grid.major.x = element_line(color = "#999999", size = .15))
        
        Tr2 <- Tr2 + ggtitle(paste(AmpliconDistribution$CumSumTotal[1] - AmpliconDistribution$CumSumTotal[index], " of ",
                                   AmpliconDistribution$CumSumTotal[1], " (", round(100*(AmpliconDistribution$CumSumTotal[1] - AmpliconDistribution$CumSumTotal[index])/AmpliconDistribution$CumSumTotal[1], 2),
                                    " %) amplicons are from SVs present in less than ", round((prevalence/100)*dim(seqtab)[1],1), " samples.", sep = ""))
        
        Tr3 <- ggplot(AmpliconDistribution, aes(x = InNumberSamples, y = CumPerCTotal))
        Tr3 <- Tr3 + geom_point(col = "#E69F00", size = 3) +
                xlab("prevalence") +
                ylab("cumulative percentage of amplicons") +
                theme_bw() +
                theme(panel.grid.minor = element_blank(),
                      panel.grid.major.y = element_line(color = "#999999", size = .15),
                      panel.grid.major.x = element_line(color = "#999999", size = .15))
        
        Tr3 <- Tr3 + 
                geom_hline(yintercept = PCKeptAtPCValue, lty =  "dashed") +
                geom_vline(xintercept = (prevalence/100)*dim(seqtab)[1], lty = 'dashed') +
                ggtitle(paste("prevalence ", prevalence, " % = ", round((prevalence/100)*dim(seqtab)[1],1), "; ", AmpliconDistribution$CumSumTotal[index], " of ", AmpliconDistribution$CumSumTotal[1], 
                              " amplicons (", round(100*AmpliconDistribution$CumSumTotal[index]/AmpliconDistribution$CumSumTotal[1], 1), " %) would remain", sep = ""))
                
        
        TrList <- list(Tr, Tr1, Tr2, Tr3, AmpliconDistribution)
        
        return(TrList)
        
}

#######################################
#### plotAlphaDiversity
#######################################
# plotAlphaDiversity is related to phyloseq::plot_richness but also offers boxplots (when group != NULL) and
# always renders each plot individually in a list instead of using facet_wraps

plotAlphaDiversity <- function(physeq, measures = NULL, x = "samples", color = NULL, shape = NULL,
                               group = NULL, defCol = "#E69F00"){
        
        
        AlphaDiv <- suppressWarnings(estimate_richness(physeq, measures = measures))
        measures = colnames(AlphaDiv)
        ses = colnames(AlphaDiv)[grep("^se\\.", colnames(AlphaDiv))]
        measures = measures[!measures %in% ses]
        if("Observed" %in% measures){
                
                measures[measures == "Observed"] <- "Richness"
                colnames(AlphaDiv)[colnames(AlphaDiv) == "Observed"] <- "Richness"
        }
        
        if (!is.null(sample_data(physeq, errorIfNULL = FALSE))) {
                DF <- data.frame(AlphaDiv, sample_data(physeq))
        } else {
                DF <- data.frame(AlphaDiv)
        }
        
        DF$samples <- sample_names(physeq)
        
        if (!is.null(x)) {
                if (x %in% c("sample", "samples", "sample_names", "sample.names", "Sample", "Samples")) {
                        x <- "samples"
                }
        } else {
                x <- "samples"
        }
        
        if(!is.null(group)){

                
                TrList <- list()
                
                for (i in 1:length(measures)) {
                        
                        aes_map = aes_string(x = group, y = measures[i], shape = shape, color = color, group = group)
                        
                        if(is.null(color)){
                                Tr = ggplot(DF, aes_map) + geom_boxplot(na.rm = TRUE, col = defCol)
                        } else {
                                
                                Tr = ggplot(DF, aes_map) + geom_boxplot(na.rm = TRUE)
                        }
                        
                        Tr = Tr + geom_jitter(width = .2, height = 0, alpha = 0.65)
                        
                        Tr <- Tr + scale_colour_manual(values = cbPalette[2:8])
                        
                        Tr <- Tr + theme_bw() + xlab("")
                        
                        Tr = Tr + theme(panel.grid.minor = element_blank(),
                                        panel.grid.major.y = element_blank(),
                                        panel.grid.major.x = element_line(color = "#999999", size = .15))
                        
                        TrList[[i]] <- Tr
                        
                }
                
                

        } else {
                
                TrList <- list()
                
                for (i in 1:length(measures)) {
                        
                        aes_map = aes_string(x = x, y = measures[i], shape = shape, color = color)
                        
                        if(is.null(color)){
                                Tr = ggplot(DF, aes_map) + geom_point(na.rm = TRUE, col = defCol)
                        } else {
                                
                                Tr = ggplot(DF, aes_map) + geom_point(na.rm = TRUE)
                        }
                        
                        
                        if (measures[i] == "Chao1") {
                                Tr = Tr + geom_errorbar(aes(ymax = Chao1 + se.chao1, ymin = Chao1 -
                                                                    se.chao1), width = 0.1)
                        }
                        if (measures[i] == "ACE") {
                                Tr = Tr + geom_errorbar(aes(ymax = ACE + se.ACE, ymin = ACE -
                                                                    se.ACE), width = 0.1)
                        }
                        
                        Tr <- Tr + scale_colour_manual(values = cbPalette[2:8])
                        
                        Tr <- Tr + theme_bw()
                        
                        Tr = Tr + theme(panel.grid.minor = element_blank(),
                                        panel.grid.major.y = element_blank(),
                                        panel.grid.major.x = element_line(color = "#999999", size = .15),
                                        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6))
                        
                        TrList[[i]] <- Tr
                        
                }
                
                
                
        }
        
        names(TrList) <- measures
        return(TrList)
        
}


#######################################
#### boxplot_alphaDiv_fromDF
#######################################
# see plotAlphaDiversity

boxplot_alphaDiv_fromDF <- function(DF, measures, color = NULL, shape = NULL,
                                    group = NULL, defCol = "#E69F00") {
        
        
        TrList <- list()
        
        if (!is.null(DF$Total) && !is.null(DF$filtered_reads)) {
                measures <- c(measures, "Total", "filtered_reads")
        }
        
        for (i in 1:length(measures)) {
                
                aes_map = aes_string(x = group, y = measures[i], shape = shape, color = color)
                
                if(is.null(color)){
                        Tr = ggplot(DF, aes_map) + geom_boxplot(na.rm = TRUE, col = defCol)
                } else {
                        
                        Tr = ggplot(DF, aes_map) + geom_boxplot(na.rm = TRUE)
                }
                
                Tr = Tr + geom_jitter(width = .2, height = 0, alpha = 0.65)
                
                Tr <- Tr + scale_colour_manual("", values = cbPalette[2:8])
                
                Tr <- Tr + theme_bw() + xlab("")
                
                Tr = Tr + theme(panel.grid.minor = element_blank(),
                                panel.grid.major.y = element_blank(),
                                panel.grid.major.x = element_line(color = "#999999", size = .15))
                
                TrList[[i]] <- Tr
                
        }
        
        names(TrList) <- measures
        return(TrList)
}



#######################################
#### plot_alphaDivVstotalAmplicons
#######################################
# Function plots the alpha diversity measures from estimate_richness against the total number of reads/amplicons per sample
# and adds a linear fit.
# if color is given the fit is individual on the different colors, but the p-value in the title is still to an overall fit!!

plot_alphaDivVstotalAmplicons <- function(physeq, measures = NULL, color = NULL, shape = NULL, defCol = "#E69F00"){
        
        # needed to get the p-value from-linear fit objects (from stackoverflow)
        lmp <- function (modelobject) {
                if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
                f <- summary(modelobject)$fstatistic
                p <- pf(f[1],f[2],f[3],lower.tail=F)
                attributes(p) <- NULL
                return(p)
        }

        AlphaDiv <- suppressWarnings(estimate_richness(physeq, measures = measures))
        measures = colnames(AlphaDiv)
        ses = colnames(AlphaDiv)[grep("^se\\.", colnames(AlphaDiv))]
        measures = measures[!measures %in% ses]
        
        if("Observed" %in% measures){
                
                measures[measures == "Observed"] <- "Richness"
                colnames(AlphaDiv)[colnames(AlphaDiv) == "Observed"] <- "Richness"
        }
        
        if (!is.null(sample_data(physeq, errorIfNULL = FALSE))) {
                DF <- data.frame(AlphaDiv, sample_data(physeq))
        } else {
                DF <- data.frame(AlphaDiv)
        }
        
        DF$samples <- sample_names(physeq)
        
        if(taxa_are_rows(physeq)){
                DF$TotalReads <- colSums(otu_table(physeq))
                
        } else {
                DF$TotalReads <- rowSums(otu_table(physeq))
        }
        
        
        TrList <- list()
        
        for (i in 1:length(measures)) {
                
                aes_map = aes_string(x = "TotalReads", y = measures[i], shape = shape, color = color)
                
                if(is.null(color)){
                        Tr = ggplot(DF, aes_map) + geom_point(na.rm = TRUE, col = defCol)
                } else {
                        
                        Tr = ggplot(DF, aes_map) + geom_point(na.rm = TRUE)
                }
                
                
                if (measures[i] == "Chao1") {
                        Tr = Tr + geom_errorbar(aes(ymax = Chao1 + se.chao1, ymin = Chao1 -
                                                            se.chao1), width = 0.1)
                }
                if (measures[i] == "ACE") {
                        Tr = Tr + geom_errorbar(aes(ymax = ACE + se.ACE, ymin = ACE -
                                                            se.ACE), width = 0.1)
                }
                
                Tr <- Tr + scale_colour_manual(values = cbPalette[2:8])
                
                # add the regression line and the p-values
                fit <- lm(DF[,measures[i]] ~ DF[,"TotalReads"])
                adjR2 <- round(summary(fit)$adj.r.squared,3)
                pfit <- lmp(fit)
                
                Tr <- Tr + geom_smooth(aes_string(x = "TotalReads", y = measures[i]), method = "lm", se = TRUE, inherit.aes = F)
                
                Tr <- Tr + ggtitle(paste("lm fit: p-value: ", format(pfit, digits = 4), ", adj.r.squared: ", adjR2, sep = ""))
                
                Tr <- Tr + theme_bw() + xlab("total amplicons (taxa_sums())")
                
                Tr = Tr + theme(panel.grid.minor = element_blank())
                
                TrList[[i]] <- Tr
                
        }
  
        names(TrList) <- measures
        return(TrList)
        
}    


#######################################
#### plot_alphaDivVstotalAmplicons
#######################################
# Function plots the alpha diversity measures from estimate_richness against the total number of reads/amplicons per sample
# and adds a linear fit.
# if color is given the fit is individual on the different colors, but the p-value in the title is still to an overall fit!!

plot_alphaDivVstotalAmplicons_fromList <- function(DF_List, measures = NULL, color = NULL, shape = NULL, defCol = "#E69F00"){
        
        TrList <- list()
        
        DF <- DF_List[[1]]
        fitlist <- DF_List[["fitlist"]]
        
        for (i in 1:length(measures)) {
                
                aes_map = aes_string(x = "Total", y = measures[i], shape = shape, color = color)
                
                if(is.null(color)){
                        Tr = ggplot(DF, aes_map) + geom_point(na.rm = TRUE, col = defCol)
                } else {
                        
                        Tr = ggplot(DF, aes_map) + geom_point(na.rm = TRUE)
                }
                
                
                if (measures[i] == "Chao1") {
                        Tr = Tr + geom_errorbar(aes(ymax = Chao1 + se.chao1, ymin = Chao1 -
                                                            se.chao1), width = 0.1)
                }
                if (measures[i] == "ACE") {
                        Tr = Tr + geom_errorbar(aes(ymax = ACE + se.ACE, ymin = ACE -
                                                            se.ACE), width = 0.1)
                }
                
                Tr <- Tr + scale_colour_manual("", values = cbPalette[2:8])
                
                # add the regression line and the p-values
                fit <- fitlist[[measures[i]]]
                adjR2 <- round(summary(fit)$adj.r.squared,3)
                if (adjR2 != 0){
                        pfit <- lmp(fit)
                } else {
                        pfit <- NA
                }
                
                Tr <- Tr + geom_smooth(aes_string(x = "Total", y = measures[i]), method = "lm", se = TRUE, inherit.aes = F)
                
                Tr <- Tr + ggtitle(paste("lm fit: p-value: ", format(pfit, digits = 4), ", adj.r.squared: ", adjR2, sep = ""))
                
                Tr <- Tr + theme_bw() + xlab("total amplicons (sample_sums())")
                
                Tr = Tr + theme(panel.grid.minor = element_blank())
                
                TrList[[i]] <- Tr
                
        }
        
        names(TrList) <- measures
        return(TrList)
        
}

#######################################
#### plot_alphaDivVsfilteredReads_fromList
#######################################
# Function plots the alpha diversity measures from estimate_richness against the given number of filteredReads per sample
# and adds a linear fit.
# if color is given the fit is individual on the different colors, but the p-value in the title is still to an overall fit!!

plot_alphaDivVsfilteredReads_fromList <- function(DF_List, measures = NULL, color = NULL, shape = NULL, defCol = "#E69F00"){
        
        TrList <- list()
        
        DF <- DF_List[[1]]
        fitlist <- DF_List[["fitlist_FilteredReads"]]
        
        for (i in 1:length(measures)) {
                
                aes_map = aes_string(x = "filtered_reads", y = measures[i], shape = shape, color = color)
                
                if(is.null(color)){
                        Tr = ggplot(DF, aes_map) + geom_point(na.rm = TRUE, col = defCol)
                } else {
                        
                        Tr = ggplot(DF, aes_map) + geom_point(na.rm = TRUE)
                }
                
                
                if (measures[i] == "Chao1") {
                        Tr = Tr + geom_errorbar(aes(ymax = Chao1 + se.chao1, ymin = Chao1 -
                                                            se.chao1), width = 0.1)
                }
                if (measures[i] == "ACE") {
                        Tr = Tr + geom_errorbar(aes(ymax = ACE + se.ACE, ymin = ACE -
                                                            se.ACE), width = 0.1)
                }
                
                Tr <- Tr + scale_colour_manual("", values = cbPalette[2:8])
                
                # add the regression line and the p-values
                fit <- fitlist[[measures[i]]]
                adjR2 <- round(summary(fit)$adj.r.squared,3)
                if (adjR2 != 0){
                        pfit <- lmp(fit)
                } else {
                        pfit <- NA
                }
                
                Tr <- Tr + geom_smooth(aes_string(x = "filtered_reads", y = measures[i]), method = "lm", se = TRUE, inherit.aes = F)
                
                Tr <- Tr + ggtitle(paste("lm fit: p-value: ", format(pfit, digits = 4), ", adj.r.squared: ", adjR2, sep = ""))
                
                Tr <- Tr + theme_bw() + xlab("number of filtered reads that entered dada command")
                
                Tr = Tr + theme(panel.grid.minor = element_blank())
                
                TrList[[i]] <- Tr
                
        }
        
        names(TrList) <- measures
        return(TrList)
        
}

#######################################
#### plot_alphaDivVSoriginalTotalAmplicons
#######################################
# wanted to add this plot after rarefaction

plot_alphaDivVSoriginalTotalAmplicons <- function(DF_alpha_rare, DF_alpha_no_rare, measures = NULL, color = NULL, shape = NULL, defCol = "#E69F00"){
        
        measures2 <- measures
        if ("Observed" %in% measures2) {
                measures2[measures2 == "Observed"] <- "Richness" 
        }
        
        DF_alpha_rare$originalTotal <- DF_alpha_no_rare$Total
        DF <- DF_alpha_rare
        
        TrList <- list()
        
        for (i in 1:length(measures2)) {
                
                aes_map = aes_string(x = "originalTotal", y = measures2[i], shape = shape, color = color)
                
                if(is.null(color)){
                        Tr = ggplot(DF, aes_map) + geom_point(na.rm = TRUE, col = defCol)
                } else {
                        
                        Tr = ggplot(DF, aes_map) + geom_point(na.rm = TRUE)
                }
                
                
                if (measures[i] == "Chao1") {
                        Tr = Tr + geom_errorbar(aes(ymax = Chao1 + se.chao1, ymin = Chao1 -
                                                            se.chao1), width = 0.1)
                }
                if (measures[i] == "ACE") {
                        Tr = Tr + geom_errorbar(aes(ymax = ACE + se.ACE, ymin = ACE -
                                                            se.ACE), width = 0.1)
                }
                
                Tr <- Tr + scale_colour_manual("", values = cbPalette[2:8])
                
                # add the regression line and the p-values
                fit <- lm(DF[, measures2[i]] ~ DF[,"originalTotal"])
                adjR2 <- round(summary(fit)$adj.r.squared,3)
                if (adjR2 != 0){
                        pfit <- lmp(fit)
                } else {
                        pfit <- NA
                }
                
                Tr <- Tr + geom_smooth(aes_string(x = "originalTotal", y = measures2[i]), method = "lm", se = TRUE, inherit.aes = F)
                
                Tr <- Tr + ggtitle(paste("lm fit: p-value: ", format(pfit, digits = 4), ", adj.r.squared: ", adjR2, sep = ""))
                
                Tr <- Tr + theme_bw() + xlab("total amplicons before rarefying")
                
                Tr = Tr + theme(panel.grid.minor = element_blank())
                
                TrList[[i]] <- Tr
                
        }
        
        names(TrList) <- measures2
        return(TrList)
        
}



#######################################
#### plot_alphaDivVSoriginalTotalAmplicons2
#######################################
# wanted to add this plot after rarefaction

plot_alphaDivVSoriginalTotalAmplicons2 <- function(DF_alpha_rare, DF_alpha_no_rare, measures = NULL, color = NULL, defCol = "#E69F00"){
        
        measures2 <- measures
        if ("Observed" %in% measures2) {
                measures2[measures2 == "Observed"] <- "Richness" 
        }
        
        DF_alpha_no_rare <- DF_alpha_no_rare[, c("Sample", measures2, "Total", "filtered_reads", color)]
        DF_alpha_no_rare$Type <- "Before rarefaction"
        DF_alpha_rare <- DF_alpha_rare[, c("Sample", measures2, "Total", "filtered_reads", color)]
        to_level <- DF_alpha_rare$Total[1]
        DF_alpha_rare$Type <- paste("After rarefaction to ", to_level, " amplicons", sep = "")
        DF_alpha_rare$Total <- DF_alpha_no_rare$Total
        
        DF <- rbind(DF_alpha_no_rare, DF_alpha_rare)
        DF$Type <- factor(DF$Type, levels = c("Before rarefaction", paste("After rarefaction to ", to_level, " amplicons", sep = "")),
                          ordered = TRUE)
        
        TrList <- list()
        
        for (i in 1:length(measures2)) {
                
                aes_map = aes_string(x = "Total", y = measures2[i], color = color)
                
                if(is.null(color)){
                        Tr = ggplot(DF, aes_map) + geom_point(na.rm = TRUE, col = defCol)
                } else {
                        
                        Tr = ggplot(DF, aes_map) + geom_point(na.rm = TRUE)
                }
                
                
                Tr <- Tr + scale_colour_manual("", values = cbPalette[2:8])
                
                # for the regressin lines and the p-values
                fit <- lm(DF_alpha_no_rare[, measures2[i]] ~ DF_alpha_no_rare[,"Total"])
                adjR2 <- round(summary(fit)$adj.r.squared,3)
                if (adjR2 != 0){
                        pfit <- lmp(fit)
                } else {
                        pfit <- NA
                }
                
                fit_rare <- lm(DF_alpha_rare[, measures2[i]] ~ DF_alpha_rare[,"Total"])
                adjR2_rare <- round(summary(fit_rare)$adj.r.squared,3)
                if (adjR2_rare != 0){
                        pfit_rare <- lmp(fit)
                } else {
                        pfit_rare <- NA
                }
                
                
                
                Tr <- Tr + geom_smooth(aes_string(x = "Total", y = measures2[i]), method = "lm", se = TRUE, inherit.aes = F)
                
                Tr <- Tr + facet_grid(~ Type, scales = "free_y")
                
                Tr <- Tr + ggtitle(paste("fit before: p: ", format(pfit, digits = 4), ", R^2: ", adjR2, "; after: p: ",
                                         format(pfit_rare, digits = 4), ", R^2: ", adjR2_rare, sep = ""))
                
                Tr <- Tr + theme_bw() + xlab("total amplicons before rarefaction")
                
                Tr = Tr + theme(panel.grid.minor = element_blank())
                
                TrList[[i]] <- Tr
                
        }
        
        names(TrList) <- measures2
        
        TrList2 <- list()
        
        for (i in 1:length(measures2)) {
                
                aes_map = aes_string(x = "filtered_reads", y = measures2[i], color = color)
                
                if(is.null(color)){
                        Tr = ggplot(DF, aes_map) + geom_point(na.rm = TRUE, col = defCol)
                } else {
                        
                        Tr = ggplot(DF, aes_map) + geom_point(na.rm = TRUE)
                }
                
                
                Tr <- Tr + scale_colour_manual("", values = cbPalette[2:8])
                
                # for the regressin lines and the p-values
                fit <- lm(DF_alpha_no_rare[, measures2[i]] ~ DF_alpha_no_rare[,"filtered_reads"])
                adjR2 <- round(summary(fit)$adj.r.squared,3)
                if (adjR2 != 0){
                        pfit <- lmp(fit)
                } else {
                        pfit <- NA
                }
                
                fit_rare <- lm(DF_alpha_rare[, measures2[i]] ~ DF_alpha_rare[,"filtered_reads"])
                adjR2_rare <- round(summary(fit_rare)$adj.r.squared,3)
                if (adjR2_rare != 0){
                        pfit_rare <- lmp(fit)
                } else {
                        pfit_rare <- NA
                }
                
                
                
                Tr <- Tr + geom_smooth(aes_string(x = "filtered_reads", y = measures2[i]), method = "lm", se = TRUE, inherit.aes = F)
                
                Tr <- Tr + facet_grid(~ Type, scales = "free_y")
                
                Tr <- Tr + ggtitle(paste("fit before: p: ", format(pfit, digits = 4), ", R^2: ", adjR2, "; after: p: ",
                                         format(pfit_rare, digits = 4), ", R^2: ", adjR2_rare, sep = ""))
                
                Tr <- Tr + theme_bw() + xlab("No of filtered reads (that entered dada algorithm)")
                
                Tr = Tr + theme(panel.grid.minor = element_blank())
                
                TrList2[[i]] <- Tr
                
        }
        
        names(TrList2) <- measures2
        
        TrList3 <- list()
        
        for (i in 1:length(measures2)) {
                
                aes_map = aes_string(x = "Total", y = measures2[i], color = color, shape = "Type")
                
                if(is.null(color)){
                        Tr = ggplot(DF, aes_map) + geom_point(na.rm = TRUE, col = defCol)
                } else {
                        
                        Tr = ggplot(DF, aes_map) + geom_point(na.rm = TRUE)
                }
                
                
                Tr <- Tr + scale_colour_manual("", values = cbPalette[2:8])
                
                
                # for the regressin lines and the p-values
                fit <- lm(DF_alpha_no_rare[, measures2[i]] ~ DF_alpha_no_rare[,"Total"])
                adjR2 <- round(summary(fit)$adj.r.squared,3)
                if (adjR2 != 0){
                        pfit <- lmp(fit)
                } else {
                        pfit <- NA
                }
                
                fit_rare <- lm(DF_alpha_rare[, measures2[i]] ~ DF_alpha_rare[,"Total"])
                adjR2_rare <- round(summary(fit_rare)$adj.r.squared,3)
                if (adjR2_rare != 0){
                        pfit_rare <- lmp(fit)
                } else {
                        pfit_rare <- NA
                }
                
                
                
                Tr <- Tr + geom_smooth(aes_string(x = "Total", y = measures2[i], fill = "Type"), method = "lm", se = FALSE, inherit.aes = F)
                
                
                Tr <- Tr + ggtitle(paste("fit before: p: ", format(pfit, digits = 4), ", R^2: ", adjR2, "; after: p: ",
                                         format(pfit_rare, digits = 4), ", R^2: ", adjR2_rare, sep = ""))
                
                Tr <- Tr + theme_bw() + xlab("total amplicons before rarefaction")
                
                Tr = Tr + theme(panel.grid.minor = element_blank(),
                                legend.title = element_blank(),
                                legend.position = "top")
                
                TrList3[[i]] <- Tr
                
        }
        
        names(TrList3) <- measures2
        
        TrList4 <- list()
        
        for (i in 1:length(measures2)) {
                
                aes_map = aes_string(x = "filtered_reads", y = measures2[i], color = color, shape = "Type")
                
                if(is.null(color)){
                        Tr = ggplot(DF, aes_map) + geom_point(na.rm = TRUE, col = defCol)
                } else {
                        
                        Tr = ggplot(DF, aes_map) + geom_point(na.rm = TRUE)
                }
                
                
                Tr <- Tr + scale_colour_manual("", values = cbPalette[2:8])
                
                
                # for the regressin lines and the p-values
                fit <- lm(DF_alpha_no_rare[, measures2[i]] ~ DF_alpha_no_rare[,"filtered_reads"])
                adjR2 <- round(summary(fit)$adj.r.squared,3)
                if (adjR2 != 0){
                        pfit <- lmp(fit)
                } else {
                        pfit <- NA
                }
                
                fit_rare <- lm(DF_alpha_rare[, measures2[i]] ~ DF_alpha_rare[,"filtered_reads"])
                adjR2_rare <- round(summary(fit_rare)$adj.r.squared,3)
                if (adjR2_rare != 0){
                        pfit_rare <- lmp(fit)
                } else {
                        pfit_rare <- NA
                }
                
                
                
                Tr <- Tr + geom_smooth(aes_string(x = "filtered_reads", y = measures2[i], fill = "Type"), method = "lm", se = FALSE, inherit.aes = F)
                
                
                Tr <- Tr + ggtitle(paste("fit before: p: ", format(pfit, digits = 4), ", R^2: ", adjR2, "; after: p: ",
                                         format(pfit_rare, digits = 4), ", R^2: ", adjR2_rare, sep = ""))
                
                Tr <- Tr + theme_bw() + xlab("No of filtered reads (that entered dada algorithm)")
                
                Tr = Tr + theme(panel.grid.minor = element_blank(),
                                legend.title = element_blank(),
                                legend.position = "top")
                
                TrList4[[i]] <- Tr
                
        }
        
        names(TrList4) <- measures2
        
        list(TrList_total = TrList, TrList_filtered_reads = TrList2, 
             TrList_total_one = TrList3, TrList_filtered_reads_one = TrList4)
}


#######################################
#### plot_correlations_abundance_prev_sparsity
#######################################
# df_ab_prev: data frame with "SV_ID", "total_abundance", "prevalence", "sparsity", "mean_abundance_nonzero",
# "median_abundance_nonzero"
# outputs a list with different plots and fits

plot_correlations_abundance_prev_sparsity <- function(df_ab_prev, col = NULL){ # NB: you could color by Phylum for example
        
        nsamples <- df_ab_prev$prevalence[1] + df_ab_prev$sparsity[1]
        
        Tr_ab <- ggplot(df_ab_prev, aes(x = SV_ID, y = total_abundance))
        if (is.null(col)) {
                Tr_ab <- Tr_ab + geom_point(col = cbPalette[2], size = 2, alpha = 0.7) 
        } else {
                Tr_ab <- Tr_ab + geom_point(aes_string(col = col), size = 2, alpha = 0.7)
        }
        Tr_ab <- Tr_ab +
                ylab("total abundance (taxa_sums())") +
                theme_bw(12)
        
        
        Tr_prev <- ggplot(df_ab_prev, aes(x = SV_ID, y = prevalence))
        if (is.null(col)) {
                Tr_prev <- Tr_prev + geom_point(col = cbPalette[2], size = 2, alpha = 0.7) 
        } else {
                Tr_prev <- Tr_prev + geom_point(aes_string(col = col), size = 2, alpha = 0.7)
        }
        Tr_prev <- Tr_prev +
                theme_bw(12)
        
        
        # - associations of total abundance and log10(total_abundance) to sparsity/prevalence -
        # NB: turned out: log10(total_abundance) is far better correlated with prevalence/sparsity, so I stick with the log fits
        # Also NB: since prevalence = constant - sparsity, a fit prevalence ~ abundance is equal to fit sparsity ~ abundance just with
        # opposite coefficient signs
        
        # fit_spar <- lm(formula = sparsity ~ total_abundance, data = df_ab_prev)
        # pval_spar <- lmp(fit_spar)
        # fit_prev <- lm(formula = prevalence ~ total_abundance, data = df_ab_prev)
        # pval_prev <- lmp(fit_prev)
        fit_prev_log10 <- lm(formula = prevalence ~ log10(total_abundance), data = df_ab_prev)
        pval_prev_log10 <- lmp(fit_prev_log10)
        # fit_spar_log10 <- lm(formula = sparsity ~ log10(total_abundance), data = df_ab_prev)
        # pval_spar_log10 <- lmp(fit_spar_log10)
        # identical(pval_prev_log10, pval_spar_log10) # TRUE
        
        # it comes natural that total abundance and prevalence/sparsity are correlated, but how about mean abundance in non zero samples
        # here I prepare all combinations, but stick to prevalence for now (comment sparsity out) since correlations are the same
        
        
        # fit_spar_mean <- lm(formula = sparsity ~ mean_abundance_nonzero, data = df_ab_prev)
        # pval_spar_mean <- lmp(fit_spar_mean)
        # fit_spar_mean_log10 <- lm(formula = sparsity ~ log10(mean_abundance_nonzero), data = df_ab_prev)
        # pval_spar_mean_log10 <- lmp(fit_spar_mean_log10)
        # fit_spar_median <- lm(formula = sparsity ~ median_abundance_nonzero, data = df_ab_prev)
        # pval_spar_median <- lmp(fit_spar_median)
        # fit_spar_median_log10 <- lm(formula = sparsity ~ log10(median_abundance_nonzero), data = df_ab_prev)
        # pval_spar_median_log10 <- lmp(fit_spar_median_log10)
        
        fit_prev_mean <- lm(formula = prevalence ~ mean_abundance_nonzero, data = df_ab_prev)
        pval_prev_mean <- lmp(fit_prev_mean)
        fit_prev_mean_log10 <- lm(formula = prevalence ~ log10(mean_abundance_nonzero), data = df_ab_prev)
        pval_prev_mean_log10 <- lmp(fit_prev_mean_log10)
        fit_prev_median <- lm(formula = prevalence ~ median_abundance_nonzero, data = df_ab_prev)
        pval_prev_median <- lmp(fit_prev_median)
        fit_prev_median_log10 <- lm(formula = prevalence ~ log10(median_abundance_nonzero), data = df_ab_prev)
        pval_prev_median_log10 <- lmp(fit_prev_median_log10)
        
        
        # Tr_ab_vs_prev <- ggplot(df_ab_prev, aes(x = total_abundance, y = prevalence))
        # Tr_ab_vs_prev <- Tr_ab_vs_prev +
        #         geom_point(col = cbPalette[2], size = 2, alpha = 0.7) +
        #         geom_smooth(method = "lm") +
        #         scale_y_continuous(limits = c(-1, max(df_ab_prev$prevalence) + 5)) +
        #         # annotate("text", x = max(df_ab_prev$abundance)/5, y = 5, label = paste("R.square of lm: ", as.character(round(summary(fit_prev)$r.squared,4), sep = ""))) +
        #         xlab("total abundance (taxa_sums())") +
        #         ggtitle(paste("lm fit: p.val: ", format(pval_prev, digits = 4), "R.square: ", as.character(round(summary(fit_prev)$r.squared,4), sep = ""))) +
        # theme_bw(12)
        # Tr_ab_vs_prev_75Q <- Tr_ab_vs_prev + coord_cartesian(xlim = c(-5, quantile(df_ab_prev$abundance, .75)))
        
        Tr_prev_vs_log10_ab <- ggplot(df_ab_prev, aes(x = total_abundance, y = prevalence))
        if (is.null(col)) {
                Tr_prev_vs_log10_ab <- Tr_prev_vs_log10_ab + geom_point(col = cbPalette[2], size = 2, alpha = 0.7) 
        } else {
                Tr_prev_vs_log10_ab <- Tr_prev_vs_log10_ab + geom_point(aes_string(col = col), size = 2, alpha = 0.7)
        }
        Tr_prev_vs_log10_ab <- Tr_prev_vs_log10_ab +
                scale_x_log10() +
                scale_y_continuous(limits = c(-1, max(df_ab_prev$prevalence) + 5)) +
                geom_smooth(method = "lm") +
                # annotate("text", x = max(df_ab_prev$abundance)/5, y = 5, label = paste("R.square of lm: ", as.character(round(summary(fit_prev)$r.squared,4), sep = ""))) +
                xlab("total abundance (taxa_sums())") +
                ggtitle(paste("lm fit: p.val: ", format(pval_prev_log10, digits = 4), "R.square: ", as.character(round(summary(fit_prev_log10)$r.squared,4), sep = ""))) +
                theme_bw(12)
        
        # NB: you could color or facet by phylum
        
        # Tr_prev_vs_log10_ab <- ggplot(df_ab_prev, aes(x = total_abundance, y = prevalence))
        # Tr_prev_vs_log10_ab <- Tr_prev_vs_log10_ab +
        #         geom_point(aes(col = Phylum), size = 2, alpha = 0.7) +
        #         scale_x_log10() +
        #         geom_smooth(method = "lm") +
        #         # annotate("text", x = max(df_ab_prev$abundance)/5, y = 5, label = paste("R.square of lm: ", as.character(round(summary(fit_prev)$r.squared,4), sep = ""))) +
        #         xlab("total abundance (taxa_sums())") +
        #         ggtitle(paste("lm fit: p.val: ", format(pval_prev_log10, digits = 4), "R.square: ", as.character(round(summary(fit_prev_log10)$r.squared,4), sep = ""))) +
        #         theme_bw(12)
        
        
        # Tr_spar_vs_meanab <- ggplot(df_ab_prev, aes(x = mean_abundance_nonzero, y = sparsity))
        # Tr_spar_vs_meanab <- Tr_spar_vs_meanab +
        #         geom_point(col = cbPalette[2], size = 2, alpha = 0.7) +
        #         geom_smooth(method = "lm") +
        #         scale_y_continuous(limits = c(-1, nsamples + 5)) +
        #         xlab("SVs mean abundance in non-zero samples") +
        #         ggtitle(paste("lm fit: p.val: ", format(pval_spar_mean, digits = 4), "R.square: ", as.character(round(summary(fit_spar_mean)$r.squared,4), sep = ""))) +
        #         theme_bw(12)
        # Tr_spar_vs_meanab_75Q <- Tr_spar_vs_meanab + coord_cartesian(xlim = c(-5, quantile(df_ab_prev$mean_abundance_nonzero, .75)))
        # 
        
        # Tr_spar_vs_log10_meanab <- ggplot(df_ab_prev, aes(x = mean_abundance_nonzero, y = sparsity))
        # Tr_spar_vs_log10_meanab <- Tr_spar_vs_log10_meanab +
        #         geom_point(col = cbPalette[2], size = 2, alpha = 0.7) +
        #         scale_x_log10() +
        #         geom_smooth(method = "lm") +
        #         xlab("SVs mean abundance in non-zero samples") +
        #         ggtitle(paste("lm fit: p.val: ", format(pval_spar_mean_log10, digits = 4), "R.square: ", as.character(round(summary(fit_spar_mean_log10)$r.squared,4), sep = ""))) +
        #         theme_bw(12)
        
        Tr_prev_vs_log10_meanab <- ggplot(df_ab_prev, aes(x = mean_abundance_nonzero, y = prevalence))
        if (is.null(col)) {
                Tr_prev_vs_log10_meanab <- Tr_prev_vs_log10_meanab + geom_point(col = cbPalette[2], size = 2, alpha = 0.7) 
        } else {
                Tr_prev_vs_log10_meanab <- Tr_prev_vs_log10_meanab + geom_point(aes_string(col = col), size = 2, alpha = 0.7)
        }
        Tr_prev_vs_log10_meanab <- Tr_prev_vs_log10_meanab +
                scale_x_log10() +
                geom_smooth(method = "lm") +
                scale_y_continuous(limits = c(-1, max(df_ab_prev$prevalence) + 5)) +
                xlab("mean abundance of SV in non-zero samples") +
                ggtitle(paste("lm fit: p.val: ", format(pval_prev_mean_log10, digits = 4), "R.square: ", as.character(round(summary(fit_prev_mean_log10)$r.squared,4), sep = ""))) +
                theme_bw(12)
        
        Tr_prev_vs_log10_medianab <- ggplot(df_ab_prev, aes(x = median_abundance_nonzero, y = prevalence))
        if (is.null(col)) {
                Tr_prev_vs_log10_medianab <- Tr_prev_vs_log10_medianab + geom_point(col = cbPalette[2], size = 2, alpha = 0.7) 
        } else {
                Tr_prev_vs_log10_medianab <- Tr_prev_vs_log10_medianab + geom_point(aes_string(col = col), size = 2, alpha = 0.7)
        }
        Tr_prev_vs_log10_medianab <- Tr_prev_vs_log10_medianab +
                scale_x_log10() +
                geom_smooth(method = "lm") +
                scale_y_continuous(limits = c(-1, max(df_ab_prev$prevalence) + 5)) +
                xlab("median abundance of SV in non-zero samples") +
                ggtitle(paste("lm fit: p.val: ", format(pval_prev_median_log10, digits = 4), "R.square: ", as.character(round(summary(fit_prev_median_log10)$r.squared,4), sep = ""))) +
                theme_bw(12)
        
  
        
        
        fitlist <- list(fit_prev_log10 = fit_prev_log10, fit_prev_mean = fit_prev_mean, 
                        fit_prev_mean_log10 = fit_prev_mean_log10, fit_prev_median = fit_prev_median,
                        fit_prev_median_log10 = fit_prev_median_log10)
        
        out <- list(Tr_ab = Tr_ab, 
                    Tr_prev = Tr_prev, 
                    Tr_prev_vs_log10_ab = Tr_prev_vs_log10_ab, 
                    Tr_prev_vs_log10_meanab = Tr_prev_vs_log10_meanab,
                    Tr_prev_vs_log10_medianab = Tr_prev_vs_log10_medianab,
                    fitlist = fitlist)
        
}



#######################################
#### plot_abundance_prev_filter
#######################################
# df_ab_prev: data frame with "SV_ID", "total_abundance", "prevalence", "sparsity", "mean_abundance_nonzero",
# "median_abundance_nonzero" plus tax_table

plot_abundance_prev_filter <- function(physeq, prevalence, taxa_sums_quantile){ 
        
        df_ab_prev <- data.frame(SV_ID = 1:ntaxa(physeq), 
                                 total_abundance = taxa_sums(physeq),
                                 prevalence = colSums(as(otu_table(physeq), "matrix") != 0),
                                 sparsity = colSums(as(otu_table(physeq), "matrix") == 0), 
                                 mean_abundance_nonzero = apply(as(otu_table(physeq), "matrix"), 2, function(x){mean(x[x > 0])}),
                                 median_abundance_nonzero = apply(as(otu_table(physeq), "matrix"), 2, function(x){median(x[x > 0])}))
        
        df_ab_prev <- cbind(df_ab_prev, tax_table(physeq))
        
        prev_thresh <- (prevalence/100)*nsamples(physeq)
        abund_thresh <- quantile(taxa_sums(physeq), probs = taxa_sums_quantile/100)
        
        df_ab_prev_filt <- dplyr::filter(df_ab_prev, prevalence > prev_thresh | total_abundance > abund_thresh)
        
        no_samples <- df_ab_prev$prevalence[1] + df_ab_prev$sparsity[1]
        shade_df <- data.frame(total_abundance = 0, prevalence = 0)
        
        
        Tr_prev_vs_log10_ab <- ggplot(df_ab_prev, aes(x = total_abundance, y = prevalence))
        Tr_prev_vs_log10_ab <- Tr_prev_vs_log10_ab +
                geom_point(col = cbPalette[2], size = 2, alpha = 0.7) +
                scale_x_log10() +
                geom_rect(data = shade_df, xmin = -Inf, xmax = log10(abund_thresh), ymin = -Inf, ymax = prev_thresh, fill = "#660033", alpha = 0.4) +
                geom_hline(yintercept = (prevalence/100)*nsamples(physeq), col = cbPalette[1], lty = "dashed") +
                geom_vline(xintercept = quantile(taxa_sums(physeq), probs = taxa_sums_quantile/100), col = cbPalette[1], lty = "dashed") +
                xlab("total abundance (taxa_sums())") + 
                theme_bw(12) +
                ggtitle(paste(nrow(df_ab_prev_filt), " of ", nrow(df_ab_prev), " SVs (", round(100*nrow(df_ab_prev_filt)/nrow(df_ab_prev), 1),
                              " %) and ", round(sum(df_ab_prev_filt$total_abundance)), " of ", round(sum(df_ab_prev$total_abundance)), " amplicons (",
                              round((sum(df_ab_prev_filt$total_abundance)/sum(df_ab_prev$total_abundance))*100, 1), " %) remain", sep = ""))
        
        
        Tr_prev_vs_log10_ab_col <- ggplot(df_ab_prev, aes(x = total_abundance, y = prevalence))
        Tr_prev_vs_log10_ab_col <- Tr_prev_vs_log10_ab_col +
                geom_point(aes(col = Phylum), size = 2, alpha = 0.7) +
                scale_x_log10() +
                geom_rect(data = shade_df, xmin = -Inf, xmax = log10(abund_thresh), ymin = -Inf, ymax = prev_thresh, fill = "#660033", alpha = 0.4) +
                geom_hline(yintercept = (prevalence/100)*nsamples(physeq), col = cbPalette[1], lty = "dashed") +
                geom_vline(xintercept = quantile(taxa_sums(physeq), probs = taxa_sums_quantile/100), col = cbPalette[1], lty = "dashed") +
                xlab("total abundance (taxa_sums())") +
                facet_wrap(~Phylum) +
                theme_bw(12) +
                theme(legend.position = "none")
        
        
        # phylum_df <- df_ab_prev[, c("Phylum", "total_abundance", "prevalence")]
        # phylum_df <- group_by(phylum_df, Phylum)
        # phylum_df <- dplyr::summarise(phylum_df, SVs = n(), abundance = round(sum(total_abundance)))
        # phylum_df_filt <- df_ab_prev_filt[, c("Phylum", "total_abundance", "prevalence")]
        # phylum_df_filt <- group_by(phylum_df_filt, Phylum)
        # phylum_df_filt <- dplyr::summarise(phylum_df_filt, SVs = n(), abundance = round(sum(total_abundance)))
        # phylum_df_summary <- merge(phylum_df, phylum_df_filt, by = "Phylum")
        # colnames(phylum_df_summary) <- c("Phylum", "SVs_before", "abundance_before", "SVs_after", 'abundance_after')
        # phylum_df_summary <- mutate(phylum_df_summary, SV_r_PC = round(100*SVs_after/SVs_before, 1), abundance_r_PC = round(100*abundance_after/abundance_before, 1),
        #                             SV_PC = round(100*SVs_after/sum(SVs_after), 1), abundance_PC = round(100*abundance_after/sum(abundance_after), 1))
        
        Before <- summarise(group_by(df_ab_prev, Phylum), SVs_bef = n(), PC_SV_bef = round(100*n()/nrow(df_ab_prev),1), PC_total_pre_bef = round(100*sum(prevalence)/sum(df_ab_prev$prevalence), 1),
                           PC_total_ab_bef = round(100*sum(total_abundance)/sum(df_ab_prev$total_abundance), 1), 
                           mean_pre_bef = round(100*mean(prevalence)/no_samples, 1),
                           mean_tot_ab_bef = round(mean(total_abundance)), 
                           mean_mean_ab_nonzero_bef = round(mean(mean_abundance_nonzero)),
                           med_med_ab_nonzero_bef = round(median(median_abundance_nonzero)),
                           total_ab_bef = sum(total_abundance))
        
        After <- summarise(group_by(df_ab_prev_filt, Phylum), SVs_aft = n(), PC_SV_aft = round(100*n()/nrow(df_ab_prev_filt),1), PC_total_pre_aft = round(100*sum(prevalence)/sum(df_ab_prev_filt$prevalence), 1),
                           PC_total_ab_aft = round(100*sum(total_abundance)/sum(df_ab_prev_filt$total_abundance), 1), 
                           mean_pre_aft = round(100*mean(prevalence)/no_samples, 1),
                           mean_tot_ab_aft = round(mean(total_abundance)),
                           mean_mean_ab_nonzero_aft = round(mean(mean_abundance_nonzero)),
                           med_med_ab_nonzero_aft = round(median(median_abundance_nonzero)),
                           total_ab_aft = sum(total_abundance))
        
        Before_total <- summarise(df_ab_prev, SVs_bef = n(), PC_SV_bef = round(100*n()/nrow(df_ab_prev),1), PC_total_pre_bef = round(100*sum(prevalence)/sum(df_ab_prev$prevalence), 1),
                            PC_total_ab_bef = round(100*sum(total_abundance)/sum(df_ab_prev$total_abundance), 1), 
                            mean_pre_bef = round(100*mean(prevalence)/no_samples, 1),
                            mean_tot_ab_bef = round(mean(total_abundance)), 
                            mean_mean_ab_nonzero_bef = round(mean(mean_abundance_nonzero)),
                            med_med_ab_nonzero_bef = round(median(median_abundance_nonzero)),
                            total_ab_bef = sum(total_abundance))
        
        Before_total <- data.frame(Phylum = "Total", Before_total)
        
        After_total <- summarise(df_ab_prev_filt, SVs_aft = n(), PC_SV_aft = round(100*n()/nrow(df_ab_prev_filt),1), PC_total_pre_aft = round(100*sum(prevalence)/sum(df_ab_prev_filt$prevalence), 1),
                           PC_total_ab_aft = round(100*sum(total_abundance)/sum(df_ab_prev_filt$total_abundance), 1), 
                           mean_pre_aft = round(100*mean(prevalence)/no_samples, 1),
                           mean_tot_ab_aft = round(mean(total_abundance)),
                           mean_mean_ab_nonzero_aft = round(mean(mean_abundance_nonzero)),
                           med_med_ab_nonzero_aft = round(median(median_abundance_nonzero)),
                           total_ab_aft = sum(total_abundance))
        
        After_total <- data.frame(Phylum = "Total", After_total)
        
        Before <- rbind(Before, Before_total)
        
        After <- rbind(After, After_total)
        
        
        Merged <- merge(Before, After, by = "Phylum", all = TRUE, sort = FALSE)
        
        Merged <- mutate(Merged, PC_SV_rem = round(100*SVs_aft/SVs_bef, 1), PC_ab_rem = round(100*total_ab_aft/total_ab_bef, 1))
        
        Merged <- select(Merged, 1, 20, 21, 2, 11, 3, 12, 4, 13, 5, 14, 6, 15, 7, 16, 8, 17, 9, 18)
        
        Merged$tot_ab_bef <- round(Before$total_ab_bef)
        Merged$tot_ab_aft <- round(After$total_ab_aft)
        
        out <- list(Tr_prev_vs_log10_ab = Tr_prev_vs_log10_ab,
                    Tr_prev_vs_log10_ab_col = Tr_prev_vs_log10_ab_col,
                    phylum_df_summary = Merged)
        
}


###################### OLD FUNCTIONS #############################

#######################################
#### NoReads_Steps
#######################################

NoReads_Steps <- function(QStatsList = NULL, NoFilteredReads = NULL, mergers = NULL, mergers.nochim = NULL, SampleNames, sort = TRUE) {
        
        # Input:
        # QStatsList: list of QStats data frames, such as FW_QualityStats
        # SampleNames: List of Sample names that must be names of the QStatsList, NoFilteredReads, mergers, and mergers.nochim
        # sort: if the samples should be sorted in the plot based on NoReads
        
        if (all(c(is.null(QStatsList), is.null(NoFilteredReads), is.null(mergers), is.null(mergers.nochim)))) {
                stop("QStatsList, derepFs, mergers, mergers.nochim can not all be NULL")
        }
        
        if (!all(SampleNames %in% c(names(QStatsList), names(NoFilteredReads), names(mergers), names(mergers.nochim)))) {
                stop("not all SampleNames were found in the given data")
        }
        
        cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
        
        ## Construct the plotting data frame
        
        df.all <- NULL
        
        if (!is.null(QStatsList) & any(SampleNames %in% names(QStatsList))) {
                
                QStatsList <- QStatsList[SampleNames[which(SampleNames %in% names(QStatsList))]]
                
                for (i in seq_along(QStatsList)) {
                        QStatsList[[i]]$Sample <- names(QStatsList[i])
                }
                
                df.all <- do.call("rbind",QStatsList)
                df.all <- df.all[!duplicated(df.all$Sample), c("Sample", "NoReads")]
                df.all$Type <- "all"
                
        }
        
        df.filtered <- NULL
        
        if (!is.null(NoFilteredReads) & any(SampleNames %in% names(NoFilteredReads))) {
                
                NoFilteredReads <- NoFilteredReads[SampleNames[which(SampleNames %in% names(NoFilteredReads))]]
                
                df.filtered <- data.frame(Sample = names(NoFilteredReads), NoReads = NoFilteredReads, Type = "filtered")
                
        }
        
        df.merged <- NULL
        
        if (!is.null(mergers) & any(SampleNames %in% names(mergers))) {
                
                mergers <- mergers[SampleNames[which(SampleNames %in% names(mergers))]]
                
                df.merged <- data.frame(Sample = names(mergers), NoReads = 0, Type = "merged")
                for(i in seq_along(mergers)) {
                        df.merged$NoReads[i] <- sum(mergers[[i]]$abundance)
                }
                
        } 
        
        df.nochim <- NULL
        
        if (!is.null(mergers.nochim) & any(SampleNames %in% names(mergers.nochim))) {
                
                mergers.nochim <- mergers.nochim[SampleNames[which(SampleNames %in% names(mergers.nochim))]]
                
                df.nochim <- data.frame(Sample = names(mergers.nochim), NoReads = 0, Type = "nochim")
                for(i in seq_along(mergers.nochim)) {
                        df.nochim$NoReads[i] <- sum(mergers.nochim[[i]]$abundance)
                }
                
        }
        
        df.plot <- rbind(df.all, df.filtered, df.merged, df.nochim)
        
        if (sort) {
                # sort always by highest level: all, filtered, merged, nochim
                if (!is.null(df.all)) {
                        df.all <- dplyr::arrange(df.all, desc(NoReads))
                        LevelsWant <- as.character(df.all$Sample)
                } else if (!is.null(df.filtered)) {
                        df.filtered <- dplyr::arrange(df.filtered, desc(NoReads))
                        LevelsWant <- as.character(df.filtered$Sample)
                } else if (!is.null(df.merged)) {
                        df.merged <- dplyr::arrange(df.merged, desc(NoReads))
                        LevelsWant <- as.character(df.merged$Sample)
                } else if (!is.null(df.nochim)) {
                        df.nochim <- dplyr::arrange(df.nochim, desc(NoReads))
                        LevelsWant <- as.character(df.nochim$Sample)
                }
                
                
                df.plot$Sample <- factor(df.plot$Sample)
                for (i in seq_along(LevelsWant)) {
                        df.plot$Sample <- relevel(df.plot$Sample, ref = LevelsWant[i])
                }
                
        }
        
        
        ggplot(data = df.plot, aes(x = Sample, y = NoReads, color = Type))  + 
                geom_hline(yintercept = 10000, color = 'darkred', linetype = "dashed", size = .35) +
                geom_point(size = 2) +
                scale_color_manual(values = c(cbPalette[2:4],cbPalette[7])) +
                ylab("No Reads") + 
                xlab("") +
                theme_bw() + 
                theme(panel.grid.minor = element_blank(),
                      panel.grid.major.y = element_blank(),
                      panel.grid.major.x = element_line(color = "#999999", size = .15),
                      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
                      legend.title = element_blank())
        
}




