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
#### AmpliconDistribution
#######################################
## REQUIRES dplyr
## Input:
# Seqtab: seqtab from dada2 wrapper
# PCentage: illustrates in how many amplicons you would keep if only considering amplicons that are present in at least PCentage of the samples
## Output: 
# TrList: List of four Trelis objects, plus AmpliconDistribution data.frame


AmpliconDistribution <- function(seqtab, PCentage = 10) {
        
        FinalNumbersSeq <- data.frame(Sequence = colnames(seqtab), InNumberSamples = colSums(seqtab != 0), TotalAmplicons = colSums(seqtab))
        FinalNumbersSeq <- group_by(FinalNumbersSeq, InNumberSamples)
        AmpliconDistribution <- dplyr::summarise(FinalNumbersSeq, UniqueAmplicons = n(), TotalAmplicons = sum(TotalAmplicons))
        AmpliconDistribution$CumSumUnique <- rev(cumsum(rev(AmpliconDistribution$UniqueAmplicons)))
        AmpliconDistribution$CumPerCUnique <- rev(cumsum(rev(AmpliconDistribution$UniqueAmplicons/ncol(seqtab))))
        AmpliconDistribution$CumSumTotal <- rev(cumsum(rev(AmpliconDistribution$TotalAmplicons)))
        AmpliconDistribution$CumPerCTotal <- rev(cumsum(rev(AmpliconDistribution$TotalAmplicons/sum(colSums(seqtab)))))
        
        PCValue <- ceiling((PCentage/100)*dim(seqtab)[1])
        Diff <- AmpliconDistribution$InNumberSamples - PCValue
        index <- which.max(Diff[Diff<0]) + which.min(Diff[Diff>=0])
        PCKeptAtPCValue <- AmpliconDistribution$CumPerCTotal[index]
        
        # The number of samples the unique amplicons are present in
        Tr <- ggplot(AmpliconDistribution, aes(x = InNumberSamples, y = UniqueAmplicons))
        Tr <- Tr + geom_point(col = "#E69F00", size = 3) +
                xlab("Present in Number of Samples") +
                ylab("Unique Amplicons") +
                theme_bw() +
                theme(panel.grid.minor = element_blank(),
                      panel.grid.major.y = element_blank(),
                      panel.grid.major.x = element_line(color = "#999999", size = .15))
        
        # Cumulative Percentage of Unique Amplicons
        Tr1 <- ggplot(AmpliconDistribution, aes(x = InNumberSamples, y = CumPerCUnique))
        Tr1 <- Tr1 + geom_point(col = "#E69F00", size = 3) +
                xlab("Present in Number of Samples") +
                ylab("Cumulative Percentage of unique Amplicons") +
                theme_bw() +
                theme(panel.grid.minor = element_blank(),
                      panel.grid.major.y = element_blank(),
                      panel.grid.major.x = element_line(color = "#999999", size = .15))
        
        # The number of samples the total amplicons are present in
        Tr2 <- ggplot(AmpliconDistribution, aes(x = InNumberSamples, y = TotalAmplicons))
        Tr2 <- Tr2 + geom_point(col = "#E69F00", size = 3) +
                xlab("Present in Number of Samples") +
                ylab("Total Amplicons") +
                theme_bw() +
                theme(panel.grid.minor = element_blank(),
                      panel.grid.major.y = element_blank(),
                      panel.grid.major.x = element_line(color = "#999999", size = .15))
        
        Tr3 <- ggplot(AmpliconDistribution, aes(x = InNumberSamples, y = CumPerCTotal))
        Tr3 <- Tr3 + geom_point(col = "#E69F00", size = 3) +
                xlab("Present in Number of Samples") +
                ylab("Cumulative Percentage of total Amplicons") +
                geom_hline(yintercept = PCKeptAtPCValue, lty =  "dashed") +
                geom_vline(xintercept = PCValue, lty = 'dashed') +
                ggtitle(paste("Total Amplicons", AmpliconDistribution$CumSumTotal[1], ";", PCentage, "%:", index, "samples;", round(PCKeptAtPCValue, 3), "% of Total Amplicons remain")) +
                theme_bw() +
                theme(panel.grid.minor = element_blank(),
                      panel.grid.major.y = element_line(color = "#999999", size = .15),
                      panel.grid.major.x = element_line(color = "#999999", size = .15))
        
        TrList <- list(Tr, Tr1, Tr2, Tr3, AmpliconDistribution)
        
        return(TrList)
        
}

#######################################
#### plotAlphaDiversity
#######################################
# Function related to plot_richness from phyloseq but plots each measure individually and outputs a Trlist

plotAlphaDiversity <- function(physeq, measures = NULL, x = "samples", color = NULL, shape = NULL,
                               group = NULL, defCol = "#E69F00"){
        
        cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
        
        AlphaDiv <- estimate_richness(physeq, measures = measures)
        measures = colnames(AlphaDiv)
        ses = colnames(AlphaDiv)[grep("^se\\.", colnames(AlphaDiv))]
        measures = measures[!measures %in% ses]
        
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

                # # OLD VERSION COmmented out
                # DF <- dplyr::group_by_(DF, group)
                # 
                # DF1 <- dplyr::summarise_each_(DF, funs(mean, sd), measures)
                # 
                # measures1 <- paste(measures, "_mean", sep = "")
                # measures2 <- paste(measures, "_sd", sep = "")
                # 
                # TrList <- list()
                # 
                # for (i in 1:length(measures)) {
                #         
                #         aes_map = aes_string(x = group, y = measures1[i], shape = shape, color = color)
                #         
                #         if(is.null(color)){
                #                 Tr = ggplot(DF1, aes_map) + geom_point(na.rm = TRUE, col = defCol, size = 3)
                #         } else {
                #                 
                #                 Tr = ggplot(DF1, aes_map) + geom_point(na.rm = TRUE, size = 3)
                #         }
                #         
                #         Tr = Tr + geom_errorbar(aes_string(ymax = paste(measures1[i], "+", measures2[i], sep =""), ymin = paste(measures1[i], "-", measures2[i], sep ="")), width = 0.1)
                #         
                #         Tr <- Tr + scale_colour_manual(values = cbPalette[2:8])
                #         
                #         Tr = Tr + theme(panel.grid.minor = element_blank(),
                #                         panel.grid.major.y = element_blank(),
                #                         panel.grid.major.x = element_line(color = "#999999", size = .15))
                #         
                #         TrList[[i]] <- Tr
                #         
                # }
                
                TrList <- list()
                
                for (i in 1:length(measures)) {
                        
                        aes_map = aes_string(x = group, y = measures[i], shape = shape, color = color, group = group)
                        
                        if(is.null(color)){
                                Tr = ggplot(DF, aes_map) + geom_boxplot(na.rm = TRUE, col = defCol)
                        } else {
                                
                                Tr = ggplot(DF, aes_map) + geom_boxplot(na.rm = TRUE)
                        }
                        
                        Tr = Tr + geom_jitter(width = .2, alpha = 0.65)
                        
                        Tr <- Tr + scale_colour_manual(values = cbPalette[2:8])
                        
                        Tr <- Tr + theme_bw()
                        
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
#### plotAlphaDivVsSeqDepth
#######################################
# Function plots the alpha diversity measures from estimate_richness against the total number of reads/amplicons per sample
# and adds a linear fit.
# if color is given the fit is individual on the different colors, but the p-value in the title is still to an overall fit!!

plotAlphaDivVsSeqDepth <- function(physeq, measures = NULL, color = NULL, shape = NULL, defCol = "#E69F00"){
        
        # needed to get the p-value from-linear fit objects (from stackoverflow)
        lmp <- function (modelobject) {
                if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
                f <- summary(modelobject)$fstatistic
                p <- pf(f[1],f[2],f[3],lower.tail=F)
                attributes(p) <- NULL
                return(p)
        }
        
        cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
        
        AlphaDiv <- estimate_richness(physeq, measures = measures)
        measures = colnames(AlphaDiv)
        ses = colnames(AlphaDiv)[grep("^se\\.", colnames(AlphaDiv))]
        measures = measures[!measures %in% ses]
        
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
                
                Tr <- Tr + geom_smooth(method = "lm", se = TRUE)
                
                Tr <- Tr + ggtitle(paste("Linear Fit: p-value: ", format(pfit, digits = 4), ", adj.r.squared: ", adjR2, sep = ""))
                
                Tr <- Tr + theme_bw()
                
                Tr = Tr + theme(panel.grid.minor = element_blank())
                
                TrList[[i]] <- Tr
                
        }
  
        names(TrList) <- measures
        return(TrList)
        
}    

