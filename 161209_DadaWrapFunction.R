#path <- "/Users/jvb740/MarieCurie_Work/BackgroundKnowledge/16S_Learning/DanFunD"
#F_pattern <- "*F.fastq.gz"
#R_pattern <- "*R.fastq.gz"

## Input
# path: The path to the folder with the fastq files: NB: currently only the version
# where the path folder contains folders (named after the samples) that contain the fastq files
# Further versions such as fastq files directly in the path folder should be added when
# necessary

# NB: the plots are currently assuming 250 nt

Dada2_wrap <- function(path, F_pattern, R_pattern, path2 = NULL,
                       trimLeft = c(10,10), truncLen = c(220, 160), 
                       maxN = 0, maxEE = 2, truncQ = 2,
                       NSAM.LEARN = NULL,
                       err_F = NULL,
                       err_R = NULL,
                       minOverlap = 20,
                       maxMismatch = 0) {
        
        ptm <- proc.time()
        
        ##############################
        ### call the required packages
        ##############################
        
        ## dada2:
        # source("https://bioconductor.org/biocLite.R")
        try(library(dada2), biocLite("dada2"))
        
        ## Short Read
        try(library(ShortRead), biocLite("ShortRead"))
        
        ## ggplot2
        try(library(ggplot2), install.packages("ggplot2"))
        
        ## dplyr
        try(library(dplyr), install.packages("dplyr"))
        
        ## dplyr
        try(library(tidyr), install.packages("tidyr"))
        
        
        message("*********************** All packages loaded ***********************
                ********************************************************************")
        
        
        ## tidyr
        # library(tidyr)
        ##############################
        
        ##############################
        ### save the package Versions
        ##############################
        # NB: outputs an error and stops function if a Package is not installed
        PackageVersions <- data.frame(Package = c("dada2", "ShortRead", "ggplot2", "dplyr"),
                                      Version = c(packageVersion("dada2"),
                                                  packageVersion("ShortRead"),
                                                  packageVersion("ggplot2"),
                                                  packageVersion("dplyr")))
        
        message(paste(PackageVersions$Package, ": ", PackageVersions$Version, "; ", sep= ""))
        ##############################
        
        
        
        
        ##############################
        ### Construct character vectors to the FW and RV fastq files
        ##############################
        if(is.null(path2)){
                
                path2 = path
        }
        
        
        folders <- list.dirs(path, recursive = FALSE, full.names = FALSE) 
        
        if(length(folders != 0)) {
                
                SampleNames <- folders
                
                if(sum(grepl("Dada", SampleNames)) != 0){
                        warning("**There are folders with \"Dada\" in, so you might have run the function on the folder before\n.
                                Folders without fastq files that fit the F and R patterns will cause an error below")
                }
                
                F_fastq <- character(length = length(folders))
                R_fastq <- character(length = length(folders))
                
                for (i in 1:length(SampleNames)) {
                        CurrentPath <- file.path(path, SampleNames[i])
                        
                        if(sum(grepl(F_pattern, list.files(CurrentPath))) == 0) {
                                stop(paste("F_pattern fits no file in ", CurrentPath))
                        }
                        if(sum(grepl(F_pattern, list.files(CurrentPath))) > 1) {
                                stop(paste("F_pattern fits several files in ", CurrentPath))
                        }
                        if(sum(grepl(R_pattern, list.files(CurrentPath))) == 0) {
                                stop(paste("R_pattern fits no file in ", CurrentPath))
                        }
                        if(sum(grepl(R_pattern, list.files(CurrentPath))) > 1) {
                                stop(paste("R_pattern fits several files in ", CurrentPath))
                        }
                        F_fastq[i] <- file.path(CurrentPath, list.files(CurrentPath)[grepl(F_pattern, list.files(CurrentPath))])
                        R_fastq[i] <- file.path(CurrentPath, list.files(CurrentPath)[grepl(R_pattern, list.files(CurrentPath))])
                        
                }
                
                # silly check
                if (!identical(length(F_fastq), length(R_fastq))) {
                        stop("Number of F and R fastq files differs? But why???")
                }
                
        }
        
        ##############################
        ### Start the log file
        ##############################
        
        DataFolder <- file.path(path2, "Dada_Data")
        dir.create(DataFolder, showWarnings = FALSE)
        
        LogFile <- file.path(DataFolder, "DadaWrapper.log")
        cat("Time after package installation: ", file = LogFile)
        cat(paste(Sys.time(), "\n"), file = LogFile, append = TRUE)
        cat("Package Versions: ", file = LogFile, append = TRUE)
        cat(paste(PackageVersions$Package, ": ", PackageVersions$Version, "; ", sep= ""), file = LogFile, append = TRUE)
        # cat("\n Inputs: ", file = LogFile, append = TRUE)
        # cat(paste("path: ", path, "\n"), file = LogFile, append = TRUE)
        # cat(paste("path2: ", path2, "\n"), file = LogFile, append = TRUE)
        # cat(paste("F_pattern: ", F_pattern, "\n"), file = LogFile, append = TRUE)
        # cat(paste("R_pattern: ", R_pattern, "\n"), file = LogFile, append = TRUE)
        # cat(paste("trimLeft: ", trimLeft[1], trimLeft[2], "\n"), file = LogFile, append = TRUE)
        # cat(paste("tuncLen: ", tuncLen[1], tuncLen[2], "\n"), file = LogFile, append = TRUE)
        # cat(paste("maxN: ", maxN, "\n"), file = LogFile, append = TRUE)
        # cat(paste("maxEE: ", maxEE, "\n"), file = LogFile, append = TRUE)
        # cat(paste("trunQ: ", trunQ, "\n"), file = LogFile, append = TRUE)
        # cat(paste("NSAM.LEARN: ", NSAM.LEARN, "\n"), file = LogFile, append = TRUE)
        # cat(paste("err_F: ", R_pattern, "\n"), file = LogFile, append = TRUE)
        # 
        
        ##############################
        ### Determine the quality scores and save the stats in Data folder
        ##############################
        
        ## Collect quality Score data of the fastQ files and store a data.frame for each sample in the lists FW_QualityStats and RV_QualityStats ####################
        
        
        F_QualityStats <- list()
        R_QualityStats <- list()
        
        for (i in seq_along(F_fastq)) {
                
                Current_FWfq <- F_fastq[i]
                Current_dfFW <- qa(Current_FWfq, n = 1e06)[["perCycle"]]$quality
                # df is a data frame containing for each cycle (nt) the distribution of Quality scores
                Current_dfFW <- dplyr::group_by(Current_dfFW, Cycle)
                
                Current_dfQStatsFW <- dplyr::summarise(
                        Current_dfFW,
                        NoReads = sum(Count),
                        Mean_QS = sum(Count*Score)/sum(Count),
                        SD_QS = sqrt(sum(Count*((Score-Mean_QS)^2))/(NoReads-1)),
                        Median_QS = Score[which(cumsum(Count)/sum(Count) >= .5)][[1]],
                        q25_QS = Score[which(cumsum(Count)/sum(Count) >= .25)][[1]],
                        q75_QS = Score[which(cumsum(Count)/sum(Count) >= .75)][[1]])
                
                ## Check that all reads are of same length
                x <- range(Current_dfQStatsFW$NoReads)
                if (!all.equal(x[1], x[2], tolerance = .Machine$double.eps ^ 0.5)) {
                        stop("Not all reads of same length in file", F_fastq[i])
                }
                
                F_QualityStats[[i]] <- as.data.frame(Current_dfQStatsFW)
                # as.data.frame to un-dplyr the data.frame
                
                ####### collect the same stats for the RV FastQ files
                Current_RVfq <- R_fastq[i]
                Current_dfRV <- qa(Current_RVfq, n = 1e06)[["perCycle"]]$quality
                
                Current_dfRV <- dplyr::group_by(Current_dfRV, Cycle)
                
                Current_dfQStatsRV <- dplyr::summarise(
                        Current_dfRV,
                        NoReads = sum(Count),
                        Mean_QS = sum(Count*Score)/sum(Count),
                        SD_QS = sqrt(sum(Count*((Score-Mean_QS)^2))/(NoReads-1)),
                        Median_QS = Score[which(cumsum(Count)/sum(Count) >= .5)][[1]],
                        q25_QS = Score[which(cumsum(Count)/sum(Count) >= .25)][[1]],
                        q75_QS = Score[which(cumsum(Count)/sum(Count) >= .75)][[1]])
                
                ## Check that all reads are of same length
                x <- range(Current_dfQStatsRV$NoReads)
                if (!all.equal(x[1], x[2], tolerance = .Machine$double.eps ^ 0.5)) {
                        stop("Not all reads of same length in file", R_fastq[i])
                }
                
                R_QualityStats[[i]] <- as.data.frame(Current_dfQStatsRV)
                
                rm(Current_FWfq, Current_dfFW, Current_dfQStatsFW,
                   Current_RVfq, Current_dfRV, Current_dfQStatsRV, x)
                
        }
        
        # add the sample names as names to the lists:
        names(F_QualityStats) <- SampleNames
        names(R_QualityStats) <- SampleNames
        
        save(PackageVersions, F_QualityStats, R_QualityStats, file = file.path(DataFolder, "QualityStats.RData"))
        
        message("*********************** Quality Stats Collected ***********************
                ********************************************************************")
        
        cat("\n*** Quality Stats Collected ***", file = LogFile, append = TRUE)
        
        ##############################
        ### Generate and save some quality plots
        ##############################
        
        PlotFolder <- file.path(path2, "Dada_Plots")
        dir.create(PlotFolder, showWarnings = FALSE)
        
        
        pdf(file = file.path(PlotFolder, "MedianQScore_F_AllNucleotides.pdf"), width = 7, height = 6)
        print(QS_Median_OverviewPlot(F_QualityStats, SampleNames))
        dev.off()
        
        pdf(file = file.path(PlotFolder, "MedianQScore_F_150to240.pdf"), width = 7, height = 6)
        print(QS_Median_OverviewPlot(F_QualityStats, folders, xlim_low = 150, xlim_high = 240))
        dev.off()
        
        pdf(file = file.path(PlotFolder, "MedianQScore_R_AllNucleotides.pdf"), width = 7, height = 6)
        print(QS_Median_OverviewPlot(R_QualityStats, SampleNames))
        dev.off()
        
        pdf(file = file.path(PlotFolder, "MedianQScore_R_150to240.pdf"), width = 7, height = 6)
        print(QS_Median_OverviewPlot(R_QualityStats, folders, xlim_low = 150, xlim_high = 240))
        dev.off()
        
        message("*********************** Plots generated start filtering ***********************
                ********************************************************************")
        
        cat("\n*** Quality Plots generated ***", file = LogFile, append = TRUE)
        
        ##############################
        ### Filtering
        ##############################
        
        FilteredFolder <- file.path(path2, "Dada_FilteredFastqs")
        dir.create(FilteredFolder, showWarnings = FALSE)
        
        filtFs <- file.path(FilteredFolder, paste0(SampleNames, "_F_Filtered.fastq.gz"))
        filtRs <- file.path(FilteredFolder, paste0(SampleNames, "_R_Filtered.fastq.gz"))
        names(filtFs) <- SampleNames
        names(filtRs) <- SampleNames
        
        
        for(i in seq_along(F_fastq)) {
                
                message(paste("**Filtering sample No:", i, "called:", SampleNames[i]), " **")
                fastqPairedFilter(c(F_fastq[i], R_fastq[i]), c(filtFs[i], filtRs[i]), trimLeft = trimLeft, 
                                  truncLen = truncLen, maxN = maxN, maxEE = maxEE, truncQ = truncQ, 
                                  compress=TRUE, verbose=TRUE)
        }
        
        TimePassed <- proc.time()-ptm
        
        # check if files have been created
        if(!all(file.exists(filtFs)) & !all(file.exists(filtRs))) {
                cat("\n*** ERROR: Not all filtered files were created, maybe trimming impossible***", file = LogFile, append = TRUE)
                stop("** Not all filtered files were created, maybe trimming impossible")
                
        }
        
        message("*********************** Filtering Done ***********************
                ********************************************************************")
        
        cat("\n*** Filtering Done ***", file = LogFile, append = TRUE)
        cat(paste("\nTime Passed: ", TimePassed[3]), file = LogFile, append = TRUE)
        
        
        ##############################
        ### Estimate err_F matrix
        ##############################
        
        if(is.null(err_F)){
                
                if(is.null(NSAM.LEARN)){
                        NSAM.LEARN <- length(SampleNames)
                        message("Your NSAM.LEARN was \"NULL\" all samples are used for err_F estimation")
                } else if(NSAM.LEARN >= length(SampleNames)){
                        NSAM.LEARN <- length(SampleNames)
                        message("Your NSAM.LEARN covered all samples, i.e. all samples are used for err_F estimation")
                }
                
                
                drp.learnF <- derepFastq(filtFs[1:NSAM.LEARN])
                
                ReadsForErrFEstimation <- sum(sapply(1:NSAM.LEARN, function(x) sum(drp.learnF[[x]]$uniques)))
                
                message(paste("**", ReadsForErrFEstimation, "reads will be used for err_F estimation **"))
                
                dd.learnF <- dada(drp.learnF, err=NULL, selfConsist=TRUE, multithread=TRUE)
                # ## NB: The determined error rates are all identical, so here samples were pooled that is why it took so long
                # identical(dd.learnF[[1]]$err_out, dd.learnF[[2]]$err_out, dd.learnF[[3]]$err_out)
                
                err_F <- dd.learnF[[1]]$err_out
                
                pdf(file = file.path(PlotFolder, "errorRates_F_Sample1.pdf"), width = 10, height = 10)
                print(plotErrors(dd.learnF[[1]], nominalQ=TRUE))
                dev.off()
                
                message("*********************** err_F has been estimated ***********************
                ********************************************************************")
                cat("\n*** err_F has been estimated ***", file = LogFile, append = TRUE)
                cat(paste("\nSamples used for err_F estimation: ", NSAM.LEARN), file = LogFile, append = TRUE)
                cat(paste("\nReads used for err_F estimation: ", ReadsForErrFEstimation), file = LogFile, append = TRUE)
                TimePassed <- proc.time()-ptm
                cat(paste("\nTime Passed: ", TimePassed[3]), file = LogFile, append = TRUE)
                
        }
        
        
        ##############################
        ### Estimate err_R matrix
        ##############################
        
        if(is.null(err_R)){
                
                if(is.null(NSAM.LEARN)){
                        NSAM.LEARN <- length(SampleNames)
                        message("Your NSAM.LEARN was \"NULL\" all samples are used for err_R estimation")
                } else if(NSAM.LEARN >= length(SampleNames)){
                        NSAM.LEARN <- length(SampleNames)
                        message("Your NSAM.LEARN covered all samples, i.e. all samples are used for err_R estimation")
                }
                
                
                drp.learnR <- derepFastq(filtRs[1:NSAM.LEARN])
                
                ReadsForErrREstimation <- sum(sapply(1:NSAM.LEARN, function(x) sum(drp.learnR[[x]]$uniques)))
                
                message(paste("**", ReadsForErrREstimation, "reads will be used for err_R estimation **"))
                
                dd.learnR <- dada(drp.learnR, err=NULL, selfConsist=TRUE, multithread=TRUE)
                
                err_R <- dd.learnR[[1]]$err_out
                
                pdf(file = file.path(PlotFolder, "errorRates_R_Sample1.pdf"), width = 10, height = 10)
                print(plotErrors(dd.learnR[[1]], nominalQ=TRUE))
                dev.off()
                
                message("*********************** err_R has been estimated ***********************
                        ********************************************************************")
                cat("\n*** err_R has been estimated ***", file = LogFile, append = TRUE)
                cat(paste("\nSamples used for err_R estimation: ", NSAM.LEARN), file = LogFile, append = TRUE)
                cat(paste("\nReads used for err_R estimation: ", ReadsForErrREstimation), file = LogFile, append = TRUE)
                TimePassed <- proc.time()-ptm
                cat(paste("\nTime Passed: ", TimePassed[3]), file = LogFile, append = TRUE)
                
        }
        
        
        ##############################
        ### Denoising (dada command for all samples)
        ##############################
        
        #### NB: write here an alternative if err_F and err_R were NULL and 
        ## NSAM.LEARN = length(SampleNames)
        ## AND 
        
        rm(drp.learnF, drp.learnR, dd.learnF, dd.learnR)
        
        mergers <- vector("list", length(SampleNames))
        names(mergers) <- SampleNames
        NoFilteredReads <- vector("numeric", length(SampleNames))
        names(NoFilteredReads) <- SampleNames
        bimFs <- vector("list", length(SampleNames))
        names(bimFs) <- SampleNames
        bimRs <- vector("list", length(SampleNames))
        names(bimRs) <- SampleNames
        
        for(sam in SampleNames) {
                cat("Processing:", sam, "\n")
                cat(paste("\nDenoising sample:", sam), file = LogFile, append = TRUE)
                derepF <- derepFastq(filtFs[[sam]])
                NoFilteredReads[sam] <- sum(derepF$uniques)
                ddF <- dada(derepF, err=err_F, multithread=TRUE) 
                bimFs[[sam]] <- isBimeraDenovo(ddF, verbose=TRUE)
                derepR <- derepFastq(filtRs[[sam]])
                ddR <- dada(derepR, err=err_R, multithread=TRUE)
                bimRs[[sam]] <- isBimeraDenovo(ddR, verbose=TRUE)
                merger <- mergePairs(ddF, derepF, ddR, derepR, minOverlap = minOverlap, maxMismatch = maxMismatch)
                mergers[[sam]] <- merger
        }
        
        if(!all((sapply(1:length(mergers), function(i) dim(mergers[[i]])[1]))!=0)){
                message("**Not in all samples merged amplicons could be found! maybe minOverlap too strict??**")
                cat("\n*** Not in all samples merged amplicons could be found! maybe minOverlap too strict?? ***", file = LogFile, append = TRUE)
        }
        
        if(all((sapply(1:length(mergers), function(i) dim(mergers[[i]])[1]))==0)){
                cat("\n*** ERROR: In none of the samples merged amplicons could be found! maybe minOverlap too strict?? ***", file = LogFile, append = TRUE)
                stop("**In none of the samples merged amplicons could be found! maybe minOverlap too strict??**")
        }
        
        rm(derepF, derepR, ddF, ddR, merger)
        
        message("*********************** all samples denoised ***********************
                        ********************************************************************")
        cat("\n*** all samples denoised ***", file = LogFile, append = TRUE)
        TimePassed <- proc.time()-ptm
        cat(paste("\nTime Passed: ", TimePassed[3]), file = LogFile, append = TRUE)
        
        
        ##############################
        ### Bimera removal
        ##############################
        
        mergers.nochim <- mergers
        for (i in seq_along(mergers)) {
                mergers.nochim[[i]] <- mergers[[i]][!bimFs[[i]][mergers[[i]]$forward] & !bimRs[[i]][mergers[[i]]$reverse],]
        }
        
        message("*********************** Bimeras removed ***********************
                        ********************************************************************")
        cat("\n*** bimeras removed ***", file = LogFile, append = TRUE)
        
        
        ##############################
        ### Generating sequence table
        ##############################
        
        seqtab <- makeSequenceTable(mergers.nochim)
        
        save(seqtab, mergers, mergers.nochim, bimFs, bimRs, NoFilteredReads, file = file.path(DataFolder, "DenoiseData.RData"))
        
        message("*********************** Sequence table generated, Data saved ***********************
                        ********************************************************************")
        cat("\n*** Sequence table generated, Data saved ***", file = LogFile, append = TRUE)
        
        
        ##############################
        ### Summary Plots
        ##############################
        
        pdf(file = file.path(PlotFolder, "NoReads_AllSamples.pdf"), width = 17, height = 6)
        print(NoReads_Steps(QStatsList = F_QualityStats, NoFilteredReads = NoFilteredReads, mergers = mergers, mergers.nochim = mergers.nochim, SampleNames = SampleNames, sort = TRUE))
        dev.off()
        
        FinalNumbers <- data.frame(Sample = rownames(seqtab), UniqueAmplicons = rowSums(seqtab != 0), NoAmplicons = rowSums(seqtab))
        
        Tr <- TotalandUniqueAmplicons(FinalNumbers)
        pdf(file = file.path(PlotFolder, "TotalvsUniqueAmplicons.pdf"), width = 17, height = 6)
        print(Tr)
        dev.off()
        
        
        Tr2 <- ggplot(FinalNumbers, aes(x = NoAmplicons, y = UniqueAmplicons)) +
                geom_point() +
                geom_smooth(method = "lm", se = FALSE) +
                xlab("Total Amplicons") +
                ylab("Unique Amplicons") +
                theme_bw() +
                theme(panel.grid.minor = element_blank(),
                      panel.grid.major.y = element_line(color = "#999999", size = .15),
                      panel.grid.major.x = element_line(color = "#999999", size = .15))
        
        
        pdf(file = file.path(PlotFolder, "AssociationTotaltoUniqueAmplicons.pdf"), width = 7, height = 6)
        print(Tr)
        dev.off()
        
        
        FinalNumbers2 <- data.frame(InNoSamples = colSums(seqtab != 0))
        
        Trr <- ggplot(data = FinalNumbers2, aes(x = InNoSamples))  + 
                geom_histogram(binwidth = 3, col = "black", fill = "#E69F00")+
                geom_rug() +
                xlab("Present in No Samples") + 
                ylab("Count") +
                theme_bw() + 
                ggtitle(paste("Total No of unique amplicons:", dim(seqtab)[2], "Only in 1 sample:", sum(FinalNumbers2 ==1))) +
                theme(panel.grid.minor = element_blank(),
                      panel.grid.major.y = element_blank(),
                      panel.grid.major.x = element_line(color = "#999999", size = .15)) +
                coord_cartesian(ylim = c(0,150))
        
        pdf(file = file.path(PlotFolder, "HistogramAmpliconsinNoSamples.pdf"), width = 7, height = 6)
        print(Trr)
        dev.off()
        
}
