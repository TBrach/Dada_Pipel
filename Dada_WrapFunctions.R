#######################################
### FUNCTION: Dada2_wrap
#######################################

# Function that runs the dada2 pipeline up to the sequence table (seqtab)
## Input
# path: The path to the folder containing the sample folders with the fastq files: 
# F_pattern: a regular expression to find the fastq files with the forward reads in the sample folders
# R_pattern: a regular expression to find the fastq files with the reverse reads in the sample folders
# path2: default = NULL, the function creates three folders Dada_Plots, Dada_Data and Dada_FilteredFastqs, it creates them by default in the path folder, when path2 is given 
# the folders will instead be generated in path2
# trimLeft: = trimLeft from the fastqPairedFilter command (how many nucleotids will be clipped from the beginning of the reads)
# truncLen: = truncLen from the fastqPairedFilter command (where the reads will be truncated)
# maxN: = maxN from the fastqPairedFilter command
# maxEE: = maxEE from the fastqPairedFilter command, the maximum expected error allowed to let a read pass the filter
# trunQ: = truncQ from the fastqPairedFilter command
# NSAM.LEARN: default = NULL: the number of samples used to estimate the F and R error matrixes. When NULL all samples are used. 
# Should account for a number of samples that in total have about 1 milliion filtered reads (number will be given in the log file)
# err_F: the error matrix for the dada command for the forward reads. Default NULL then the error matrix will be estimated (using SelfConsit = TRUE)
# err_R: analogous to err_F for the reverse reads
# minOverlap: = minOverlap from the mergePairs command
# maxMismatch: = maxMismatch from the mergePairs command
# F_QualityStats and R_QualityStats: NULL by default, if given e.g. from Dada_QualityCheck the quality stats collection part is jumped over
# filtFs and filtRs: NULL by default, if given and FilteredFolder with files exist, Filtering can be jumped over.
## Output
# PLOTS:
# Several Plots saved as pdfs in generated Dada_Plots folder # NB: the plots are currently assuming 250 nt
# DATA AND LOG FILE:
# Data and the file DadaWrapper.log are saved in the generated folder: Dada_Data
# DenoisedData.RData contains: seqtab (final sequence table), mergers, mergers.nochim, bimFs, bimRs, 
# ReadSummary (data frame summarising the reads and amplicons at the different stages), SamplesFor_errF, SamplesFor_errR
# QualityStats.RData contains PackageVersions, F_QualityStats, R_QualityStats
# FILTERED FASTQ files
# are stored in the generated folder Dada_FilteredFastqs
# NB the filtered fastq files are currently named "SampleName"_F_Filtered.fastq.gz and "SampleName"_R_Filtered.fastq.gz

Dada2_wrap <- function(path, F_pattern, R_pattern, path2 = NULL,
                       trimLeft = c(10,10), truncLen = c(220, 160), 
                       maxN = 0, maxEE = 2, truncQ = 2,
                       NSAM.LEARN = NULL,
                       err_F = NULL,
                       err_R = NULL,
                       minOverlap = 20,
                       maxMismatch = 0,
                       F_QualityStats = NULL,
                       R_QualityStats = NULL,
                       filtFs = NULL,
                       filtRs = NULL) {
        
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
        ### Construct character vectors to the FW and RV fastq files
        ##############################
        
        if(is.null(path2)){
                
                path2 = path
        }
        
        
        folders <- list.dirs(path, recursive = FALSE, full.names = FALSE)
        
        if(length(folders) != 0) {
                
                SampleNames <- folders
                
                if(sum(grepl("^Dada", SampleNames)) != 0){
                        warning("**There are folders starting with \"Dada\" in your path folder, maybe you have run the function on this path folder before\n.
                                The folders starting with Dada will not be considered sample folders!!\n Files within Dada_Data Dada_FilteredFastqs and Dada_Plots will be overwritten**")
                }
                
                # exclude Folders starting with "Dada" from the folders considered as sample fodlers
                if(length(grep("^Dada", SampleNames))!=0) {
                        SampleNames <- SampleNames[-grep("^Dada", SampleNames)]
                }
                
                F_fastq <- character(length = length(SampleNames))
                R_fastq <- character(length = length(SampleNames))
                
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
                
        } else {
                
                stop("No sample folders were found in the given path! Currently the Dada2_wrap function can only handle the situation where the fastq files are in separate folders for each sample.
                     These folders have to be in the \"path\" folder. Other situations have to be added.")
        }
        
        ##############################
        ### Start the log file
        ##############################
        
        DataFolder <- file.path(path2, "Dada_Data")
        
        ## ATTENTION: FUNCTION DELETES ALL FILES ALREADY PRESRENT in DataFolder
        if(file.exists(DataFolder)){
                file.remove(list.files(DataFolder, full.names = TRUE))
        }
        
        
        dir.create(DataFolder, showWarnings = FALSE)
        
        ptm <- proc.time()
        
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
        
        if(is.null(F_QualityStats) | is.null(R_QualityStats)) {
                
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
                
        } else {
                
                if(names(F_QualityStats) != SampleNames | names(R_QualityStats) != SampleNames) {
                        
                        stop("The given F_QualityStats or R_QualityStats do not fit to the SampleNames in path!")
                }
                
        }
        
        save(PackageVersions, F_QualityStats, R_QualityStats, file = file.path(DataFolder, "QualityStats.RData"))
        
        message("*********************** Quality Stats Collected ***********************
                ********************************************************************")
        
        cat("\n*** Quality Stats Collected ***", file = LogFile, append = TRUE)
        TimePassed <- proc.time()-ptm
        cat("\nTime after Quality Stats collection: ", file = LogFile)
        cat(paste(Sys.time(), "\n"), file = LogFile, append = TRUE)
        cat(paste("\nTime Passed in total: ", TimePassed[3]), file = LogFile, append = TRUE)
        cat("\n*** Start generating Quality Plots ***", file = LogFile, append = TRUE)
        
        ##############################
        ### Generate and save some quality plots
        ##############################
        
        PlotFolder <- file.path(path2, "Dada_Plots")
        
        ## ATTENTION: FUNCTION DELETES ALL FILES ALREADY PRESRENT in PlotFolder
        if(file.exists(PlotFolder)){
                file.remove(list.files(PlotFolder, full.names = TRUE))
        }
        
        dir.create(PlotFolder, showWarnings = FALSE)
        
        
        pdf(file = file.path(PlotFolder, "MedianQScore_F_AllNucleotides.pdf"), width = 7, height = 6)
        print(QS_Median_OverviewPlot(F_QualityStats, SampleNames))
        dev.off()
        
        pdf(file = file.path(PlotFolder, "MedianQScore_F_150to240.pdf"), width = 7, height = 6)
        print(QS_Median_OverviewPlot(F_QualityStats, SampleNames, xlim_low = 150, xlim_high = 240))
        dev.off()
        
        pdf(file = file.path(PlotFolder, "MedianQScore_R_AllNucleotides.pdf"), width = 7, height = 6)
        print(QS_Median_OverviewPlot(R_QualityStats, SampleNames))
        dev.off()
        
        pdf(file = file.path(PlotFolder, "MedianQScore_R_150to240.pdf"), width = 7, height = 6)
        print(QS_Median_OverviewPlot(R_QualityStats, SampleNames, xlim_low = 150, xlim_high = 240))
        dev.off()
        
        message("*********************** Plots generated start filtering ***********************
                ********************************************************************")
        
        cat("\n*** Quality Plots generated ***", file = LogFile, append = TRUE)
        cat("\n*** Start Filtering ***", file = LogFile, append = TRUE)
        
        ##############################
        ### Filtering
        ##############################
        
        FilteredFolder <- file.path(path2, "Dada_FilteredFastqs")
        
        
        if(is.null(filtFs) | is.null(filtRs)){
                
                ## Not sure if this is wanted, deleting all files that are already in the DataFolder folder
                if(file.exists(FilteredFolder)){
                        file.remove(list.files(FilteredFolder, full.names = TRUE))
                }
                
                dir.create(FilteredFolder, showWarnings = FALSE)
                
                filtFs <- file.path(FilteredFolder, paste0(SampleNames, "_F_Filtered.fastq.gz"))
                filtRs <- file.path(FilteredFolder, paste0(SampleNames, "_R_Filtered.fastq.gz"))
                names(filtFs) <- SampleNames
                names(filtRs) <- SampleNames
                
                
                for(i in seq_along(F_fastq)) {
                        
                        message(paste("**Filtering sample No:", i, "called:", SampleNames[i]), " **")
                        cat(paste("\n**Filtering sample No:", i, "called:", SampleNames[i]), file = LogFile, append = TRUE)
                        fastqPairedFilter(c(F_fastq[i], R_fastq[i]), c(filtFs[i], filtRs[i]), trimLeft = trimLeft, 
                                          truncLen = truncLen, maxN = maxN, maxEE = maxEE, truncQ = truncQ, 
                                          compress=TRUE, verbose=TRUE)
                }
                
                # check if files have been created
                if(!all(file.exists(filtFs)) & !all(file.exists(filtRs))) {
                        cat("\n*** ERROR: Not all filtered files were created, maybe trimming impossible***", file = LogFile, append = TRUE)
                        stop("** Not all filtered files were created, maybe trimming impossible")
                        
                }
                
                message("*********************** Filtering Done ***********************
                ********************************************************************")
                cat("\n*** Filtering Done ***", file = LogFile, append = TRUE)
                
        } else {
                
                if(names(filtFs) != SampleNames | names(filtRs) != SampleNames) {
                        
                        cat("\n*** ERROR: The given filtFs or filtRs do not have sample names!***", file = LogFile, append = TRUE)
                        stop("The given filtFs or filtRs do not have sample names!")
                }
                
                if(!all(file.exists(filtFs)) & !all(file.exists(filtRs))) {
                        cat("\n*** ERROR: Not all files in the given filtFs or filtRs existed***", file = LogFile, append = TRUE)
                        stop("** Not all files in the given filtFs or filtRs existed")
                        
                }
                
                message("*********************** Filtered files were given ***********************
                ********************************************************************")
                cat("\n*** Filtered files were given ***", file = LogFile, append = TRUE)
                
                
        }
        
        
        save(PackageVersions, F_QualityStats, R_QualityStats, filtFs, filtRs, file = file.path(DataFolder, "QualityStats.RData"))
        
        cat("\nTime after filtering step: ", file = LogFile)
        cat(paste(Sys.time(), "\n"), file = LogFile, append = TRUE)
        TimePassed <- proc.time()-ptm
        cat(paste("\nTime Passed in total: ", TimePassed[3]), file = LogFile, append = TRUE)
        cat("\n*** Start estimating err_F if not given ***", file = LogFile, append = TRUE)
        
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
                
                # Random selection of samples used for err_F estimation unless all samples are used
                
                if(NSAM.LEARN == length(SampleNames)){
                        SamplesFor_errF <- filtFs
                } else {
                        prng <- .Random.seed
                        SamplesFor_errF <- sample(filtFs, NSAM.LEARN)
                        attr(SamplesFor_errF, "seed") <- prng
                }
                
                drp.learnF <- derepFastq(SamplesFor_errF)
                
                if(class(drp.learnF) == "list") {
                        ReadsForErrFEstimation <- sum(sapply(1:length(drp.learnF), function(x) sum(drp.learnF[[x]]$uniques)))
                } else {
                        ReadsForErrFEstimation <- sum(drp.learnF$uniques)
                }
                
                message(paste("**", ReadsForErrFEstimation, "reads will be used for err_F estimation **"))
                
                dd.learnF <- dada(drp.learnF, err=NULL, selfConsist=TRUE, multithread=TRUE)
                # ## NB: The determined error rates are all identical, so here samples were pooled that is why it took so long
                # identical(dd.learnF[[1]]$err_out, dd.learnF[[2]]$err_out, dd.learnF[[3]]$err_out)
                
                if(class(dd.learnF) == "list") {
                        err_F <- dd.learnF[[1]]$err_out
                } else {
                        err_F <- dd.learnF$err_out
                }
                
                pdf(file = file.path(PlotFolder, "errorRates_F_Sample1.pdf"), width = 10, height = 10)
                if(class(dd.learnF) == "list") {
                        print(plotErrors(dd.learnF[[1]], nominalQ=TRUE))
                } else {
                        print(plotErrors(dd.learnF, nominalQ=TRUE))
                }
                dev.off()
                
                save(err_F, PackageVersions, F_QualityStats, R_QualityStats, filtFs, filtRs, file = file.path(DataFolder, "QualityStats.RData"))
                
                message("*********************** err_F has been estimated ***********************
                ********************************************************************")
                cat("\n*** err_F has been estimated ***", file = LogFile, append = TRUE)
                cat(paste("\nNumber of samples used for err_F estimation: ", NSAM.LEARN), file = LogFile, append = TRUE)
                cat("\nSamples used for err_F estimation: \n", file = LogFile, append = TRUE)
                cat(paste(SampleNames[sort(match(names(SamplesFor_errF), SampleNames))]), file = LogFile, append = TRUE)
                cat(paste("\nReads used for err_F estimation: ", ReadsForErrFEstimation), file = LogFile, append = TRUE)
                TimePassed <- proc.time()-ptm
                cat("\nTime after err_F estimation: ", file = LogFile)
                cat(paste(Sys.time(), "\n"), file = LogFile, append = TRUE)
                cat(paste("\nTime Passed in total: ", TimePassed[3]), file = LogFile, append = TRUE)
                cat("\n*** Start estimating err_R if not given ***", file = LogFile, append = TRUE)
                
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
                
                # Random selection of samples used for err_F estimation unless all samples are used
                if(NSAM.LEARN == length(SampleNames)){
                        SamplesFor_errR <- filtRs
                } else {
                        prng <- .Random.seed
                        SamplesFor_errR <- sample(filtRs, NSAM.LEARN)
                        attr(SamplesFor_errR, "seed") <- prng
                }
                
                drp.learnR <- derepFastq(SamplesFor_errR)
                
                if(class(drp.learnR) == "list") {
                        ReadsForErrREstimation <- sum(sapply(1:length(drp.learnR), function(x) sum(drp.learnR[[x]]$uniques)))
                } else {
                        ReadsForErrREstimation <- sum(drp.learnR$uniques)
                }
                
                message(paste("**", ReadsForErrREstimation, "reads will be used for err_R estimation **"))
                
                dd.learnR <- dada(drp.learnR, err=NULL, selfConsist=TRUE, multithread=TRUE)
                
                if(class(dd.learnR) == "list") {
                        err_R <- dd.learnR[[1]]$err_out
                } else {
                        err_R <- dd.learnR$err_out
                }
                
                pdf(file = file.path(PlotFolder, "errorRates_R_Sample1.pdf"), width = 10, height = 10)
                if(class(dd.learnR) == "list") {
                        print(plotErrors(dd.learnR[[1]], nominalQ=TRUE))
                } else {
                        print(plotErrors(dd.learnR, nominalQ=TRUE))
                }
                dev.off()
                
                save(err_F, err_R, PackageVersions, F_QualityStats, R_QualityStats, filtFs, filtRs, file = file.path(DataFolder, "QualityStats.RData"))
                
                message("*********************** err_R has been estimated ***********************
                        ********************************************************************")
                cat("\n*** err_R has been estimated ***", file = LogFile, append = TRUE)
                cat(paste("\nNumber of samples used for err_R estimation: ", NSAM.LEARN), file = LogFile, append = TRUE)
                cat("\nSamples used for err_R estimation: \n", file = LogFile, append = TRUE)
                cat(paste(SampleNames[sort(match(names(SamplesFor_errR), SampleNames))]), file = LogFile, append = TRUE)
                cat(paste("\nReads used for err_R estimation: ", ReadsForErrREstimation), file = LogFile, append = TRUE)
                TimePassed <- proc.time()-ptm
                cat("\nTime after err_F estimation: ", file = LogFile)
                cat(paste(Sys.time(), "\n"), file = LogFile, append = TRUE)
                cat(paste("\nTime Passed in total: ", TimePassed[3]), file = LogFile, append = TRUE)
                cat("\n*** Start denoising data, bimera detection, and merging of reads into amplicons***", file = LogFile, append = TRUE)
                
        }
        
        
        ##############################
        ### Denoising (dada command for all samples) and bimeara identification
        ##############################
        
        #### NB: if all samples have been used to generate drp.learnR, dd.learnR, drp.learnF, dd.learnF, these should of course be used
        # instead of running again through all samples, therefore the following if check
        
        if (NSAM.LEARN == length(SampleNames) & exists("drp.learnR", inherits = FALSE) & exists("drp.learnF", inherits = FALSE) &
            exists("dd.learnR", inherits = FALSE) & exists("dd.learnF", inherits = FALSE)) {
                # if drp.learnF and R have been generated in this function and all samples where used for doing so, then:
                
                #NB: the follwoing demands that the samples were denoised in the given order, i.e. just 1:NSAM.LEARN was used not the sample command!!
                names(dd.learnF) <- SampleNames
                names(drp.learnF) <- SampleNames
                names(dd.learnR) <- SampleNames
                names(drp.learnR) <- SampleNames
                
                ## only for the ReadSummary later
                Uniques_F <- sapply(1:length(SampleNames), function(i) length(drp.learnF[[i]]$uniques))
                Uniques_R <- sapply(1:length(SampleNames), function(i) length(drp.learnR[[i]]$uniques))
                Denoised_F <- sapply(1:length(SampleNames), function(i) length(dd.learnF[[i]]$denoised))
                Denoised_R <- sapply(1:length(SampleNames), function(i) length(dd.learnR[[i]]$denoised))
                names(Uniques_F) <- SampleNames
                names(Uniques_R) <- SampleNames
                names(Denoised_F) <- SampleNames
                names(Denoised_R) <- SampleNames
                ##
                
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
                        #cat(paste("\nDenoising sample:", sam), file = LogFile, append = TRUE)
                        derepF <- drp.learnF[[sam]]
                        NoFilteredReads[sam] <- sum(derepF$uniques)
                        ddF <- dd.learnF[[sam]] 
                        bimFs[[sam]] <- isBimeraDenovo(ddF, verbose=TRUE)
                        derepR <- drp.learnR[[sam]]
                        ddR <- dd.learnR[[sam]]
                        bimRs[[sam]] <- isBimeraDenovo(ddR, verbose=TRUE)
                        merger <- mergePairs(ddF, derepF, ddR, derepR, minOverlap = minOverlap, maxMismatch = maxMismatch)
                        mergers[[sam]] <- merger
                        
                        rm(derepF, derepR, ddF, ddR, merger)
                }
                
                rm(drp.learnF, drp.learnR, dd.learnF, dd.learnR)
                
        } else {
                
                rm(drp.learnF, drp.learnR, dd.learnF, dd.learnR)
                
                mergers <- vector("list", length(SampleNames))
                names(mergers) <- SampleNames
                NoFilteredReads <- vector("numeric", length(SampleNames))
                names(NoFilteredReads) <- SampleNames
                bimFs <- vector("list", length(SampleNames))
                names(bimFs) <- SampleNames
                bimRs <- vector("list", length(SampleNames))
                names(bimRs) <- SampleNames
                Uniques_F <- vector("numeric", length(SampleNames))
                names(Uniques_F) <- SampleNames
                Uniques_R <- vector("numeric", length(SampleNames))
                names(Uniques_R) <- SampleNames
                Denoised_F <- vector("numeric", length(SampleNames))
                names(Denoised_F) <- SampleNames
                Denoised_R <- vector("numeric", length(SampleNames))
                names(Denoised_R) <- SampleNames
                
                for(sam in SampleNames) {
                        cat("Processing:", sam, "\n")
                        cat(paste("\nDenoising sample:", sam), file = LogFile, append = TRUE)
                        derepF <- derepFastq(filtFs[[sam]])
                        NoFilteredReads[sam] <- sum(derepF$uniques)
                        ddF <- dada(derepF, err=err_F, multithread=TRUE) 
                        bimFs[[sam]] <- isBimeraDenovo(ddF, verbose=TRUE)
                        Uniques_F[sam] <- length(derepF$uniques)
                        Denoised_F[sam] <- length(ddF$denoised)
                        
                        derepR <- derepFastq(filtRs[[sam]])
                        ddR <- dada(derepR, err=err_R, multithread=TRUE)
                        bimRs[[sam]] <- isBimeraDenovo(ddR, verbose=TRUE)
                        Uniques_R[sam] <- length(derepR$uniques)
                        Denoised_R[sam] <- length(ddR$denoised)
                        merger <- mergePairs(ddF, derepF, ddR, derepR, minOverlap = minOverlap, maxMismatch = maxMismatch)
                        mergers[[sam]] <- merger
                        
                        rm(derepF, derepR, ddF, ddR, merger)
                }
                
        }
        
        if(!all((sapply(1:length(mergers), function(i) dim(mergers[[i]])[1]))!=0)){
                message("**Not in all samples merged amplicons could be found! maybe minOverlap too strict??**")
                cat("\n*** Not in all samples merged amplicons could be found! maybe minOverlap too strict?? ***", file = LogFile, append = TRUE)
        }
        
        if(all((sapply(1:length(mergers), function(i) dim(mergers[[i]])[1]))==0)){
                cat("\n*** ERROR: In none of the samples merged amplicons could be found! maybe minOverlap too strict?? ***", file = LogFile, append = TRUE)
                stop("**In none of the samples merged amplicons could be found! maybe minOverlap too strict??**")
        }
        
        message("*********************** all samples denoised, bimeras identified, mergerd amplicons generated ***********************
                        ********************************************************************")
        cat("\n*** all samples denoised, bimeras identified, mergerd amplicons generated ***", file = LogFile, append = TRUE)
        TimePassed <- proc.time()-ptm
        cat("\nTime after denoising: ", file = LogFile)
        cat(paste(Sys.time(), "\n"), file = LogFile, append = TRUE)
        cat(paste("\nTime Passed in total: ", TimePassed[3]), file = LogFile, append = TRUE)
        cat("\n*** Start removing bimera ***", file = LogFile, append = TRUE)
        
        
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
        cat("\n*** Start generating Sequence table ***", file = LogFile, append = TRUE)
        
        
        ##############################
        ### Generating sequence table
        ##############################
        
        seqtab <- makeSequenceTable(mergers.nochim)
        
        ## generate also a read summary data frame
        QStatsList <- F_QualityStats
        for (i in seq_along(QStatsList)) {
                QStatsList[[i]]$Sample <- names(QStatsList[i])
        }
        ReadSummary <- do.call("rbind",QStatsList)
        ReadSummary <- ReadSummary[!duplicated(ReadSummary$Sample), c("Sample", "NoReads")]
        ReadSummary$Filtered <- NoFilteredReads
        ReadSummary$Merged <- sapply(1:length(SampleNames), function(i) sum(mergers[[i]]$abundance))
        ReadSummary$NoChimera <- sapply(1:length(SampleNames), function(i) sum(mergers.nochim[[i]]$abundance))
        ReadSummary$Uniques_F <- Uniques_F
        ReadSummary$Denoised_F <- Denoised_F
        ReadSummary$bimera_F <- sapply(1:length(SampleNames), function(i) sum(bimFs[[i]]))
        ReadSummary$Uniques_R <- Uniques_R
        ReadSummary$Denoised_R <- Denoised_R
        ReadSummary$bimera_R <- sapply(1:length(SampleNames), function(i) sum(bimRs[[i]]))
        ReadSummary$Unique_Amplicons <- sapply(1:length(SampleNames), function(i) dim(mergers[[i]])[1])
        ReadSummary$Unique_Amplicons_nochim <- sapply(1:length(SampleNames), function(i) dim(mergers.nochim[[i]])[1])
        rownames(ReadSummary) <- NULL
        
        if(!exists("SamplesFor_errF")){SamplesFor_errF = NULL}
        if(!exists("SamplesFor_errR")){SamplesFor_errR = NULL}
        
        save(seqtab, mergers, mergers.nochim, bimFs, bimRs, ReadSummary, SamplesFor_errF,
             SamplesFor_errR, file = file.path(DataFolder, "DenoisedData.RData"))
        
        message("*********************** Sequence table generated, Data saved ***********************
                        ********************************************************************")
        cat("\n*** Sequence table generated, Data saved ***", file = LogFile, append = TRUE)
        cat("\n*** Start generating Summary Plots ***", file = LogFile, append = TRUE)
        
        ##############################
        ### Summary Plots
        ##############################
        
        # the width of the plots probably needs adjustment
        width = 5 + 0.5*(length(SampleNames)/10)
        
        # Plot the number of reads at the different steps for each sample (see also ReadSummary)
        pdf(file = file.path(PlotFolder, "NoReads_AllSamples.pdf"), width = width, height = 6)
        print(NoReads_Steps(QStatsList = F_QualityStats, NoFilteredReads = NoFilteredReads, mergers = mergers, mergers.nochim = mergers.nochim, SampleNames = SampleNames, sort = TRUE))
        dev.off()
        
        # PLot the total number of amplicons against the number of unique amplicons for each sample
        FinalNumbers <- data.frame(Sample = rownames(seqtab), UniqueAmplicons = rowSums(seqtab != 0), NoAmplicons = rowSums(seqtab))
        
        Tr <- TotalandUniqueAmplicons(FinalNumbers = FinalNumbers, seqtab = seqtab)
        pdf(file = file.path(PlotFolder, "TotalvsUniqueAmplicons.pdf"), width = width, height = 6)
        print(Tr)
        dev.off()
        
        # Plot a linear regression line to illustrate the possible association between the total number of amplicons and the number of 
        # unique amplicons
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
        print(Tr2)
        dev.off()
        
        # Plot a histogram illustrating in how many samples the amplicons are present
        if(length(SampleNames) < 10) {binwidth = 1}
        if(length(SampleNames) > 10 & length(SampleNames) < 100) {binwidth = 2}
        if(length(SampleNames) > 100) {binwidth = 3}
        
        FinalNumbers2 <- data.frame(InNoSamples = colSums(seqtab != 0))
        
        Trr <- ggplot(data = FinalNumbers2, aes(x = InNoSamples))  + 
                geom_histogram(binwidth = binwidth, col = "black", fill = "#E69F00")+
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
        cat("\n*** Summary Plots generated, Function done ***", file = LogFile, append = TRUE)
}



#######################################
### FUNCTION: Dada2_QualityCheck
#######################################

# runs the first part of the Dada2_wrap function, creating the quality plots and data. Based on these one can decide on the filtering
# parameters when using the Dada2_wrap function subsequently
## Input
# path: The path to the folder containing the sample folders with the fastq files: 
# F_pattern: a regular expression to find the fastq files with the forward reads in the sample folders
# R_pattern: a regular expression to find the fastq files with the reverse reads in the sample folders
# path2: default = NULL, the function creates the folders Dada_Plots and Dada_Data, it creates them by default in the path folder, when path2 is given 
# the folders will instead be generated in path2
## Output
# PLOTS:
# Quality Plots saved as pdfs in the generated Dada_Plots folder # NB: the plots are currently assuming 250 nt
# DATA:
# saved in the Dada_Data folder:
# QualityStats.RData contains PackageVersions, F_QualityStats, and R_QualityStats

Dada2_QualityCheck <- function(path, F_pattern, R_pattern, path2 = NULL) {
        
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
        ### Construct character vectors to the FW and RV fastq files
        ##############################
        if(is.null(path2)){
                
                path2 = path
        }
        
        
        folders <- list.dirs(path, recursive = FALSE, full.names = FALSE)
        
        if(length(folders) != 0) {
                
                SampleNames <- folders
                
                if(sum(grepl("^Dada", SampleNames)) != 0){
                        warning("**There are folders starting with \"Dada\" in your path folder, maybe you have run the function on this path folder before\n.
                                The folders starting with Dada will not be considered sample folders!!\n Files within Dada_Data Dada_FilteredFastqs and Dada_Plots will be overwritten**")
                }
                
                # exclude Folders starting with "Dada" from the folders considered as sample fodlers
                if(length(grep("^Dada", SampleNames))!=0) {
                        SampleNames <- SampleNames[-grep("^Dada", SampleNames)]
                }
                
                F_fastq <- character(length = length(SampleNames))
                R_fastq <- character(length = length(SampleNames))
                
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
                else {
                        
                        stop("No sample folders were found in the given path! Currently the Dada2_wrap function can only handle the situation where the fastq files are in separate folders for each sample.
                             These folders have to be in the \"path\" folder. Other situations have to be added.")
                }
                
        }
        
        
        ##############################
        ### Determine the quality scores and save the stats in Data folder
        ##############################
        
        DataFolder <- file.path(path2, "Dada_Data")
        
        ## Not sure if this is wanted, deleting all files that are already in the DataFolder folder
        if(file.exists(DataFolder)){
                file.remove(list.files(DataFolder, full.names = TRUE))
        }
        
        dir.create(DataFolder, showWarnings = TRUE)
        
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
        
        ##############################
        ### Generate and save some quality plots
        ##############################
        
        PlotFolder <- file.path(path2, "Dada_Plots")
        
        ## Not sure if this is wanted, deleting all files that are already in the PlotFolder folder
        if(file.exists(PlotFolder)){
                file.remove(list.files(PlotFolder, full.names = TRUE))
        }
        
        dir.create(PlotFolder, showWarnings = FALSE)
        
        
        pdf(file = file.path(PlotFolder, "MedianQScore_F_AllNucleotides.pdf"), width = 7, height = 6)
        print(QS_Median_OverviewPlot(F_QualityStats, SampleNames))
        dev.off()
        
        pdf(file = file.path(PlotFolder, "MedianQScore_F_150to240.pdf"), width = 7, height = 6)
        print(QS_Median_OverviewPlot(F_QualityStats, SampleNames, xlim_low = 150, xlim_high = 240))
        dev.off()
        
        pdf(file = file.path(PlotFolder, "MedianQScore_R_AllNucleotides.pdf"), width = 7, height = 6)
        print(QS_Median_OverviewPlot(R_QualityStats, SampleNames))
        dev.off()
        
        pdf(file = file.path(PlotFolder, "MedianQScore_R_150to240.pdf"), width = 7, height = 6)
        print(QS_Median_OverviewPlot(R_QualityStats, SampleNames, xlim_low = 150, xlim_high = 240))
        dev.off()
        
        message("*********************** Quality plots generated ***********************
                ********************************************************************")
}



