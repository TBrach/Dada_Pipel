path <- "/Users/jvb740/MarieCurie_Work/MouseProject/ResultsAndProtocols/ManiAging_Results/16S_Sequencing/2017-07-13_DK_age_ManiAging/Dada2_Analysis/Dada_FilteredFastqs"
pathsave <- "/Users/jvb740/MarieCurie_Work/MouseProject/ResultsAndProtocols/ManiAging_Results/16S_Sequencing/2017-07-13_DK_age_ManiAging/Dada2_Analysis/Dada_FilteredFastqs_rarified"


full_files <- dir(path)
#no_reads <- vector(mode = "numeric", length = length(full_files))
min_no_reads <- 40389

for (i in 1:length(full_files)) {
        current_fast <-  ShortRead::readFastq(file.path(path, full_files[i]))
        indexes <- sample(1:length(current_fast), size = min_no_reads)
        rared_fast <- current_fast[indexes]
        ShortRead::writeFastq(rared_fast, file.path(pathsave, full_files[i]))
        
}
