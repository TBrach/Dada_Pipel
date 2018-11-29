path <- "/Users/jvb740/MarieCurie_Work/Project_Normalization/DK_Healthy_Normalization_20181023_Reads/Dada_Analysis/Dada_FilteredFastqs/"
pathsave <- "/Users/jvb740/MarieCurie_Work/Project_Normalization/DK_Healthy_Normalization_20181023_Reads/Dada_Analysis_filteredReads42500/Dada_FilteredFastqs_42500/"


full_files <- dir(path)
#no_reads <- vector(mode = "numeric", length = length(full_files))
min_no_reads <- 42500 # min as determined from ReadSummary

for (i in seq(1, to = length(full_files), by = 2)) {
        current_fast <-  ShortRead::readFastq(file.path(path, full_files[i]))
        indexes <- sort(sample(1:length(current_fast), size = min_no_reads))
        rared_fast <- current_fast[indexes]
        current_fast_r <- ShortRead::readFastq(file.path(path, full_files[i+1]))
        rared_fast_r <- current_fast_r[indexes]
        ShortRead::writeFastq(rared_fast, file.path(pathsave, full_files[i]))
        ShortRead::writeFastq(rared_fast_r, file.path(pathsave, full_files[i+1]))
        
}
