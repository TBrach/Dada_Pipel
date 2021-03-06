# --
###########################
### FUNCTION: areColors ###
############################
# Found here: #https://stackoverflow.com/questions/13289009/check-if-character-string-is-a-valid-color-representation?utm_medium=organic&utm_source=google_rich_qa&utm_campaign=google_rich_qa
# is used to test color_levels input, i.e. tests if each entry in color_levels is a character
# representation of a color

areColors <- function(x) {
        sapply(x, function(X) {
                tryCatch(is.matrix(grDevices::col2rgb(X)), 
                         error = function(e) FALSE)
        })
}
# --



# --
#######################################
### lmp (to get p_value from lm fit)
#######################################
# needed to get the p-value from-linear fit objects (from stackoverflow)

lmp <- function (modelobject) {
        if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
        f <- summary(modelobject)$fstatistic
        p <- pf(f[1],f[2],f[3],lower.tail=F)
        attributes(p) <- NULL
        return(p)
}
# --


# --
#######################################
### gm_own: calculate geometric mean
#######################################
# see commented below, this function comes from 
# <http://stackoverflow.com/questions/2602583/geometric-mean-is-there-a-built-in>

## Input:
# x numeric vector
# na.rm: if FALSE you get NA as soon as an NA is in your data, if TRUE the NA get basically treated as 0 (but NOTE these zeros always count also when zero.count = FALSE)
# zeros.count: This is IMPORTANT, if TRUE 0s count, so if x contains a 0 the GM will be lower than when zero.count = FALSE, in 
# which case the gm is calculated for the same vector in which all 0 have been removed.
## Output:
# the geometric mean of x
gm_own = function(x, na.rm=FALSE, zeros.count = TRUE){
        if(any(x < 0, na.rm = TRUE)){
                return(NaN)
        }
        if(zeros.count){
                exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
        } else {
                exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x[x > 0]))
        }
}
# --



# --
#######################################
### get_unique_facLevel_combinations
#######################################
# often needed to set the comparisons attribute in stat_compare_means from ggplot, it excludes comparison with each other and also removes 2 to 3 vs 3 to 2

get_unique_facLevel_combinations <- function(fac_levels){
        
        fac_levels_num <- setNames(seq_along(fac_levels), fac_levels) 
        i_s <- outer(fac_levels_num, fac_levels_num, function(ivec, jvec){
                sapply(seq_along(jvec), function(x){
                        i <- jvec[x]
                })
        })
        j_s <- outer(fac_levels_num, fac_levels_num, function(ivec, jvec){
                sapply(seq_along(ivec), function(x){
                        j <- ivec[x]
                })
        })
        i_s <- i_s[lower.tri(i_s)]
        j_s <- j_s[lower.tri(j_s)]
        
        comparisonList <- list()
        for (k in seq_along(i_s)){
                comparisonList[[k]] <- c(fac_levels[i_s[k]], fac_levels[j_s[k]])
        }
        
        comparisonList
        
}
# --


# --
#######################################
### get_unique_facLevel_combis
#######################################
# as get_unique_facLevel_combinations but outputting i_s and j_s

get_unique_facLevel_combis <- function(fac_levels){
        
        fac_levels_num <- setNames(seq_along(fac_levels), fac_levels) 
        i_s <- outer(fac_levels_num, fac_levels_num, function(ivec, jvec){
                sapply(seq_along(jvec), function(x){
                        i <- jvec[x]
                })
        })
        j_s <- outer(fac_levels_num, fac_levels_num, function(ivec, jvec){
                sapply(seq_along(ivec), function(x){
                        j <- ivec[x]
                })
        })
        i_s <- i_s[lower.tri(i_s)]
        j_s <- j_s[lower.tri(j_s)]
        
        list(i_s = i_s, j_s = j_s)
        
}
# --



