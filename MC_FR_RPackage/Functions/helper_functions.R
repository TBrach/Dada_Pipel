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