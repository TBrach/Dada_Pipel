function (D = EucDist, correction = "none", rn = NULL) 
{
        centre <- function(D, n) {
                One <- matrix(1, n, n)
                mat <- diag(n) - One/n
                mat.cen <- mat %*% D %*% mat
        }
        bstick.def <- function(n, tot.var = 1, ...) {
                res <- rev(cumsum(tot.var/n:1)/n)
                names(res) <- paste("Stick", seq(len = n), sep = "")
                return(res)
        }
        D <- as.matrix(D)
        n <- nrow(D) # would be equal to ncol(D)
        epsilon <- sqrt(.Machine$double.eps)
        if (length(rn) != 0) {
                names <- rn
        } else {
                names <- rownames(D)
        }
        CORRECTIONS <- c("none", "lingoes", "cailliez")
        correct <- pmatch(correction, CORRECTIONS)
        if (is.na(correct)) 
                stop("Invalid correction method")
        delta1 <- centre((-0.5 * D^2), n)
        ## understand this a bit better:
        if(is.null(5)) {
                D1 <- -.5*D^2 # NB all values except the 0 are negative
                One <- matrix(1, n, n) # a matrix full of ones
                mat <- diag(n) - One/n # 0.9166 on the diagonal, aotherwise -1/12
                mat.cen <- mat %*% D1 %*% mat # a double matrix multiplication
                
                apply(mat.cen, 2, mean) # basically all 0, so centering puts the mean to 0
                apply(mat.cen, 2, sd) # are different
                # NB: it is not just subtracting the colon means:
                mat.comp <- scale(D1, scale = FALSE)
                mat.comp2 <- scale(D1, scale = TRUE)
                apply(mat.comp, 2, mean) # also all 0 but the values are different
                apply(mat.comp2, 2, mean)
                apply(mat.comp2, 2, sd) # all 1 because scaled
                apply(mat.cen, 1, mean) # also all 0!!!
                # NOTE mat.cen is still fully symmetrical
                all.equal(mat.cen, t(mat.cen))
                mat.cen[4,]/mat.cen[,4]
                # so centre function leaves matrix symmetrical but centers columns (and also rows) in a way that the mean becomes 0.
                # is it a double centering?
                mat.comp <- scale(D1, scale = FALSE)
                mat.comp.comp <- t(scale(t(mat.comp), scale = F))
                all.equal(mat.comp.comp, t(mat.comp.comp)) # TRUE
                # AND INDEED
                all.equal(mat.comp.comp, mat.cen, check.attributes = FALSE) # TRUE!!!
                ## So centre function is just centering with the mean first the columns then the rows
                # BUT NO CLUE WHY NOT centering D and instead .5*D^2
        }
        
        ###
        trace <- sum(diag(delta1))
        D.eig <- eigen(delta1)
        ## NB: here ----
        svddelta1 <- svd(delta1) # singular value decomposition
        all.equal(svddelta1$d, D.eig$values) # TRUE
        all.equal(svddelta1$v, D.eig$vectors) # Not equal but differences only in sign (Vorzeichen), see
        all.equal(abs(svddelta1$v), abs(D.eig$vectors))
        # -------
        
        min.eig <- min(D.eig$values)
        zero.eig <- which(abs(D.eig$values) < epsilon) # here the last three eigen values are basically 0
        D.eig$values[zero.eig] <- 0
        # NB: all.equal(svddelta1$d, D.eig$values) # not True anymore
        if (min.eig > -epsilon) {
                correct <- 1
                eig <- D.eig$values
                k <- length(which(eig > epsilon))
                rel.eig <- eig[1:k]/trace
                cum.eig <- cumsum(rel.eig)
                vectors <- sweep(D.eig$vectors[, 1:k], 2, sqrt(eig[1:k]), 
                                 FUN = "*")
                bs <- bstick.def(k)
                cum.bs <- cumsum(bs)
                res <- data.frame(eig[1:k], rel.eig, bs, cum.eig, cum.bs)
                colnames(res) <- c("Eigenvalues", "Relative_eig", "Broken_stick", 
                                   "Cumul_eig", "Cumul_br_stick")
                rownames(res) <- 1:nrow(res)
                rownames(vectors) <- names
                colnames(vectors) <- colnames(vectors, do.NULL = FALSE, 
                                              prefix = "Axis.")
                note <- paste("There were no negative eigenvalues. No correction was applied")
                out <- (list(correction = c(correction, correct), note = note, 
                             values = res, vectors = vectors, trace = trace))
        }
        else {
                k <- n
                eig <- D.eig$values
                rel.eig <- eig/trace
                rel.eig.cor <- (eig - min.eig)/(trace - (n - 1) * min.eig)
                rel.eig.cor = c(rel.eig.cor[1:(zero.eig[1] - 1)], rel.eig.cor[(zero.eig[1] + 
                                                                                       1):n], 0)
                cum.eig.cor <- cumsum(rel.eig.cor)
                k2 <- length(which(eig > epsilon))
                k3 <- length(which(rel.eig.cor > epsilon))
                vectors <- sweep(D.eig$vectors[, 1:k2], 2, sqrt(eig[1:k2]), 
                                 FUN = "*")
                if ((correct == 2) | (correct == 3)) {
                        if (correct == 2) {
                                c1 <- -min.eig
                                note <- paste("Lingoes correction applied to negative eigenvalues: D' = -0.5*D^2 -", 
                                              c1, ", except diagonal elements")
                                D <- -0.5 * (D^2 + 2 * c1)
                        }
                        else if (correct == 3) {
                                delta2 <- centre((-0.5 * D), n)
                                upper <- cbind(matrix(0, n, n), 2 * delta1)
                                lower <- cbind(-diag(n), -4 * delta2)
                                sp.matrix <- rbind(upper, lower)
                                c2 <- max(Re(eigen(sp.matrix, symmetric = FALSE, 
                                                   only.values = TRUE)$values))
                                note <- paste("Cailliez correction applied to negative eigenvalues: D' = -0.5*(D +", 
                                              c2, ")^2, except diagonal elements")
                                D <- -0.5 * (D + c2)^2
                        }
                        diag(D) <- 0
                        mat.cor <- centre(D, n)
                        toto.cor <- eigen(mat.cor)
                        trace.cor <- sum(diag(mat.cor))
                        min.eig.cor <- min(toto.cor$values)
                        zero.eig.cor <- which((toto.cor$values < epsilon) & 
                                                      (toto.cor$values > -epsilon))
                        toto.cor$values[zero.eig.cor] <- 0
                        if (min.eig.cor > -epsilon) {
                                eig.cor <- toto.cor$values
                                rel.eig.cor <- eig.cor[1:k]/trace.cor
                                cum.eig.cor <- cumsum(rel.eig.cor)
                                k2 <- length(which(eig.cor > epsilon))
                                vectors.cor <- sweep(toto.cor$vectors[, 1:k2], 
                                                     2, sqrt(eig.cor[1:k2]), FUN = "*")
                                bs <- bstick.def(k2)
                                bs <- c(bs, rep(0, (k - k2)))
                                cum.bs <- cumsum(bs)
                        }
                        else {
                                if (correct == 2) 
                                        cat("Problem! Negative eigenvalues are still present after Lingoes", 
                                            "\n")
                                if (correct == 3) 
                                        cat("Problem! Negative eigenvalues are still present after Cailliez", 
                                            "\n")
                                rel.eig.cor <- cum.eig.cor <- bs <- cum.bs <- rep(NA, 
                                                                                  n)
                                vectors.cor <- matrix(NA, n, 2)
                        }
                        res <- data.frame(eig[1:k], eig.cor[1:k], rel.eig.cor, 
                                          bs, cum.eig.cor, cum.bs)
                        colnames(res) <- c("Eigenvalues", "Corr_eig", "Rel_corr_eig", 
                                           "Broken_stick", "Cum_corr_eig", "Cum_br_stick")
                        rownames(res) <- 1:nrow(res)
                        rownames(vectors) <- names
                        colnames(vectors) <- colnames(vectors, do.NULL = FALSE, 
                                                      prefix = "Axis.")
                        out <- (list(correction = c(correction, correct), 
                                     note = note, values = res, vectors = vectors, 
                                     trace = trace, vectors.cor = vectors.cor, trace.cor = trace.cor))
                }
                else {
                        note <- "No correction was applied to the negative eigenvalues"
                        bs <- bstick.def(k3)
                        bs <- c(bs, rep(0, (k - k3)))
                        cum.bs <- cumsum(bs)
                        res <- data.frame(eig[1:k], rel.eig, rel.eig.cor, 
                                          bs, cum.eig.cor, cum.bs)
                        colnames(res) <- c("Eigenvalues", "Relative_eig", 
                                           "Rel_corr_eig", "Broken_stick", "Cum_corr_eig", 
                                           "Cumul_br_stick")
                        rownames(res) <- 1:nrow(res)
                        rownames(vectors) <- names
                        colnames(vectors) <- colnames(vectors, do.NULL = FALSE, 
                                                      prefix = "Axis.")
                        out <- (list(correction = c(correction, correct), 
                                     note = note, values = res, vectors = vectors, 
                                     trace = trace))
                }
        }
        class(out) <- "pcoa"
        out
}


# --------------- Steps from the distance matrix
# Step 1
# Multiply squared distance matrix with -.5
DistMat <-  as.matrix(EucDist)
DM <- -.5*DistMat^2
# Step 2 center DM
# either Center function above or
DMC <- scale(DM, scale = FALSE)
DMC <- t(scale(t(DMC), scale = F))
# or
DMC <- centre(DM, nrow(DM))
# Step 3 Compute eigenvectors
D.eig <- eigen(DMC)
# compare
svd1 <- svd(DMC)
# all.equal(svd1$d, D.eig$values)
all.equal(svd1$v, D.eig$vectors)
all.equal(abs(svd1$v), abs(D.eig$vectors)) # absolute values are the same but some sign differences
all.equal(abs(svd1$u), abs(D.eig$vectors)) # same applies for u here, probably because DMC was double centered
all.equal(svd1$u, D.eig$vectors)
# usually svd1$u the right singular vectors are the principal components, let's check
pc <- prcomp(DMC)
all.equal(svd1$u, pc$rotation, check.attributes = FALSE) # TRUE
# also note 
# would be variance explained by the principal components
svd1$d^2/sum(svd1$d^2))
# play with plot here a bit