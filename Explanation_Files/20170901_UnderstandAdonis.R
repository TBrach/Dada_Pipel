#  -- Adonis and co background to check for differences between groups in a distance matrix --

# - adonis only tells you if the group_var as a whole is important to explain variance
# in the distance matrix -

# NB adonis gives you no pairwise comparisons if you have more than two levels in your group_var factor, so adonis gives you only the overall significance of the group, see pairwise.perm.manova {RVAideMemoire} example 
# <https://www.researchgate.net/post/How_can_I_conduct_multiple_comparisons_when_I_use_veganadonis>
# <https://www.mothur.org/forum/viewtopic.php?t=3897>
# maybe also <http://thebiobucket.blogspot.dk/2011/08/two-way-permanova-adonis-with-custom.html>

require(vegan)
data(iris)

# adonis works with either calculating distance on the go, or 
euclidDist <- vegdist(iris[,1:4], method = "euclidian")
# permutation MANOVA
set.seed(1)
adonis(iris[,1:4]~Species, data=iris, method="euclidian")
set.seed(1)
adonis(formula = euclidDist ~ Species, data = iris)
# fewer permutations gives higher p-values
set.seed(1)
adonis(formula = euclidDist ~ Species, data = iris, permutations = 25)
# adonis strata on species is gives same coefficients but p-value = 1
set.seed(1)
adonis(formula = euclidDist ~ Species, data = iris, strata = iris$Species)
# adonis2 gives same result as adonis but somewhat less info, so stick with adonis, even though it might take longer
# as mentioned in the help
set.seed(1)
adonis2(formula = euclidDist ~ Species, data = iris)

# so adonis said Species has a p-value of 0.001 and explained 86 % of variance (R2)
# there are three Species
length(unique(iris$Species))

# pairwise.perm.manova from RVAideMemoire package can give you the p-values between each level
pairwise.perm.manova(euclidDist, iris$Species, nperm = 49, p.method = "none")
# NB: more permutations = lower p values, just as for adonis
pairwise.perm.manova(euclidDist, iris$Species, nperm = 999, p.method = "none")


# see the code of pairwise.perm.manova, in case of resp = a dist object, it in
# fact loops simply through vegan::adonis!

# based on it I wrote pairwise.perm.manova.own, so see this code


        

     