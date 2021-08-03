#########################
#########################
#### Clear Workspace ####
#########################
#########################

rm(list = ls()) 
# clear global environment to remove all loaded data sets, functions and so on.

###################
###################
#### Libraries ####
###################
###################

library(easypackages) # enables the libraries function
suppressPackageStartupMessages(
  libraries("corpcor",
            "GPArotation",
            "psych"))

#############
#############
# Functions #
#############
#############

#####################################################
# KMO Kaiser-Meyer-Olkin Measure of Sampling Adequacy
#####################################################

# Function by G. Jay Kerns, Ph.D., Youngstown State University (http://tolstoy.newcastle.edu.au/R/e2/help/07/08/22816.html)

kmo = function( data ){
  library(MASS) 
  X <- cor(as.matrix(data)) 
  iX <- ginv(X) 
  S2 <- diag(diag((iX^-1)))
  AIS <- S2%*%iX%*%S2                      # anti-image covariance matrix
  IS <- X+AIS-2*S2                         # image covariance matrix
  Dai <- sqrt(diag(diag(AIS)))
  IR <- ginv(Dai)%*%IS%*%ginv(Dai)         # image correlation matrix
  AIR <- ginv(Dai)%*%AIS%*%ginv(Dai)       # anti-image correlation matrix
  a <- apply((AIR - diag(diag(AIR)))^2, 2, sum)
  AA <- sum(a) 
  b <- apply((X - diag(nrow(X)))^2, 2, sum)
  BB <- sum(b)
  MSA <- b/(b+a)                        # indiv. measures of sampling adequacy
  AIR <- AIR-diag(nrow(AIR))+diag(MSA)  # Examine the anti-image of the correlation matrix. That is the  negative of the partial correlations, partialling out all other variables.
  kmo <- BB/(AA+BB)                     # overall KMO statistic
  # Reporting the conclusion 
  if (kmo >= 0.00 && kmo < 0.50){test <- 'The KMO test yields a degree of common variance unacceptable for FA.'} 
  else if (kmo >= 0.50 && kmo < 0.60){test <- 'The KMO test yields a degree of common variance miserable.'} 
  else if (kmo >= 0.60 && kmo < 0.70){test <- 'The KMO test yields a degree of common variance mediocre.'} 
  else if (kmo >= 0.70 && kmo < 0.80){test <- 'The KMO test yields a degree of common variance middling.' } 
  else if (kmo >= 0.80 && kmo < 0.90){test <- 'The KMO test yields a degree of common variance meritorious.' }
  else { test <- 'The KMO test yields a degree of common variance marvelous.' }
  
  ans <- list( overall = kmo,
               report = test,
               individual = MSA,
               AIS = AIS,
               AIR = AIR )
  return(ans)
} 

##################
# Residual Stats #
##################

# Includes all the factor analysis residual commands into one function

residual.stats <- function(matrix){
  residuals <- as.matrix(matrix[upper.tri(matrix)])
  large.resid <- abs(residuals) > 0.05
  numberLargeResids <- sum(large.resid)
  propLargeResid <- numberLargeResids/nrow(residuals)
  rmsr <- sqrt(mean(residuals^2))
  
  cat("Root means squared residual = ", rmsr, "\n")
  cat("Number of absolute residuals > 0.05 = ", numberLargeResids, "\n")
  cat("Proportion of absolute residuals > 0.05 = ", propLargeResid, "\n")
  hist(residuals)
}


setwd("C:/R Portfolio/Factor Analysis")
raqData <- read.delim("raq.dat", header = T)


# create a correlation matrix
raqMatrix <- cor(raqData)
raqMatrix
round(raqMatrix, 2)

###############
###############
# Assumptions #
###############
###############

# For factor analysis to be appropriate, relationships between variables must be evident, patterns of correlation must be compact and multicollinearity must not be evident. Tests for each assumption are completed next.

############################################
############################################
# Test for relationships between variables #
############################################
############################################

# For factor analysis to work, some relationships between variables must be evident and this can be tested using the bartlett test. A significant result means that there are some relationships between variables intended to be used in the analysis. The output is significant, and so, factor analysis is appropriate.

# Bartlett's test

cortest.bartlett(raqData)
cortest.bartlett(raqMatrix, n = 2571)

####################################
####################################
# Compact Patterns of Correlations #
####################################
####################################

# The purpose of this test is to indicate adequacy of sample size and data set for factor analysis. Generally, the following is a rule of thumb:

# A value below .5 should cause the researcher to think again regarding suitability of factor analysis. Here, more data could be collected or rethink which variables to include.
# Values between .5 and .7 are mediocre. between .7 and .8 are good, between .8 and .9 are great and above .9 are superb. The highest value possible is 1

# This can be tested using the Kaiser-Meyer-Olkin test, as follows:

kmo(raqData)

# The overall value is 0.9302245, so we can be confident the sample size and data are adequate for factor analysis

#####################
#####################
# Multicollinearity #
#####################
#####################

# Determinants of the Correlation Matrix

# Here, the critical value is 0.00001, above which there is nothing problematic in the determinant. Our value is 0.0005271037, which is above this. Therefore, multicollinearity is not a problem

det(raqMatrix)

# Alternatively 

det(cor(raqData))

#####################
#####################
# Factor Extraction #
#####################
#####################

# We use principal components analysis to extract factors, but pca generates results similar to factor analysis

pca <- principal(raqData, 
                 nfactors = length(raqData), # use data length as factor number
                 rotate = "none")
pca

# with explicit factor analysis, the number of factors must be lower than the number of variables.

# in the pca matrix, the column h2 refers to all the communalities(sometimes called h2). These communalities are equal to 1 as we have extracted 23 factors. When fewer factors (or components) have been extracted, lower communalities will be evident. 
# The unique variance column is headed by u2. This refers to the unique variance or communality for each variable and it is calculated by subtracting the h2 communality from 1. So, all our u2 values are zero.
# The eigenvalues are values associated with each factor that represent  the variance explained by the particular linear of the squared loadings. They are called sum of square loadings  and can be accessed by using the term values attathed to each component or factor name, for instance pca$values list th values for each component or factor. The proportion of variance row explains the proportion that each factor or principal component explains. For instance, pc1 explains 7.29 units out of 23, which is 32% of the total variance.

##############
# Scree Plot #
##############

plot(pca$values, 
     type = "b") # indicates the line and the dot for each factor or component along the x-axis.
# The plot begins to tail off after three factors and then there is another drop after four factors, before it begins to plateau afterwards. Therefore, between 2 - four factors could be plausibly retained.Given the large sample, it is probably safe to assume Kaiser's criterion. 
# The scree plot an eigenvalues evidence suggests a 4-component may be best
# Therefore, w can rerun the analysis with the optimal number of factors (4 in this case) in mind to be extracted.

###########################
###########################
# Factor Extraction Rerun #
###########################
###########################

# both lines of code below produce the same results

pca2 <- principal(raqData, nfactors = 4, rotate = "none")
pca2 <- principal(raqMatrix, nfactors = 4, rotate = "none")

pca2
# The only changes from the original pca model are found in the h2 and u2 columns, The SS loadings and proportion of variances remain the same.
# When factors have been extracted, we can see that question 1 explains 43% of the common variance. When certain factors have been extracted, we can get a better idea of how much variance is actually common.
# We can return to Kaisers criterion to determine if such factor extraction is justified. This criterion is accurate when there are fewer than 30 variables and communalities after extraction that are greater than .07 or when sample size exceeds 250 rows and the average communality is greater than .6.
# In our example, only one h2 value is above .7 and the average is .503, so on both counts, the Kaiser criterion does not justify extraction of the four factors. Yet, Kaiser criterion is based on a much smaller sample size than we have and the scree plot indicated 4 factors would be a good threshold.

# Compare with two factors #

pca3 <- principal(raqData, nfactors = 2, rotate = "none")
pca3

# Performance against the Kaiser criterion is even worse than with 4 factors.

##############################################################################
# Compare reproduced correlation matrix against correlation matrix in the data
##############################################################################

factor.model(pca2$loadings)

# The diagonal of the matrix contains the communalities after extraction for each variable. The diagonal is the uniqueness.

# When we subtract all the differences between the observed residuals from the model residuals down the diagonal, the final figure is 0.96 and any figure above 0.95 is a good fit, which indicates that four factors is sufficient.

# Another way to indicate goodness of fit using residuals is determine the count of residuals with absolute values greater than 0.05

# Here, it is first necessary to create an object called residuals that contain all the factor residuals

residuals <- factor.residuals(raqMatrix, pca2$loadings)

# Then extract the upper triangle using the upper.tri() function
# This extracts only the elements above the diagonal (so the elements below are disregared)
residuals <- as.matrix(residuals[upper.tri(residuals)])

# residuals with absolute count greater than 0.05

large.resid <- abs(residuals) > 0.05

# sum the large residuals

sum(large.resid)


# proportion of the total number of residuals

sum(large.resid) / nrow(residuals)
residuals 
# It is 36% and any figure greater than 50% is a cause for worry, so we are fine

# Alternatively, we can also explore the residuals throug the sqrt of the mean

sqrt(mean(residuals ^ 2))

# The figure is 0.055, which as it remains under 0.08 is still okay.

# Finally, there is also te histogram

hist(residuals)
# The plot indicates approximate normality

resids <- factor.residuals(raqMatrix, pca2$loadings)

residual.stats(resids)

# Feed this straight into the residual.stats function

residual.stats(factor.residuals(raqMatrix, pca2$loadings))

############
############
# Rotation #
############
############

# Rotation makes it clearer which variables relate to which factors
# Rotation changes factors to distribute variance differently,  but it cannot account for more or less variance in the variables than it could before rotation
# Therefore, the values in the h2 and u2 columns will be the same.
# If factors are independent (unrelated), then one of the orthogonal rotations should be selected. Perhaps varimax
# If factors correlate, then one of the oblique should be chosen, maybe oblimin or promax

#################################
# Orthogonal Rotation - Varimax #
#################################

pca3 <- principal(raqData, nfactors = 4, rotate = "varimax")

print.psych(pca3, sort = T)

# The results indicate that there are four subscales in the initial questionnaire: fear of statistics, fear of maths and fear of negative peer review

####################
# Oblique Rotation #
####################

# In our context, this is perhaps the better choice of rotation as the variables are correlated.

pca4 <- principal(raqData, nfactors = 4, rotate = "oblimin")
print.psych(pca4, sort = T)

# Factor one represents fear of computers, factor two represents fear of peer evaluation, factor 3 indicates fear of statistics and factor represents fear of mathematics.
# The correlation table in the output indicates correlation, which means we cant assume independence and therefore, the orthogonal rotation should not be trusted. Instead, the oblique rotated solution is perhaps more meaningful.