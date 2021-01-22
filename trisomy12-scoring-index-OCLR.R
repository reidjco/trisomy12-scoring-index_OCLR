# Human pluripotent stem cells identify molecular targets of trisomy 12 in chronic lymphocytic leukemia patients
# Jennifer C. Reid et al. Cell Reports, 2021

# Load packages
library(gelnet)
library(dplyr)

# TRAINING DATASET ANALYSIS
X <- read.delim("Rmatrix_training.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE, row.names=1)
X <- as.matrix(X)
m <- apply ( X, 1, mean )
X <- X - m

# Identify the tri12 samples and break up all samples into 2 groups
X.tr <- X[,1:3]	# assign the tri12 columns as the training dataset
X.bk <- X[,4:5]	# assign the “not tri12” columns as the background dataset
mm <- gelnet( t(X.tr), NULL, 0, 1 )

# Store the signature to a file, referred to as Table S6 in manuscript
write.table(mm$w, file = "Tri12sig.txt", sep = "\t", quote = FALSE, col.names = FALSE)

# Perform leave-one-out cross-validation
auc <- c()
for( i in 1:ncol(X.tr) )
{
  
  X1 <- X.tr[,-i]
  m1 <- gelnet( t(X1), NULL, 0, 1 )
  
  sbk <- apply( X.bk, 2, function(z) {cor( m1$w, z, method="sp" )} )
  s1 <- cor( m1$w, X.tr[,i], method="sp" )
  
  auc[i] <- sum( s1 > sbk ) / length(sbk)
  cat( "Current AUC: ", auc[i], "\n" )
  cat( "Average AUC: ", mean(auc), "\n" )
}
# Tri12 scoring index is provided as Table S6 in Reid et al. Cell Reports 2021 and also in this repository

# TESTING DATASET ANALYSIS
# predict Tri12 Score of 159 patient samples
# "Tri12sig.txt file refers to Table S5
w <- read.delim("Tri12sig.txt", sep="\t", header=FALSE, stringsAsFactors = FALSE, row.names=1)
w <- as.matrix(w)
w <- drop(w)
Z <- read.table("Rmatrix_testing.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE, row.names=1)
Z <- as.matrix(Z)
s <- apply( Z, 2, function(z) {cor( z, w, method="sp", use="complete.obs" )} )

# Scale the scores to be between 0 and 1 and then output the Tri12 Scores of the 159 patients
s <- s - min(s)
s <- s / max(s)
write.table( cbind(s), file="Tri12score_patients.txt", sep="\t", quote=FALSE, col.names=TRUE )
#This data is provided as part of Table S4 in Reid et al. Cell Reports 2021