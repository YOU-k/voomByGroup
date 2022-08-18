simscRNAseq4 <- function(mu,phi,subject,correlation,n1)
  # Simulate correlated negative-binomial single-cell RNA-seq counts from multiple subjects
  # Rows are genes and columns are cells.
  # Simulated counts are NB(mu, phi)
  # Counts from the same subject are correlated.
  #
  # Arguments
  # mu = matrix of expected values for the counts, rows are genes and columns are cells
  # phi = negative binomial dispersion. Can be a scalar or gene-vector or matrix.
  # subject = factor or integer vector assigning cells to subjects.
  # correlation = intra-subject correlation. Can be a scalar or gene-vector.
  #
# Gordon Smyth
# Created 15 Jan 2022.
{
  # Check arguments
  NGenes <- nrow(mu)
  NCells <- ncol(mu)
  subject <- as.factor(subject)
  NSubjects <- length(levels(subject))
  
  # Simulate correlated standard normal random variables
  z.cell <- matrix(rnorm(NGenes*NCells),NGenes,NCells)
  z.subject <- matrix(rnorm(NGenes*NSubjects),NGenes,NCells)
  z.subject.cell <- z.subject[,unclass(subject),drop=FALSE]
  z <- sqrt(correlation)*z.subject.cell + z.cell*sqrt(1-correlation)
  
  # Convert from normal to gamma with mean 1 and squared-CV equal to phi
  up <- (z>0)
  dn <- !up
  
  newData <- as.data.frame(log10(rowSums(mu[,1:n1])/4))
  colnames(newData) <- "predictor"
  stats::predict(formula1,newData) -> phi1
  phi1 <- matrix(phi1,NGenes,n1)
  
  
  newData <- as.data.frame(log10(rowSums(mu[,(n1+1):(NCells)])/4))
  colnames(newData) <- "predictor"
  stats::predict(formula2,newData) -> phi2
  phi2 <- matrix(phi2,NGenes,NCells-n1)
  
  cbind(phi1,phi2) -> phi
  z[up] <- qgamma(pnorm(z[up],lower.tail=FALSE,log.p=TRUE),shape=1/phi[up],scale=phi[up],lower.tail=FALSE,log.p=TRUE)
  z[dn] <- qgamma(pnorm(z[dn],lower.tail=TRUE,log.p=TRUE),shape=1/phi[dn],scale=phi[dn],lower.tail=TRUE,log.p=TRUE)
  
  # Expected counts
  mu <- mu * z
  
  # Actual counts
  z <- mu
  z[] <- rpois(NGenes*NCells, lambda=mu)
  z
}
