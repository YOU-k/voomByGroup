simscRNAseq_group4 <- function(mu,phi,subject,correlation,group)
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
  n1=table(group)[1]
  n2=table(group)[2]
  n3=table(group)[3]
  n4=table(group)[4]
  
  # Simulate correlated standard normal random variables
  z.cell <- matrix(rnorm(NGenes*NCells),NGenes,NCells)
  z.subject <- matrix(rnorm(NGenes*NSubjects),NGenes,NCells)
  z.subject.cell <- z.subject[,unclass(subject),drop=FALSE]
  z <- sqrt(correlation)*z.subject.cell + z.cell*sqrt(1-correlation)
  
  # Convert from normal to gamma with mean 1 and squared-CV equal to phi
  up <- (z>0)
  dn <- !up
  
  phi1 <- matrix(phi[1],NGenes,n1)
  phi2 <- matrix(phi[2],NGenes,n2)
  phi3 <- matrix(phi[3],NGenes,n3)
  phi4 <- matrix(phi[4],NGenes,n4)
  
  cbind(phi1,phi2,phi3,phi4) -> phi
  z[up] <- qgamma(pnorm(z[up],lower.tail=FALSE,log.p=TRUE),shape=1/phi[up],scale=phi[up],lower.tail=FALSE,log.p=TRUE)
  z[dn] <- qgamma(pnorm(z[dn],lower.tail=TRUE,log.p=TRUE),shape=1/phi[dn],scale=phi[dn],lower.tail=TRUE,log.p=TRUE)
  
  # Expected counts
  mu <- mu * z
  
  # Actual counts
  z <- mu
  z[] <- rpois(NGenes*NCells, lambda=mu)
  z
}
