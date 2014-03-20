luke.skew.test <- function(p.values, reps = 100000)
{
  p.values <- p.values[which(p.values > 0.03 & p.values < 0.05)]
  n <- length(p.values)

  observed.skew <- (mean(p.values) - median(p.values)) / sd(p.values)
 
 # Now find distribution of the expected skew, giving sampling error with n samples
  rand <- matrix(runif(reps*n, min = 0.03, max = 0.05), ncol = n, nrow = reps) # fully vectorised!
  means <- rowMeans(rand)
  medians <- apply(rand,1,median)
  sds <- apply(rand, 1, sd)
  boot.skew <- (means - medians) / sds
  quantiles <- quantile(boot.skew, probs = c(0.025, 0.975))
 
 # One-tailed p value
  p.null <- length(boot.skew[boot.skew < observed.skew]) / reps

  print(paste("Expected skew is 0.5, with 95% CIs of ", quantiles[1], " to ", quantiles[2], ", given that we have ", n, "p values in the range 0.03-0.05.", sep=""))
  print("Observed skew is ", observed.skew, ", which has a probability of ", p.null, "under the null hypothesis of no p-hacking, with an effect size of zero." sep="")
  
  # To do: add a histogram of the boot values with a little arrow on it showing the observed skew level, like this: http://stackoverflow.com/questions/11122002/how-do-i-draw-an-arrow-on-a-histogram-drawn-using-ggplot2
}

# New stuff added on 20 March. Sorry about the mess! 

\name{luke.skew.test2}
\alias{luke.skew.test2}
\title{Test whether p values in the range 0.03-0.05 are clustered near 0.05.}
\usage{
  luke.skew.test2(p, reps)
}
\arguments{
  \item{p}{A vector of p values. Only those in the range 0.03-0.05 inclusive will be used.}
  \item{reps}{Number of replicates to use when generating the null distribution. Defaults to 10,000.}
  \item{quantiles}{A vector of length 2: quantiles to be returned for the null distribution. Defaults to 0.025 and 0.975, i.e. 95% of the data.}
}
\value{
  A dataframe giving the observed and expected mean and median of the p values, as well as the specified quantiles of the latter. The respective one-tailed p values are also given.
}
\description{
  This function tests whether the mean and the median are significantly higher than predicted under the null hypothesis that the true average effect size underlying the p values in the sample is zero, and there is no p-hacking (sensu Simonsohn et al 2013). Under the this null model, the distribution of p values from 0-1 should be flat, and hence the mean and median of the p values in the range 0.03-0.05 should be 0.04. The function computes the expected distribution of means and medians under the null hypothesis, and calculates a one-tailed p that the observed mean and median are higher than expected by chance (indicating p hacking).
}
\examples{
  p <- c(0.00001, 0.024, 0.002, 0.045, 0.00003, 0.021, 0.0001, 0.0000948, 0.0000002)
  luke.skew.test2(p)
}
\references{
  Simonsohn, Uri and Nelson, Leif D. and Simmons, Joseph P.,
  P-Curve: A Key to the File Drawer (April 24, 2013). Journal
  of Experimental Psychology: General, Forthcoming. Available
  at SSRN: http://ssrn.com/abstract=2256237
}
\seealso{
  \code{\link{binomial.bias.test}},
  \code{\link{binomial.sns.test}}
}
\keyword{curve}
\keyword{p}

luke.skew.test2 <- function(p, reps = 10000, quantiles = c(0.025, 0.975))
{
    n <- length(p)
  
    if(n == 0) 
      {
          print("Argument 'p' has length of zero; enter multiple p values for best results.")
          return(NA)
      }
    
    if(reps == 0) 
      {
          print("Argument 'reps' is zero; must be 1 or greater (numbers 10,000 and above give good accuracy).")
          return(NA)
      }
    
    if(length(reps) != 1) 
    {
      print("Argument 'reps' must have a length of 1.")
            return(NA)
    }
    
    if(length(quantiles) != 2) 
    {
      print("Argument 'quantiles' must have a length of 2.")
            return(NA)
    }
    
    if(min(quantiles) <0 || max(quantiles) > 1) 
    {
      print("Argument 'quantiles' must be a vector of length 2, in the range 0-1.")
      return(NA)
    }
    
    p <- p[which(p >= 0.03 & p <= 0.05)]
  
    observed.mean <- mean(p) # Should be 0.04 or lower under the null hypothesis of no p-hacking
    observed.median <- median(p) # Should be 0.04 or lower under the null hypothesis of no p-hacking
    
    # Now find distribution of the expected means and medians, given sampling error with n samples
    rand <- matrix(runif(reps*n, min = 0.03, max = 0.05), ncol = n, nrow = reps)
    means <- apply(rand, 1, mean) #luke.skew.test2.means
    medians <- apply(rand, 1, median) #luke.skew.test2.medians
    
    # Find the requested quantiles for the null distribution of means and medians
    mean.quantiles <- quantile(means, probs = quantiles)
    median.quantiles <- quantile(medians, probs = quantiles)
    
    # One-tailed p value
    p.null.mean <- length(means[means > observed.mean]) / reps
    p.null.median <- length(medians[means > observed.median]) / reps
    
    # Write the output dataframe
    output <- data.frame(observed = c(observed.mean,observed.median), expected = c(0.04,0.04), lowerCIs = c(mean.quantiles[1],median.quantiles[1]), upperCIs =  c(mean.quantiles[2],median.quantiles[2]), one.tailed.p = c(p.null.mean,p.null.median))
    rownames(output) <- c("Mean", "Median")
    return(output)
}
##### End of new stuff

We can use a similar approach to compare the skew of two distributions. Basically resample the shit out of the two sets of p value, find the difference in skew each time, and then see if the 95% CIs of the distribution of bootstrapped differences includes zero:

luke.compare.skew <- function(p1, p2, reps = 100000)
{
  p1 <- p1[which(p1 > 0.03 & p1 < 0.05)]
  p2 <- p2[which(p2 > 0.03 & p2 < 0.05)]
  
  n1 <- length(p1)
  n2 <- length(p2)
  
  rand1 <- sample(n1,n1*reps,replace=T)
  rand2 <- sample(n2,n2*reps,replace=T)
  
  rand1 <- matrix(n1[rand1], ncol = n1, nrow = reps) # resample the shit out of p1 and p2
  rand2 <- matrix(n2[rand2], ncol = n2, nrow = reps)
  
  means1 <- rowMeans(rand1)                # Find the skew of each boot for both p1 and p2
  medians1 <- apply(rand1,1,median)
  sds1 <- apply(rand1, 1, sd)
  boot.skew1 <- (means1 - medians1) / sds1
  
  means2 <- rowMeans(rand2)
  medians2 <- apply(rand2,1,median)
  sds2 <- apply(rand2, 1, sd)
  boot.skew2 <- (means2 - medians2) / sds2
  
  difference.in.skew <- boot.skew1 - boot.skew2           # Find the difference. It should be zero on average if they have the same skew
  quantiles <- quantile(difference.in.skew, probs = c(0.025, 0.975))   # 95% CIs on the difference. Should not overlap zero if there is a difference
  
  print(paste("The bootstrapped difference in skew between p curves 1 and 2 is ", difference.in.skew, ", with 95% CIs of ", quantiles[1], "-", quantiles[2], ".", sep=""))
}
