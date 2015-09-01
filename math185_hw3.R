# Pr1 --------------------------------------------
ComputeFrac <- function(n, B = 10000, sample.dist, lambda){
	# Compute the fraction of times the confidence interval contains the true mean 
	# for different population distribution.
	#
	# Args: 
	#		n: sample size
	# 		B: repetition
	# 		sample.dist: 1--standard normal, 2--double exponential, 3--exponential
	#		lambda: parameter for (double) exponential distribution.
	#
	# Returns: frac: fraction.
	
	
	if(!(sample.dist %in% c(1, 2, 3))){
		stop("Please specify a sample distribution.")
	}
	
	if(sample.dist == 1){
		conf.itv <- matrix(NA, B, 2) # initialize a matrix to store conf itv.
		for(i in 1:B){
			x.sample <- rnorm(n) 
			x.mean   <- mean(x.sample)
			x.error  <- qt(0.975, df = n-1) * sd(x.sample) / sqrt(n) # posiive number
			conf.itv[i, ] <- c(x.mean - x.error, x.mean + x.error)
		}
		containing.true <- sum( (conf.itv[, 1] < 0) & (conf.itv[, 2] > 0) )
    # count the number of conf itvs containing the true mean--0
		frac <- containing.true / B
		return(frac)
	} else if(sample.dist == 2){
		conf.itv <- matrix(NA, B, 2)
		for(i in 1:B){
			x.sample <- sample(c(-1, 1), n, TRUE) * rexp(n, rate = lambda)
			x.mean   <- mean(x.sample)
			x.error  <- qt(0.975, df = n-1) * sd(x.sample) / sqrt(n)
			conf.itv[i, ] <- c(x.mean - x.error, x.mean + x.error)
		}
		containing.true <- sum( (conf.itv[, 1] < 0) & (conf.itv[, 2] > 0) )
		frac <- containing.true / B
		return(frac)
	} else {
		conf.itv <- matrix(NA, B, 2)
		for(i in 1:B){
			x.sample <- rexp(n, rate = lambda)
			x.mean   <- mean(x.sample)
			x.error  <- qt(0.975, df = n-1) * sd(x.sample) / sqrt(n)
			conf.itv[i, ] <- c(x.mean - x.error, x.mean + x.error)
		}
		containing.true <- sum( (conf.itv[, 1] < 1) & (conf.itv[, 2] > 1) )
		frac <- containing.true / B
		return(frac)
	}
}

# Part a: 
ComputeFrac(n = 15, sample.dist = 1)
ComputeFrac(n = 150, sample.dist = 1)
# The fraction is close to 0.95, as it should be, and increasing sample size won't make a difference.  This is expected since the level is _exact_ no matter what the sample size n is.

# part b
ComputeFrac(n = 15, sample.dist = 2, lambda = sqrt(2))
ComputeFrac(n = 150, sample.dist = 2, lambda = sqrt(2))
# The fraction is close to 0.95, despite the fact that the distribution is not normal.  This is because of the central limit theorem and the fact that the double-exponential distribution has no skewness (b/c it is symmetric) and the tails are fairly light (the decay is exponential).

# part c
ComputeFrac(n = 15, sample.dist = 3, lambda = 1)
ComputeFrac(n = 150, sample.dist = 3, lambda = 1)
# The fraction is close to 0.95 for n = 150, showing the central limit theorem at play again.  However, the fraction is off when n = 15, the central limit theorem has not "kicked in" yet.  The difference here is that the exponential distribution has substantial skewness (b/c it is asymmetric).



# Pr2 --------------------------------------------
# part a ----------------------
sd.bootCI.percentile <- function(x, level = 0.95, B = 200){
	# Compute bootstrap percentile conf itv for population standard deviation.
	#
	# Args: 
	#		   x: input sample
	# 		 level: level of conf itv
	# 		 B: number of repetitions
	#
	# Returns: sd.bootCI.perc: boot percentile conf itv for population std dev
	
	alpha   <- 1 - level
	sd.boot <- numeric(B)
	for(i in 1:B){
		x.boot <- sample(x, length(x), replace = TRUE)
		sd.boot[i] <- sd(x.boot)
	}
	sd.bootCI.perc <- quantile(sd.boot, c(alpha / 2, 1 - alpha / 2), names = FALSE)
	return(sd.bootCI.perc)
}

# part b ----------------------
sd.bootCI.pivotal <- function(x, level = 0.95, B = 200){
	# Compute bootstrap pivotal conf itv for population standard deviation.
	#
	# Args: 
	#		   x: input sample
	# 		 level: level of conf itv
	# 		 B: number of repetitions
	#
	# Returns: sd.bootCI.pivo: boot pivotal conf itv for population std dev.
	
	alpha   <- 1 - level
	sd.hat  <- sd(x)
	sd.boot <- numeric(B)
	for(i in 1:B){
		x.boot <- sample(x, length(x), replace = TRUE)
		sd.boot[i] <- sd(x.boot)
	}
	boot.quantile  <- quantile(sd.boot, c(1 - alpha / 2, alpha / 2), names = FALSE)
	sd.bootCI.pivo <- rep(2 * sd.hat, 2) - boot.quantile
	return(sd.bootCI.pivo)
}
	
# part c ----------------------
CompareAB <- function(N = 100, level = 0.95, B = 200, rep = 1000){
	# Compare methods A and B by computing the fraction that boot conf itv 
	# contains the true standard deviation under standard normal sample.
	# 
	# Args: 
	# 	    N: sample size
	#		level
	#		B: number of bootstrap samples (repetitions)
	#		rep: how many times do we repete the comparison procedure.
	#
	# Returns: fraction of each method.

	conf.itv.perc <- matrix(NA, rep, 2)
	conf.itv.pivo <- matrix(NA, rep, 2)
	for(i in 1:rep){
		x <- rnorm(N)
		conf.itv.perc[i, ] <- sd.bootCI.percentile(x, level, B) 
		conf.itv.pivo[i, ] <- sd.bootCI.pivotal(x, level, B)
	}
	contain.true.perc <- sum( (conf.itv.perc[, 1] < 1) & (conf.itv.perc[, 2] > 1) )
	frac.perc <- contain.true.perc / rep
	contain.true.pivo <- sum( (conf.itv.pivo[, 1] < 1) & (conf.itv.pivo[, 2] > 1) )
	frac.pivo <- contain.true.pivo / rep

	return(c(frac.perc, frac.pivo))
}	

set.seed(185)
CompareAB()	# Run the comparison
            # [1] 0.919 0.923
	
# Note that with B = 2000 an rep = 10000 (which takes a while to run) the results with the same seed are 0.9303 and 0.9355.  So the coverage for the pivotal method seems slightly better.  However, this difference is not statistically significant!  How would you test that? 	
	
	

	
