


# Pr1 ------------------------------------------------
norm.ks.test <- function(x, B = 2000){
	# Returns p value for the KS test for normality
	# 
	# Args: x: a sample 
	#       B: repetitions
	#
	# Returns: p.value
	
	n <- length(x)
	mu.hat <- mean(x)
	sigma.hat <- sqrt(sum((x - mean(x))^2)/n)
	
	f.hat <- ecdf(x)
	g.hat <- pnorm(x, mean = mu.hat, sd = sigma.hat)
	D <- max(abs(f.hat(x) - g.hat))
	
	D.boot <- numeric(B)
	for(b in 1:B){
		x.boot <- rnorm(n, mean = mu.hat, sd = sigma.hat)
		# parameters are from original sample
		mu.hat.boot <- mean(x.boot)
		sigma.hat.boot <- sqrt(sum((x.boot - mean(x.boot))^2)/n)
	
		f.hat.boot <- ecdf(x.boot)
		g.hat.boot <- pnorm(x.boot, mean = mu.hat.boot, sd = sigma.hat.boot) 
    	# parameters are based on boot sample
		D.boot[b]  <- max(abs(f.hat.boot(x.boot) - g.hat.boot))
	}
	
	R <- sum(D.boot < D)
	p.value <- 1 - R / (B + 1)	
	return(p.value)
}



require(nortest)
# Comparison A: standard normal sample
set.seed(185)
x <- rnorm(500)
norm.ks.test(x)
lillie.test(x)$p.value
# Both yield large p values

# Comparison B: double-exponential sample
set.seed(185)
x <- sample (c (-1, 1), 500, TRUE) * rexp (500, rate = 1)
norm.ks.test(x)
lillie.test(x)$p.value
# lillie test yields smaller p value. 
# However, larger B yields smaller p value in norm.ks.test()
norm.ks.test(x, B = 20000)
# In fact, p value gets smaller as B increases (up to a point).

# Comparison B: exponential sample
set.seed(185)
x <- rexp (500, rate = 1)
norm.ks.test(x)
lillie.test(x)$p.value
# same observation as above.

# NOTE: It happens that the Lilliefors test uses the sample standard deviation, instead of the MLE for the standard deviation, which explains the fact that the statistic that ks.test() and lillie.test() return are slightly different.


# Pr2 ------------------------------------------------
# Set FUN = "unif" below to get unif.cm.test().

cm.test <- function(x, FUN, B = 2000, ...){
	# Return p value for Cramer-von Mises test
	#
	# Args:    x: a sample
	#		 FUN: character string, name for distribution function. 
	#              eg: "norm", "unif", "exp", etc.
	#		   B: number of repetitions
	#       ... : additional arguments to pass to FUN
	#
	# Returns: p value
	
	pFUN <- get(paste("p", FUN, sep = ""), mode = "function", envir = parent.frame())
	# assemble letter "p" with FUN using paste(). e.g. "pnorm"
	# then convert assembled character to type closure (function name) using get(). e.g. pnorm	
	n <- length(x)
	
	first.term <- ((2 * seq(1, n)) - 1) / (2 * n) 
	# store it in a vector to make computation vectorized.
	D <- sqrt(((1 / (12 * n)) + sum((first.term - pFUN(sort(x), ...))^2)) / n)
	# the argument of pFUN is the sorted x's.
	
	rFUN <- get(paste("r", FUN, sep = ""), mode = "function", envir = parent.frame())
	# assemble letter "r" with FUN using paste(). e.g. "rnorm"
	# and convert assembled character to closure (function name) using get(). e.g. rnorm
	D.mc <- numeric(B) 
	# initialize vector to store monte carlo simulated statistics
	for(b in 1:B){
		x.mc <- rFUN(n, ...) 
		# generate random sample with rFUN. "..." are additional parameters for rFUN.
		D.mc[b] <- sqrt(((1 / (12 * n)) +	sum((first.term - pFUN(sort(x.mc), ...))^2)) / n)	
	}
	
	R <- sum(D.mc < D)
	p.value <- 1 - R / (B + 1)	
	return(list(pval = p.value, D = D, R = R, Summary.D.mc = summary(D.mc)))
}


# Comparison with ks.test

# uniform sample
set.seed(185)
x = runif(500)

ks.test(x = x, y = punif)$p.value
cm.test(x, "unif")

ks.test(x = x, y = pnorm)$p.value
cm.test(x, "norm")
cm.test(x, "norm", B = 10000)

# normal sample
set.seed(185)
x = rnorm(500)

ks.test(x = x, y = punif)$p.value
cm.test(x, "unif")
cm.test(x, "unif", B = 10000)

ks.test(x = x, y = pnorm)$p.value
cm.test(x, "norm")

# The performance of cm.test depends on B.



# Pr3 ------------------------------------------------
# Load the data 
dat <- dat[, 3:5]
dat <- matrix(as.numeric(as.character(as.matrix(dat))), ncol = 3)
dat <- dat %*% c(3600, 60, 1)
dat <- (dat - min(dat)) / (max(dat) - min(dat)) # scale to [0,1]
cm.test(dat, "unif", B = 1000)
# p value is 0.000999001, we reject the null at 1% level.
# So earthquakes don't happen uniformly throughout the day.


