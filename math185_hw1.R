
# Problem 1-------------------------------------------------------------------------

## Part a-----------------------------
set.seed(185)
chisq.sim.hist <- function(n, k, B){
  # Compare the histogram of chisq statistic with chisq distribution
  #
  # Args:
  #      n: sample size
  #      k: a number indicating the sample space, which is 1,...,k
  #      B: number of repeats
  #
  # Returns: the histogram of D overlaid with density of chisq(k-1)
  
	D <- numeric(B) # initialize a vector to store chisq statistics
	for(b in 1:B){ # run simulation for B times
		x.sim <- sample(1:k, n, replace = TRUE) # simulate data
		D[b]  <- sum((table(x.sim) - n/k)^2 / (n/k)) # compute chisq statistics.
		# alternatively, you can use chisq.test()$stat
	}
	
	hist(D, prob = TRUE, col = "blue", breaks = 2*ceiling(max(D)))
	curve(dchisq(x, k-1), 0, max(D), col = "red", add = TRUE, lwd = 2) # overlay density of chisq(k-1) 
}

## Part b-----------------------------
par(mfrow = c(2, 2))
chisq.sim.hist(10, 5, 2000)
chisq.sim.hist(20, 5, 2000)
chisq.sim.hist(50, 5, 2000)
chisq.sim.hist(100, 5, 2000)

## Part c-----------------------------
chisq.sim.critical <- function(n, k, B){
  # Compute critical value for testing at 5% level
  #
  # Args:
  #      n: sample size
  #      k: a number indicating the sample space, which is 1,...,k
  #      B: number of repeats
  #
  # Returns: critical.value
  
	D <- numeric(B)
	for(b in 1:B){
		x.sim <- sample(1:k, n, replace = TRUE)
		D[b]  <- sum((table(x.sim) - n/k)^2/(n/k))
	}
	critical.value <- sort(D)[floor(0.95*B)]
	return(critical.value)
}

## Part d-----------------------------
set.seed(185)
n <- seq(10, 100, by = 10)
critical.value.sim <- numeric(10) 
for(i in 1:10){
	critical.value.sim[i] <- chisq.sim.critical(n[i], 5, 2000)
}	
par(mfrow = c(1, 1))
plot(n, critical.value.sim, type = "o", pch = 16, col = "blue",
     main = "critical value: simulated v.s. theorerical (k = 5)", 
     xlab = "sample size", ylab = "critical value")
abline(h = qchisq(0.95, 4), lty = 5, col = "red")
legend("bottomright", c("simulated", "theoretical"), 
	   lty = c(1, 5), col = c("blue", "red"))
	
	
	
		

# Problem 2------------------------------------------------------------------------
install.packages("reshape2")
require(reshape2)
dat  <- read.table(url("http://www.math.ucsd.edu/~xiz113/math185/w1/w1_data.txt"),
                   header = T)
data <- acast(dat[, c("Gender", "Month.Code", "Births")], 
              Gender ~ Month.Code)

# Formalize the problem: Let X be a random variable denoting birth 
# month code of female, and define Y similarly for male. Then we test:
# H_0: X and Y have the same distribution
# V.S.
# H_1: X and Y have different distribution

# Summary statistics and plots
data.with.margin <- addmargins(data)
barplot(data, beside = TRUE, legend = TRUE, xlab="MONTHS", args.legend = list(x = "topleft", title="GENDER", inset=c(0,-0.1)))

# Perform chisq test
p.est <- data.with.margin["Sum", - 13]/data.with.margin["Sum", "Sum"] # estimate p
female.expected <- data.with.margin["Female", "Sum"]*p.est # estimated expected counts for female
male.expected   <- data.with.margin["Male", "Sum"]*p.est # estimated expected counts for male
D <- sum( ((data["Female", ] - female.expected)^2/female.expected) + 
          ((data["Male", ] - male.expected)^2/male.expected )) # chisq stat
1 - pchisq(D, 11) # p value = 0.01847684

# Alternativly, we can use chisq.test(data)


# Conclusion:
# We reject the null at 5% significance level and conclude that the birth month of female and male follow the same distribution.  (You can choose other significance level and reach you own conclusion.)
