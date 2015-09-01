
# Pr1----------------------------------------
perm.var.test <- function(x, y, B = 2000){
	# goodness-of-fit perm test based on the ratio of the sample variances
	# 
	# Args: x, y: sample; B: repeats
	#
	# Return: p value
		
	m <- length(x)
	n <- length(y)
	D <- var(x)/var(y)
	
	z <- c(x, y)
	D.perm <- numeric(B)
	for(b in 1:B){
		ind <- sample(m + n)
		D.perm[b] <- var(z[ind[1 : m]])/var(z[ind[(m + 1) : (m + n)]])
	}
	p1 <- (sum(D.perm >= D) +1)/(B+1)
	p2 <- (sum(D.perm <= D) + 1)/(B+1)
	p.value <- 2 * min(p1, p2) # p-value for the two-sided test 
	return(p.value)
}

# Alternatively, one can use the test statistic D = abs(log(var(x)/var(y))) and reject for large values, so that the p-value in this case is simply (sum(D.perm >= D) +1)/(B+1)


# Example:
x = rnorm(500)
y = rnorm(300)
perm.var.test(x, y)
perm.var.test(x, 2*y)




# Pr2----------------------------------------
var0.test <- function(x, y){
	# Perform var0.test
	# Args: x, y: two samples
	# Returns: p value 
	
	D <- (sum(y^2)/length(y))/(sum(x^2)/length(x))
	p.value <- 1 - pf(D, length(y), length(x))
	return(p.value)
}

# Example:
x = rnorm(500)
y = rnorm(300)
var0.test(x, y)
var0.test(x, 2*y)


# Pr3--------------------------------------
# note that the application of the log transform is optional as it does not change anything
Pr3 <- function(theta = 1.1, repeats = 200, 
			    sample.size = c(5, 10, 15, 20, 50, 100)){
  # Compare rank sum test (w/ log) with var0.test
  #
  # Args: repeats
  #       sample.size: a vector of sample sizes
  #       
  # Returns: a dataframe with two columns, recording median p values 
  #          of two methods at each sample size.
  
  result.1 <- result.2 <- numeric(length(sample.size))
  # vector for storing median p value of all repeats
  for(m in sample.size){
    pval.1 <- pval.2 <- numeric(repeats) 
    # vector for storing p value of all repeats
    for(i in 1:repeats){
      x <- rnorm(m)
      y <- rnorm(m, sd = theta)
      pval.1[i] <- wilcox.test(log(abs(x)), log(abs(y)), 
                               alternative = "less", 
                               paired = FALSE)$p.value
      pval.2[i] <- var0.test(x, y)
    }
    result.1[match(m, sample.size)] <- median(pval.1)
    result.2[match(m, sample.size)] <- median(pval.2)
  }
  return(data.frame("log.rank.sum" = result.1, "var0.test" = result.2, 
                    row.names = sample.size))
}

set.seed(185)
Pr3()

set.seed(185)
Pr3(theta = 2) # try other theta value



# Pr4-----------------------------------------
boot.oneway.test <- function(dat, B = 200){
	# Calibrate one way F-test by bootstrap
	#
	# Args: dat: dataframe with two columns, observations in the first
	#            and group memberships in the second
	#       B: boot repeats
	# 
	# Returns: p value
	
    colnames(dat) = c("values", "ind")
    f <- oneway.test(values ~ ind, data = dat, var.equal = T)$statistic
    
    dat <- stack(by(dat[, "values"], dat[, "ind"], function(x) x - mean(x))) 
	f.boot <- numeric(B)			
	for(b in 1:B){
		dat.boot <- stack(by(dat[, "values"], dat[, "ind"], sample, replace = TRUE))
		f.boot[b] <- oneway.test(values ~ ind, data = dat.boot, 
                             var.equal = T)$statistic
	}	
	p.value <- 1 - sum(f.boot <= f)/(B + 1)		
	return(p.value)	 
}

boot.oneway.test(stack(smokers))



# Pr5----------------------------------------

perm.r.measure <- function(mat, B = 200){
  # Testing whether different measures are exchangeable
  #
  # Args: 
  # 	mat: a matrix with blocks indexing rows and measures
  #            indexing columns
  #	B: repeats
  #
  # Returns: p value
  
  if(!is.matrix(mat)){
    mat <- as.matrix(mat)
  }
  
  y.. <- mean(mat)
  g <- nrow(mat) * sum((apply(mat, 2, function(x) mean(x) - y..))^2)
  
  g.perm <- numeric(B)
  for(b in 1:B){
    mat.perm  <- t(apply(mat, 1, sample, replace = FALSE))
    g.perm[b] <- nrow(mat.perm) * sum((apply(mat.perm, 2, function(x) mean(x) - y..))^2)
  }
  p.value <- sum(g.perm >= g)/B
  return(list(p.value = p.value, g = g, g.perm = summary(g.perm)))
}

dat <- read.table(url("http://www.statsci.org/data/general/energy.txt"),
                  header = T)
perm.r.measure(dat[, -1]) # exclude "Subject" column


