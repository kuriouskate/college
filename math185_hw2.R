
# Pr1------------------------------
chisq.boot.test <- function(tab, B = 2000){
  # Perform chisq GOF test on "tab" and return p value 
  # computed by bootstrap.
  #
  # Args:
  #      tab: a table of counts;
  #      B: number of bootstrap samples we simulate.
  #
  # Returns: p.value  
  
  D <- chisq.test(tab)$stat # compute chisq stat from data
  
  n <- sum(tab)
  p.hat <- rowSums(tab)/n # estimate p.hat and q.hat
  q.hat <- colSums(tab)/n
  table.expected <- n * (p.hat %*% t(q.hat))
  # or  table.expected <- chisq.test(tab)$expected
  
  D.boot <- numeric(B) # initialize D.boot  
  for(i in 1:B){ # simulate boot samples and compute chisq stat
    x.boot.sample <- sample(1 : nrow(tab), n, replace = TRUE, 
                            prob = p.hat)
    y.boot.sample <- sample(1 : ncol(tab), n, replace = TRUE, 
                            prob = q.hat)
    table.boot <- table(x.boot.sample, y.boot.sample)
    
    D.boot[i] <- sum( ( (table.boot - table.expected) / 
                         sqrt(table.expected) )^2 )
  }
  
  p.value <- 1 - sum(D.boot <= D)/(B + 1)
  return(p.value)
}



# Pr2---------------------------
# 1. What percentage of people on the AHA diet had some sort of 
#    illness or death?
(303 - 239) / 303 # = 0.2112211
# 2. What percentage of people on the Mediterranean diet had 
#    some sort of illness or death?
(302 - 273) / 302 # = 0.09602649
# 3. Conduct a Pearson Chi-Square test to determine 
#    if there is any relationship between diet and outcome.
data <- read.table(url('http://www.math.ucsd.edu/~xiz113/math185/w2/w2_hw_data.txt'), header = TRUE)

chisq.test(data[, -1])$p.value
# theoretical p value is 0.0008726344
chisq.test(data[, -1])$expected
# we can rely on it bcz expected counts are greater than 5.

set.seed(185)
chisq.boot.test(data[,-1], 200000)
# bootstrap p value is 0.02248876

# experiment
# chisq.test(data[, -1], sim = TRUE, B = 5000000)$p.value
# [1] 0.0006623999
