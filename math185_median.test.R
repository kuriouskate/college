median.test = function(x, y, alternative = c("two.sided", "less", "greater"))
{
	if (!is.vector(x)) stop('x needs to be a vector')
	if (!is.vector(y)) stop('y needs to be a vector')
	alternative = match.arg(alternative)

	d = median(c(x,y))
	u = sum(x > d)
	u0 = sum(x == d)
	z = (u + u0/2 - length(x)/2)/sqrt(length(x)/4)
	
	STATISTIC = z
	names(STATISTIC) = 'Z'
	PVAL = switch(alternative, two.sided = 2*(1 - pnorm(abs(z))), less = pnorm(z), greater = 1 - pnorm(z))
 	METHOD = 'Two-sample median test'
	nm_alternative <- switch(alternative, two.sided = "two-sided", less = "true difference in medians is less than 0", greater = "true difference in medians is greater than 0")
	DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
	RVAL = list(statistic = STATISTIC, p.value = PVAL, alternative = nm_alternative, method = METHOD, data.name = DNAME)
    class(RVAL) = "htest"
    return(RVAL)
}
