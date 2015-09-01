#PR1-----------------------------------------------------------------------------------------------------------

collisions<-function(x) {
dat<-as.data.frame(table(x)) #creates data frame where column 1 is unique values of x, and column 2 is their frequencies
C=0
for (i in 1:nrow(dat)) { #for each unique value of x
n<- dat[i,2] #store value's frequency
n<-((n*(n-1))/2) #number of collisions for that value
C<- C+n #adds up total number of collisions
}
return(C)
}

collision.test<- function (x, B=999) {
dat<- as.data.frame(table(x))
val<- dat[,1] #stores unique values of x
R<- 0
C.boot<-0
C.vec<-c()
collisions(x) -> C #stores number of collisions in empirical sample

for (b in 1:B) {
x.boot<- sample(val, length(x), replace=T) #sample uniformly from unique values of x
collisions(x.boot) -> C.boot #stores collisions in sample
C.vec[b]<- C.boot #creates vector of collision numbers
}

R<- sum(C.vec<-C) #sum number of times empirical sample has more collisions
p.value<- R/(B+1) #calculates pvalue
return(p.value)
}
 
 collision.test2<- function (x, B=999) {
 dat<- as.data.frame(table(x))
 val<- dat[,1]
 C=0
 for (i in 1:nrow(dat)) {
 n<-dat[i,2]
 n<- ((n*(n-1))/2)
 C<- C+n
 }
 
 R<- 0
 C.boot<- 0
 C.vec<- c()
 for (b in 1:B) {
 x.boot<- sample(val, length(x), replace=T)
 for (r in 1:length(dat)) {
 n.boot<- dat.boot[r,2]
 n.boot<-((n.boot*(n.boot-1))/2)
 C.boot<- C.boot+n.boot
 }
 C.vec[b]<- C.boot
 }
 R<- sum(C.vec<=C)
 p.value<-R/(B+1)
 return(p.value)
}

#PR2---------------------------------------------------------------------------------------------------------------------

boot.median.test<- function (dat, B=999) {
  med.tot<- median(dat[,1]) #median for all values
  dat<-unstack(dat) #creates columns for different groups
  group.num<-ncol(dat) #number of groups
  
  tab<- matrix(0,2,group,num) #table for storing values above and below median
  for (i in 1:group.num) {
  tab[1,i]<- sum(dat[,i]>med.tot) + ((sum(dat[,i]==med.tot))/2) #first row is # of observations above median
  tab[2,i]<- (nrow(dat)-sum(dat[,i]>med.tot))+((sum(dat[,i]==med.tot))/2) #second row is # of observations below median
  }
  n<- sum(tab) #total # of observations
  p.hat<- rowSums(tab)/n #stores probability of being in each group
  q.hat<- colSums(tab)/n #stores probablity of being above or below median
  
  expected<- n*(p.hat %*% t(q.hat)) #expected frequencies of contingency table
  D<- sum(((tab-expected)^2)/expected) #chi squared test statistic
  D.boot<- c()
  
  for (i in 1:B) {
  x.boot.sample<- sample(1:nrow(tab), n, replace=T, prob=p.hat)
  y.boot.sample<- sample(1:ncol(tab), n, replace=T, prob=q.hat)
  table.boot<- (y.boot.sample %*% x.boot.sample)/n
  D.boot[i]<- sum(((table.boot-expected)/sqrt(expected))^2)
  }
  R<- sum(D.boot<=D)
  p.value<- (R/B)
  return(p.value)
  }
  
  #PR3----------------------------------------------------------------------------------------------------------------
  
  
  
C.vec[b]<- C.boot #creates vector of collision numbers
