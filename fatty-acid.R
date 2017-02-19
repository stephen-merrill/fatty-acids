olive2 <- read.table("https://tofu.byu.edu/stat666/datasets/oliver2a",header=T)
olive4 <- read.table("https://tofu.byu.edu/stat666/datasets/oliver4a",header=T)

em.impute <- function(x,tolerance) {
  miss <- is.na(x)
  mu <- matrix(apply(x,2,mean,na.rm=T),nrow=nrow(x),ncol=ncol(x),byrow=T)
  x[miss] <- 0
  x <- x + miss * mu
  sigma <- cov(x)
  tol <- 1
  iter <- 0
  while(tol > tolerance) {
    old.x <- x
    mu <- matrix(apply(x,2,mean,na.rm=T),nrow=nrow(x),ncol=ncol(x),byrow=T)
    sigma <- cov(x)
    for(i in 1:nrow(x)) {
      sig22 <- sigma[!miss[i,],!miss[i,]]
      sig12 <- sigma[miss[i,],!miss[i,]]
      B <- sig12%*%solve(sig22)
      xstar <- mu[1,miss[i,]] + B%*%as.numeric(x[i,!miss[i,]]-mu[1,!miss[i,]])
      x[i,miss[i,]] <- xstar
    }
    iter <- iter + 1
    tol <- max((x-old.x)/old.x)
  }
  x
}

olive2.em <- em.impute(olive2,0.0000001)
olive4.em <- em.impute(olive4,0.0000001)

library(xtable)
xtable(olive2.em)
xtable(olive4.em)

#add error
library(expm)
mult.impute <- function(x,miss) {
  sig <- cov(x)
  for(i in 1:nrow(x)) {
    if(!all(!miss[i,])) {
      sigstarhat <- sig[miss[i,],miss[i,]]-sig[miss[i,],!miss[i,]]%*%solve(sig[!miss[i,],!miss[i,]])%*%sig[!miss[i,],miss[i,]]
      if(sum(miss[i,])==1)
        e <- 1/sqrt(sigstarhat)*rnorm(1)
      else
        e <- solve(sqrtm(sigstarhat))%*%rnorm(sum(miss[i,]))
      nummiss <- 1
      for(j in 1:ncol(x)) {
        if(miss[i,j]) {
          x[i,j] <- x[i,j]+e[nummiss]
          nummiss <- nummiss + 1
        }
      }
    }
  }
  x
}

num.imputes <- 500
p <- ncol(olive2)
betas2 <- matrix(0,nrow=num.imputes,ncol=p)
betas4 <- matrix(0,nrow=num.imputes,ncol=p)
covs2 <- matrix(0,nrow=p,ncol=p)
covs4 <- matrix(0,nrow=p,ncol=p)

for(i in 1:num.imputes) {
  o2 <- mult.impute(olive2.em,is.na(olive2))
  o4 <- mult.impute(olive4.em,is.na(olive4))
  betas2[i,]  <- colMeans(o2)
  betas4[i,]  <- colMeans(o4)
  covs2 <- covs2 + cov(o2)
  covs4 <- covs4 + cov(o4)
}

sum2 <- matrix(0,nrow=p,ncol=p)
sum4 <- matrix(0,nrow=p,ncol=p)
for(i in 1:num.imputes) {
  sum2 <- sum2 + (betas2[i,]-colMeans(betas2))%*%t(betas2[i,]-colMeans(betas2))
  sum4 <- sum4 + (betas4[i,]-colMeans(betas4))%*%t(betas4[i,]-colMeans(betas4))
}

B2 <- (1/num.imputes-1)*sum2
B4 <- (1/num.imputes-1)*sum4

var.xbar2 <- (covs2/num.imputes)+(1+1/num.imputes)*B2
var.xbar4 <- (covs4/num.imputes)+(1+1/num.imputes)*B4

imputed <- TRUE

# T2 for diff in means
n2 <- nrow(olive2.em)
n4 <- nrow(olive4.em)
if(!imputed) {
  Spl <- 1/(n2+n4-2)*((n2-1)*cov(olive2.em)+(n4-1)*cov(olive4.em))
  ybar2 <- colMeans(olive2.em)
  ybar4 <- colMeans(olive4.em)
  
} else {
  Spl <- 1/(n2+n4-2)*((n2-1)*var.xbar2+(n4-1)*var.xbar4)
  ybar2 <- colMeans(betas2)
  ybar4 <- colMeans(betas4)
  
}
T2.regions <- (ybar2-ybar4)%*%solve((1/n2+1/n4)*Spl)%*%t(t((ybar2-ybar4)))
T2.regions
Fstat <- T2.regions*(n2+n4-p-1)/((n2+n4-2)*p)
Fstat
pval <- pf(Fstat,p,n2+n4-p-1,lower.tail=F)
pval

#individual differences
results <- matrix(0,nrow=8,ncol=2)
colnames(results) <- c("t-score","p-value")
rownames(results) <- colnames(olive2)
t1 <- (ybar2[1]-ybar4[1])/sqrt(((n2+n4)/(n2*n4))*Spl[1,1])
results[1,1] <- t1
pval1 <- pt(t1,n2+n4-2,lower.tail = F)
results[1,2] <- pval1
t2 <- (ybar2[2]-ybar4[2])/sqrt(((n2+n4)/(n2*n4))*Spl[2,2])
results[2,1] <- t2
pval2 <- pt(t2,n2+n4-2,lower.tail = F)
results[2,2] <- pval2
t3 <- (ybar2[3]-ybar4[3])/sqrt(((n2+n4)/(n2*n4))*Spl[3,3])
results[3,1] <- t3
pval3 <- pt(t3,n2+n4-2,lower.tail = T)
results[3,2] <- pval3
t4 <- (ybar2[4]-ybar4[4])/sqrt(((n2+n4)/(n2*n4))*Spl[4,4])
results[4,1] <- t4
pval4 <- pt(t4,n2+n4-2,lower.tail = T)
results[4,2] <- pval4
t5 <- (ybar2[5]-ybar4[5])/sqrt(((n2+n4)/(n2*n4))*Spl[5,5])
results[5,1] <- t5
pval5 <- pt(t5,n2+n4-2,lower.tail = T)
results[5,2] <- pval5
t6 <- (ybar2[6]-ybar4[6])/sqrt(((n2+n4)/(n2*n4))*Spl[6,6])
results[6,1] <- t6
pval6 <- pt(t6,n2+n4-2,lower.tail = F)
results[6,2] <- pval6
t7 <- (ybar2[7]-ybar4[7])/sqrt(((n2+n4)/(n2*n4))*Spl[7,7])
results[7,1] <- t7
pval7 <- pt(t7,n2+n4-2,lower.tail = T)
results[7,2] <- pval7
t8 <- (ybar2[8]-ybar4[8])/sqrt(((n2+n4)/(n2*n4))*Spl[8,8])
results[8,1] <- t8
pval8 <- pt(t8,n2+n4-2,lower.tail = T)
results[8,2] <- pval8

#discriminant function for diff in means
a <- solve(Spl)%*%(t(t((ybar2-ybar4))))
D.1o2 <- diag(sqrt(diag(Spl)))
astar <- D.1o2%*%a

# T2 for comp to historical data
hist2 <- c(1300,120,265,7310,820,45,65,28)
hist4 <- c(1230,105,275,7360,830,41,75,38)
if(!imputed) {
  S2 <- cov(olive2.em)
  S4 <- cov(olive4.em)
} else {
  S2 <- var.xbar2
  S4 <- var.xbar4
}

T2.2 <- n2*(ybar2-hist2)%*%solve(S2)%*%t(t((ybar2-hist2)))
T2.2
nu2 <- n2-1
Fstat2 <- T2.2*(nu2-p+1)/(nu2*p)
Fstat2
pval2 <- pf(Fstat2,p,nu2-p+1,lower.tail=F)
pval2

T2.4 <- n4*(ybar4-hist4)%*%solve(S4)%*%t(t((ybar4-hist4)))
T2.4
nu4 <- n4-1
Fstat4 <- T2.4*(nu4-p+1)/(nu4*p)
Fstat4
pval4 <- pf(Fstat4,p,nu4-p+1,lower.tail=F)
pval4

#test for additional info
q <- 2
SplReduced <- Spl[-c(5,7),-c(5,7)]
T2Reduced <- n2*n4*t(ybar2[-c(5,7)]-ybar4[-c(5,7)]) %*% solve(SplReduced) %*% (ybar2[-c(5,7)]-ybar4[-c(5,7)]) / (n2 + n4)
T2add <- ((n2+n4-2)-p)*(T2.regions-T2Reduced)/((n2+n4-2) + T2Reduced)
FstatAdd <- ((n2+n4-2)-p-1+1)*(T2.regions-T2Reduced)/(q*(n2+n4-2+T2Reduced))
pvalReduced <- pf(FstatAdd, q, (n2+n4-2)-p-q+1, lower.tail = F)
