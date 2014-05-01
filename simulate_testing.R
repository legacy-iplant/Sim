## Get and transform data
qtrait <- read.table(file="test_qtrait.txt", header=FALSE)
genome <- read.csv(file="test_genomes.csv", header=TRUE)
genome <- genome[,-c(1,2)]

actual <- c(rep(2,3),rep(0,27))

gen <- matrix(ncol=ncol(genome)/2,nrow=nrow(genome))
get <- function(x) {
  num <- x+(x-1)
  return(c(num,num+1))
}
for (i in 1:ncol(gen)) {
  for (j in 1:nrow(gen)) {
    gen[j,i] <- sum(genome[j,get(i)])
  }
}
total <- cbind(gen,qtrait)
names(total) <- c(names(total[,-c(ncol(total))]),"qtrait")

## Basic Linear Model
fit1 <- lm(qtrait~., data=total)
summary(fit1)
error1 <- mean(abs(coef(fit1)[2:31]-actual))

## Lasso Model
require(genlasso)
X = as.matrix(total[,-c(ncol(total))])
D = diag(1,ncol(X))
fit2 <- genlasso(total$qtrait,X,D)
error2 <- mean(abs(coef(fit2)$beta[,ncol(coef(fit2)$beta)]-actual))