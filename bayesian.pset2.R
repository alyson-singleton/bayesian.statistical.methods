# Chapter 3, Problem 2
library(ggplot2)

alpha.1 <- rbeta(5000, 295, 308)
alpha.2 <- rbeta(5000, 289, 333)

parameter3.2 <- data.frame(alpha.2 - alpha.1)
colnames(parameter3.2) <- c("dif")

ggplot(data=parameter3.2, aes(dif)) +
  geom_histogram(stat = "bin", binwidth = 0.003) + 
  xlab("alpha.2 - alpha.1") + geom_vline(xintercept = 0, linetype="dashed")

mean(parameter3.2>0)

###______________________________________________________________________________________________
###______________________________________________________________________________________________
###______________________________________________________________________________________________

# Chapter 3, Problem 3

mu.c <- 1.013 + (0.24/sqrt(32))*rt(5000,31)
mu.t <- 1.173 + (0.20/sqrt(36))*rt(5000,35)

parameter3.3 <- data.frame(mu.t - mu.c)
colnames(parameter3.3) <- c("dif")

ggplot(data=parameter3.3, aes(dif)) +
  geom_histogram(stat = "bin", binwidth = 0.005) + 
  xlab("mu.t - mu.c") 

quantile(parameter3.3$dif, c(0.025, 0.975))

###______________________________________________________________________________________________
###______________________________________________________________________________________________
###______________________________________________________________________________________________


# Chapter 3, Problem 5
## Part C

# to construct contours of the posterior when we assume EXACT
#function for posterior when measurements are exact
exact.posterior <- function(mu,sd,y){
  ldens <- 0
  for (i in 1:length(y)) 
    ldens <- ldens + log(dnorm(y[i],mu,sd))
  ldens
}

#samples
y <- c(10,10,12,11,9)
#no of samples
n <- length(y)
#sample mean
ybar <- mean(y)
#sample variance
s2 <- sum((y-mean(y))^2)/(n-1)

#grids for potential mu and logsigma values
mu.vec <- seq(3,18,length=200)
logsd.vec <- seq(-2,4,length=200)
#storage for density values
dens.storage <- matrix(NA, 200, 200)
#pull values and store
for (i in 1:200) {
  for (j in 1:200) {
    dens.storage[i,j] <- exact.posterior(mu.vec[i], exp(logsd.vec[j]), y)
  }
}
#build contour plots
dens <- exp(dens.storage - max(dens.storage))
contours <- c(.0001,.001,.01, seq(.05,.95,.05))
contour(mu.vec, logsd.vec, dens, levels=contours, xlab="mu", ylab="log sigma", labex=0, cex=2)

#calculation of the posterior density from the Section 3.2 informed
#   conditional posterior density and the marginal posterior density 
nsim <- 2000
sd <- sqrt((n-1)*s2/rchisq(nsim,4))
mu <- rnorm(nsim,ybar,sd/sqrt(n))

#function to print mean, variance, quantiles of our samples
for.comparison <- function(x){c(mean(x),sqrt(var(x)), quantile(x, c(.025,.25,.5,.75,.975)))}
print(round((rbind(for.comparison(mu),for.comparison(sd))),2))

#function for posterior when measurements are rounded
rounded.posterior <- function(mu,sd,y){
  ldens <- 0
  for (i in 1:length(y)) 
    ldens <- ldens + log(pnorm(y[i]+0.5,mu,sd) - pnorm(y[i]-0.5,mu,sd))
  ldens}

# to construct contours of the posterior when we assume ROUNDED
#storage for density values
dens.storage <- matrix(NA, 200, 200)
#pull values and store
for (i in 1:200) {
  for (j in 1:200) {
    dens.storage[i,j] <- rounded.posterior(mu.vec[i], exp(logsd.vec[j]), y)
  }
}

#build contour plots
dens <- exp(dens.storage - max(dens.storage))
contour (mu.vec, logsd.vec, dens, levels=contours, xlab="mu", ylab="log sigma", labex=0, cex=2)

#calculation of the posterior, sample from the grid approximation to find
#   conditional posterior density and the marginal posterior density 

#number of samples
nsim <- 2000
#sum rowise for conditional prob on sigma
dens.mu <- apply(dens,1,sum)
#sample using probs found and store indices
mu.indices <- sample(1:length(mu.vec), nsim, replace=T, prob=dens.mu)
#pull corresponding mus
mu <- mu.vec[mu.indices]
#build sd vec
sd <- rep(NA,nsim)
#pull sds conditioning on mu
for (i in (1:nsim)) {
  sd[i] <- exp(sample(logsd.vec, 1, prob=dens[mu.indices[i],]))
}
#print mean, variance, quantiles of our samples
print(round((rbind(for.comparison(mu),for.comparison(sd))),2))


# Part D
z <- matrix (NA, nsim, length(y))
for (i in 1:length(y)){
  lower <- pnorm (y[i]-.5, mu, sd)
  upper <- pnorm (y[i]+.5, mu, sd)
  z[,i] <- qnorm (lower + runif(nsim)*(upper-lower), mu, sd)}
mean((z[,1]-z[,2])^2)

###______________________________________________________________________________________________
###______________________________________________________________________________________________
###______________________________________________________________________________________________

# Chapter 3, Problem 8
## Part A
y.samples <- c(16/(16+58), 9/(9+90), 10/(10+48), 
               13/(13+57), 19/(19+103), 20/(20+57), 
               18/(18+86), 17/(17+112), 35/(35+273), 55/(55+64))
z.samples <- c(12/(12+113), 1/(1+18), 2/(2+14), 4/(4+44), 
               9/(9+208), 7/(7+67), 9/(9+29), 8/(8+154))
mean(y.samples)
var(y.samples)

mean(z.samples)
var(z.samples)
mean(y.samples)-mean(z.samples)

#y.dist <- rbeta(10, alpha.y, beta.y)
#z.dist <- rbeta(8, alpha.z, beta.z)

## Part B
alpha.y <- runif(1, 0, 2.709)
beta.y <- runif(1, 0, 10.84)

alpha.z <- runif(1, 0, 2.027)
beta.z <- runif(1, 0, 18.87)

###______________________________________________________________________________________________
###______________________________________________________________________________________________

## Part C

#function for theta_y's joint posterior
y.joint.posterior <- function(alpha.val,beta.val,y){
  ldens <- 1
  for (i in 1:length(y)) 
    ldens <- ldens * ((y[i]^alpha.val)*((1-y[i])^beta.val))
  return(ldens*(dunif(1, 0, 2.709))*(dunif(1, 0, 10.84)))}

# to calculate posterior at different values of alpha and beta
alpha.y.vec <- seq(0,2.709,length=2000)
beta.y.vec <- seq(0,10.84,length=2000)
#storage for density values
dens.storage <- matrix(NA, 2000, 2000)
#pull values and store
for (i in 1:2000) {
  for (j in 1:2000) {
    dens.storage[i,j] <- y.joint.posterior(alpha.y.vec[i], beta.y.vec[j], y.samples)
  }
}

#build contour plots
dens<-dens.storage
contours <- c(10^(-30:0))
contour (alpha.y.vec, beta.y.vec, dens, levels=contours, xlab="alpha", ylab="beta", labex=3, cex=10)

#calculation of the posterior, sample from a grid approximation to find
#   conditional posterior density and the marginal posterior density 
#number of samples

nsim <- 1000
dens.alpha <- apply(dens,1,sum)
alpha.indices <- sample(1:length(alpha.y.vec), nsim, replace=T, prob=dens.alpha)
alpha.y <- alpha.y.vec[alpha.indices]
beta.y <- rep(NA,nsim)
for (i in (1:nsim)) {
  beta.y[i] <- sample(beta.y.vec, 1, prob=dens[alpha.indices[i],])
}
y.finals <- alpha.y/(alpha.y+beta.y)


###______________________________________________________________________________________________

#function for theta_z's joint posterior
z.joint.posterior <- function(alpha.val,beta.val,y){
  ldens <- 1
  for (i in 1:length(y)) 
    ldens <- ldens * ((y[i]^alpha.val)*((1-y[i])^beta.val))
  return(ldens*(dunif(1, 0, 2.027))*(dunif(1, 0, 18.87)))}

# to calculate posterior at different values of alpha and beta
alpha.z.vec <- seq(0,2.027,length=2000)
beta.z.vec <- seq(0,18.87,length=2000)
#storage for density values
dens.storage <- matrix(NA, 2000, 2000)
#pull values and store
for (i in 1:2000) {
  for (j in 1:2000) {
    dens.storage[i,j] <- z.joint.posterior(alpha.z.vec[i], beta.z.vec[j], z.samples)
  }
}
dens.z<-dens.storage
#build contour plots
#dens <- exp(dens.storage - max(dens.storage))
contours <- c(10^(-30:0))
contour (alpha.z.vec, beta.z.vec, dens.z, levels=contours, xlab="alpha", ylab="beta", labex=3, cex=10)

#calculation of the posterior, sample from a grid approximation to find
#   conditional posterior density and the marginal posterior density 
#number of samples
nsim <- 1000
dens.alpha <- apply(dens.z,1,sum)
alpha.indices <- sample(1:length(alpha.z.vec), nsim, replace=T, prob=dens.alpha)
alpha.y <- alpha.z.vec[alpha.indices]
beta.y <- rep(NA,nsim)
for (i in (1:nsim)) {
  beta.y[i] <- sample(beta.z.vec, 1, prob=dens.z[alpha.indices[i],])
}
z.finals <- alpha.y/(alpha.y+beta.y)

## final difference of the means histogram
diff.finals.df <- data.frame(c(y.finals-z.finals))
colnames(diff.finals.df) <- c("z")
mean(diff.finals.df$z)
ggplot(data=diff.finals.df, aes(z)) +
  geom_histogram(stat = "bin", binwidth = 0.02) + 
  xlab("diff.final") + 
  xlab("mu.y - mu.z") + geom_vline(xintercept = 0, linetype="dashed")


###______________________________________________________________________________________________
###______________________________________________________________________________________________
###______________________________________________________________________________________________


# Chapter 3, Problem 12
## Part B
storage <- matrix(NA, 100, 100)
alpha.vec <- seq(0, 50, length=100)
beta.vec <- seq(-30, 20, length=100)

for (i in 1:100) {
  for (j in 1:100) {
    storage[i,j] <- dnorm(alpha.vec[i], 30, 5)*dnorm(beta.vec[j], -3, 2)
  }
}

contour(alpha.vec, beta.vec, storage, xlab="alpha", ylab="beta", xlim=c(20, 40), ylim=c(-10, 10))

## Part E
accidents <- c(24, 25, 31, 31, 22, 21, 26, 20, 16, 22)
time <- c(1:10)

lm(accidents ~ time)
summary(lm(accidents ~ time))

## Part F

storage <- matrix(NA, 100, 100)
#create alpha and beta grids
alpha.vec <- seq(20, 40, length=100)
beta.vec <- seq(-3, 1, length=100)
for (i in 1:100){
  for (j in 1:100){
    sum <- 0
    #calculate sum (log(product)) term
    for (ii in (1:10)) {
      sum <- sum+log(((alpha.vec[i]+beta.vec[j]*time[ii])^(accidents[ii])))
    }
    #calculate and sum with the rest of the expression
    storage[i,j] <- sum + log(exp(-((length(time)*alpha.vec[i])+(beta.vec[j]*sum(time)))))
  }
}
#scale density
dens <- exp(storage - max(storage))
#plot contours
contour(alpha.vec, beta.vec, dens, xlim=c(20, 40), ylim=c(-2.5, 0.5), xlab="alpha", ylab="beta")


## Part G
#calculation of the posterior, sample from a grid approximation to find
#   conditional posterior density and the marginal posterior density 
#number of samples
nsim <- 1000
dens.alpha <- apply(dens,1,sum)
alpha.indices <- sample(1:length(alpha.vec), nsim, replace=T, prob=dens.alpha)
alpha <- alpha.vec[alpha.indices]
beta <- rep(NA,nsim)
for (i in (1:nsim)) {
  beta[i] <- sample(beta.vec, 1, prob=dens[alpha.indices[i],])
}

a.11.b.df <- data.frame(c(alpha+11*beta))
colnames(a.11.b.df) <- c("val")
mean(a.11.b.df$val)
ggplot(data=a.11.b.df, aes(val)) +
  geom_histogram(stat = "bin", binwidth = 0.5) + 
  xlab("alpha + 11(beta)") 

## Part H
num.accidents=rpois(1000, alpha+beta*11)
quantile(num.accidents, c(0.025, 0.975))

## MOM
stor <- numeric(2)
stor[1] <- x[1]/(x[1]+x[2]) - mean(y.values/n.values)
stor[2] <- x[1]*x[2]/(((x[1]+x[2])^2)*(x[1]+x[2]+1)) - sd(y.values/n.values)^2
stor
}
solution <- nleqslv(c(1,1), mom.func) 
mom1 <- solution$x[1]
mom2 <- solution$x[2]

