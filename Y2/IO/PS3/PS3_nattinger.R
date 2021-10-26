# Implements BLP to solve problem set
# Nattinger 10/23/21

# init libraries
rm(list=ls())
library(MASS)
library(expm)
library(Matrix)
library(dummies)
library(plyr)
library(FixedPoint)
setwd("C:/Users/micha/OneDrive/Documents/HOMEWORK/Y2/IO/PS3")
source("blp.R")
source("blpmerger.R")

# read-in and process data
x <- read.csv("cereal_ps3.csv", header = T, sep = ";")
d <- read.csv("demog_ps3.csv", header = T, sep = ";")
d$t <- paste(d$city, d$quarter, sep="_")

x$t <- paste(x$city, x$quarter, sep="_")
x$brand <- paste(x$firm, x$brand, sep="_")
x$brand <- as.factor(x$brand)
x$t <- as.factor(x$t)
x$id <- as.factor(x$id)
x$firm <- as.factor(x$firm)
x <- merge(x, d, by="t")

x$constant <- rep(1, nrow(x))
n <- 20

#######################
# OLS without brand f.e
theta_nlin <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
A <- 1
blp(theta_nlin=theta_nlin, x=x, brand=F, iv=F, supply=F, n=n, A=A)
summary(mean_u)
c(mean(b), quantile(b, 0.5), sd(b))
c(mean(mc), quantile(mc, 0.5), sd(mc))

#Post-Nabisco merger
merger <- x
merger$xi <- residuals(mean_u)
merger$mc <- mc
merger[which(merger$firm==3), "firm"] <- 6
merger$firm <- as.character(merger$firm)
merger$firm <- as.factor(merger$firm)
theta_lin=coefficients(mean_u)[1:4]
mer <- blpmerger(theta_nlin=theta_nlin, theta_lin=theta_lin, x=merger, va=va, n=n, brand=F)
c(mean(mer[[1]]), quantile(mer[[1]], 0.5), sd(mer[[1]]))
c(mean(mer[[2]]), quantile(mer[[2]], 0.5), sd(mer[[2]]))

#GM-Quaker merger
merger <- x
merger$xi <- residuals(mean_u)
merger$mc <- mc
merger[which(merger$firm==2), "firm"] <- 4
merger$firm <- as.character(merger$firm)
merger$firm <- as.factor(merger$firm)
mer <- blpmerger(theta_nlin=theta_nlin, theta_lin=theta_lin, x=merger, va=va, n=n, brand=F)
c(mean(mer[[1]]), quantile(mer[[1]], 0.5), sd(mer[[1]]))
c(mean(mer[[2]]), quantile(mer[[2]], 0.5), sd(mer[[2]]))

#######################
# OLS with brand f.e
theta_nlin <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
Z <- data.matrix(dummy.data.frame(as.data.frame(x$brand)))
A <- (t(Z)%*%Z)
blp(theta_nlin=theta_nlin, x=x, n=n, brand=T, iv=F, supply=F, A=A)
coefs
se
Rsq
Rsq_G
c(mean(b), quantile(b, 0.5), sd(b))
c(mean(mc), quantile(mc, 0.5), sd(mc))

nam <- names(coefficients(mean_u)[2:length(coefficients(mean_u))])
nam <- gsub("brand", "", nam)
resmd <- data.frame(nam, resmd)
colnames(resmd) <- c("brand", "resmd")

#Post-Nabisco merger
merger <- x
merger$xi <- residuals(mean_u)
merger$mc <- mc
merger[which(merger$firm==3), "firm"] <- 6
merger$firm <- as.character(merger$firm)
merger$firm <- as.factor(merger$firm)
merger <- join(merger, resmd, by="brand")
theta_lin=coefs
mer <- blpmerger(theta_nlin=theta_nlin, theta_lin=theta_lin, x=merger, va=va, n=n, brand=T)
c(mean(mer[[1]]), quantile(mer[[1]], 0.5), sd(mer[[1]]))
c(mean(mer[[2]]), quantile(mer[[2]], 0.5), sd(mer[[2]]))

#GM-Quaker merger
merger <- x
merger$xi <- residuals(mean_u)
merger$mc <- mc
merger[which(merger$firm==2), "firm"] <- 4
merger$firm <- as.character(merger$firm)
merger$firm <- as.factor(merger$firm)
merger <- join(merger, resmd, by="brand")
theta_lin=coefs
mer <- blpmerger(theta_nlin=theta_nlin, theta_lin=theta_lin, x=merger, va=va, n=n, brand=T)
c(mean(mer[[1]]), quantile(mer[[1]], 0.5), sd(mer[[1]]))
c(mean(mer[[2]]), quantile(mer[[2]], 0.5), sd(mer[[2]]))

#######################
# IV without brand f.e
theta_nlin <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
vars <- seq(1, 20, 1)
vars <- paste("z", vars, sep="")  
Z <- data.matrix(x[,vars])
A <- t(Z)%*%Z
blp(theta_nlin=theta_nlin, x=x, n=n, brand=F, iv=T, supply=F, A=A)
summary(mean_u)
c(mean(b), quantile(b, 0.5), sd(b))
c(mean(mc), quantile(mc, 0.5), sd(mc))

#Post-Nabisco merger
merger <- x
merger$xi <- residuals(mean_u)
merger$mc <- mc
merger[which(merger$firm==3), "firm"] <- 6
merger$firm <- as.character(merger$firm)
merger$firm <- as.factor(merger$firm)
theta_lin=coefficients(mean_u)
mer <- blpmerger(theta_nlin=theta_nlin, theta_lin=theta_lin, x=merger, va=va, n=n, brand=F)
c(mean(mer[[1]]), quantile(mer[[1]], 0.5), sd(mer[[1]]))
c(mean(mer[[2]]), quantile(mer[[2]], 0.5), sd(mer[[2]]))

#GM-Quaker merger
merger <- x
merger$xi <- residuals(mean_u)
merger$mc <- mc
merger[which(merger$firm==2), "firm"] <- 4
merger$firm <- as.character(merger$firm)
merger$firm <- as.factor(merger$firm)
theta_lin=coefficients(mean_u)
mer <- blpmerger(theta_nlin=theta_nlin, theta_lin=theta_lin, x=merger, va=va, n=n, brand=F)
c(mean(mer[[1]]), quantile(mer[[1]], 0.5), sd(mer[[1]]))
c(mean(mer[[2]]), quantile(mer[[2]], 0.5), sd(mer[[2]]))

#######################
# IV with brand f.e
theta_nlin <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
vars <- seq(1, 20, 1)
vars <- paste("z", vars, sep="")  
Z <- data.matrix(cbind(x[,vars], data.matrix(dummy.data.frame(as.data.frame(x$brand)))))
A <- t(Z)%*%Z
blp(theta_nlin=theta_nlin, x=x, n=n, brand=T, iv=T, supply=F, A=A)
coefs
se
Rsq
Rsq_G
c(mean(b), quantile(b, 0.5), sd(b))
c(mean(mc), quantile(mc, 0.5), sd(mc))

nam <- names(coefficients(mean_u)[2:length(coefficients(mean_u))])
nam <- gsub("brand", "", nam)
resmd <- data.frame(nam, resmd)
colnames(resmd) <- c("brand", "resmd")

#Post-Nabisco merger
merger <- x
merger$xi <- residuals(mean_u)
merger$mc <- mc
merger[which(merger$firm==3), "firm"] <- 6
merger$firm <- as.character(merger$firm)
merger$firm <- as.factor(merger$firm)
merger <- join(merger, resmd, by="brand")
theta_lin=coefs
mer <- blpmerger(theta_nlin=theta_nlin, theta_lin=theta_lin, x=merger, va=va, n=n, brand=T)
c(mean(mer[[1]]), quantile(mer[[1]], 0.5), sd(mer[[1]]))
c(mean(mer[[2]]), quantile(mer[[2]], 0.5), sd(mer[[2]]))

#GM-Quaker merger
merger <- x
merger$xi <- residuals(mean_u)
merger$mc <- mc
merger[which(merger$firm==2), "firm"] <- 4
merger$firm <- as.character(merger$firm)
merger$firm <- as.factor(merger$firm)
merger <- join(merger, resmd, by="brand")
theta_lin=coefs
mer <- blpmerger(theta_nlin=theta_nlin, theta_lin=theta_lin, x=merger, va=va, n=n, brand=T)
c(mean(mer[[1]]), quantile(mer[[1]], 0.5), sd(mer[[1]]))
c(mean(mer[[2]]), quantile(mer[[2]], 0.5), sd(mer[[2]]))

#######################
# Full model
theta_nlin <- c(0, 2, 0, 0, 0.3, 2.2, 0.01, 0.2, 5, 13, -0.2, 1.3, 0, -1, 0, 0, 0.2, 0, 0.3, -0.8)
vars <- seq(1, 20, 1)
vars <- paste("z", vars, sep="")  
Z <- data.matrix(cbind(x[,vars], dummy.data.frame(as.data.frame(x$brand))))
A <- t(Z)%*%Z
params <- optim(par=theta_nlin, fn=blp, x=x, n=n, brand=T, iv=T, supply=F, A=A,  
                method="Nelder-Mead", control=list(reltol=0.1, trace=T))
theta_nlin <- params$par
#load("theta_nlin_v1.Rdata")
blp(theta_nlin=theta_nlin1, x=x, n=n, brand=T, iv=T, supply=F, A=A)
coefs
se
Rsq
Rsq_G
c(mean(b), quantile(b, 0.5), sd(b))
c(mean(mc), quantile(mc, 0.5), sd(mc))

nam <- names(coefficients(mean_u)[2:length(coefficients(mean_u))])
nam <- gsub("brand", "", nam)
resmd <- data.frame(nam, resmd)
colnames(resmd) <- c("brand", "resmd")

#Post-Nabisco merger
merger <- x
merger$xi <- residuals(mean_u)
merger$mc <- mc
merger[which(merger$firm==3), "firm"] <- 6
merger$firm <- as.character(merger$firm)
merger$firm <- as.factor(merger$firm)
merger <- join(merger, resmd, by="brand")
theta_lin=coefs
mer <- blpmerger(theta_nlin=theta_nlin1, theta_lin=theta_lin, x=merger, va=va, n=n, brand=T)
c(mean(mer[[1]]), quantile(mer[[1]], 0.5), sd(mer[[1]]))
c(mean(mer[[2]], na.rm=T), quantile(mer[[2]], 0.5, na.rm=T), sd(mer[[2]], na.rm=T))

#GM-Quaker merger
merger <- x
merger$xi <- residuals(mean_u)
merger$mc <- mc
merger[which(merger$firm==2), "firm"] <- 4
merger$firm <- as.character(merger$firm)
merger$firm <- as.factor(merger$firm)
merger <- join(merger, resmd, by="brand")
theta_lin=coefs
mer <- blpmerger(theta_nlin=theta_nlin1, theta_lin=theta_lin, x=merger, va=va, n=n, brand=T)
c(mean(mer[[1]]), quantile(mer[[1]], 0.5), sd(mer[[1]]))
c(mean(mer[[2]], na.rm=T), quantile(mer[[2]], 0.5, na.rm=T), sd(mer[[2]], na.rm=T))

#vcov
source("jacob.R")
source("blpvcov.R")
theta_lin=coefs
x <- join(x, resmd, by="brand")
x$xi <- residuals(mean_u)
jacobian <- jacob(theta_nlin=theta_nlin1, theta_lin=theta_lin, x=x, va=va, n=n, Z=Z)
varcov <- blpvcov(theta_nlin=theta_nlin1, theta_lin=theta_lin, x=x, va=va, n=n, 
                  jacob=jacobian, nmkt=94, nprod=24, deg=44, Z=Z)
sqrt(diag(varcov))
