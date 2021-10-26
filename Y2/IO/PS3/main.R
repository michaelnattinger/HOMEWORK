# main code file

setwd("C:/Users/micha/OneDrive/Documents/HOMEWORK/Y2/IO/PS3")

rm(list=ls())
library(readxl)
library(tidyverse)
library(MASS)
library(expm)
library(Matrix)
library(dummies)
library(plyr)
library(FixedPoint)

source("blp.R")
source("blpmerger.R")

# import data
x <- read_excel("cereal_ps3.xls")
d <- read_excel("demog_ps3.xls")

x$t <- paste(x$city, x$quarter, sep="_")
d$t <- paste(d$city, d$quarter, sep="_")
x$brand <- paste(x$firm, x$brand, sep="_")

x$brand <- as.factor(x$brand)
x$t <- as.factor(x$t)
x$id <- as.factor(x$id)
x$firm <- as.factor(x$firm)
x <- merge(x, d, by="t")

x$constant <- rep(1, nrow(x))
n <- 20
theta_nlin=rep(0,20)

###########################################################
# OLS without brand FE (1-3)
blp(theta_nlin=theta_nlin, x=x, brandFE=F, iv=F, supply=F, n=n, A=1)
summary(mean_u)    # regression results
theta_lin=coefficients(mean_u)[1:4]
c(mean(b), quantile(b, 0.5), sd(b))     # markup statistics
c(mean(mc), quantile(mc, 0.5), sd(mc))    # marginal cost statistics

# Post-Nabisco merger
PNmerger <- x
PNmerger$xi <- residuals(mean_u)
PNmerger$mc <- mc
PNmerger[which(PNmerger$firm==3), "firm"] <- 6
PNmerger$firm <- as.factor(as.character(PNmerger$firm))
mer <- blpmerger(theta_nlin=theta_nlin, theta_lin=theta_lin, x=PNmerger, va=va, n=n, brandFE=F)
c(mean(mer[[1]]), quantile(mer[[1]], 0.5), sd(mer[[1]]))
c(mean(mer[[2]]), quantile(mer[[2]], 0.5), sd(mer[[2]]))

#GM-Quaker merger
GQmerger <- x
GQmerger$xi <- residuals(mean_u)
GQmerger$mc <- mc
GQmerger[which(GQmerger$firm==2), "firm"] <- 4
GQmerger$firm <- as.factor(as.character(GQmerger$firm))
mer <- blpmerger(theta_nlin=theta_nlin, theta_lin=theta_lin, x=GQmerger, va=va, n=n, brandFE=F)
c(mean(mer[[1]]), quantile(mer[[1]], 0.5), sd(mer[[1]]))
c(mean(mer[[2]]), quantile(mer[[2]], 0.5), sd(mer[[2]]))
###########################################################

###########################################################
# OLS with brand FE (1-3)
Z <- data.matrix(dummy.data.frame(as.data.frame(x$brand)))
blp(theta_nlin=theta_nlin, x=x, brandFE=T, iv=F, supply=F, n=n, A=(t(Z)%*%Z))
coefs    # regression results (coefficient estimates)
se    # regression results (standard errors)
theta_lin=coefs
c(mean(b), quantile(b, 0.5), sd(b))    # markup statistics
c(mean(mc), quantile(mc, 0.5), sd(mc))    # marginal cost statistics

nam <- names(coefficients(mean_u)[2:length(coefficients(mean_u))])
nam <- gsub("brand", "", nam)
resmd <- data.frame(nam, resmd)
colnames(resmd) <- c("brand", "resmd")

#Post-Nabisco merger
PNmerger <- x
PNmerger$xi <- residuals(mean_u)
PNmerger$mc <- mc
PNmerger[which(PNmerger$firm==3), "firm"] <- 6
PNmerger$firm <- as.factor(as.character(PNmerger$firm))
PNmerger <- join(PNmerger, resmd, by="brand")
mer <- blpmerger(theta_nlin=theta_nlin, theta_lin=theta_lin, x=PNmerger, va=va, n=n, brandFE=T)
c(mean(mer[[1]]), quantile(mer[[1]], 0.5), sd(mer[[1]]))
c(mean(mer[[2]]), quantile(mer[[2]], 0.5), sd(mer[[2]]))

#GM-Quaker merger
GQmerger <- x
GQmerger$xi <- residuals(mean_u)
GQmerger$mc <- mc
GQmerger[which(GQmerger$firm==2), "firm"] <- 4
GQmerger$firm <- as.factor(as.character(GQmerger$firm))
GQmerger <- join(GQmerger, resmd, by="brand")
mer <- blpmerger(theta_nlin=theta_nlin, theta_lin=theta_lin, x=GQmerger, va=va, n=n, brandFE=T)
c(mean(mer[[1]]), quantile(mer[[1]], 0.5), sd(mer[[1]]))
c(mean(mer[[2]]), quantile(mer[[2]], 0.5), sd(mer[[2]]))
###########################################################

###########################################################
# IV without brand FE (1-3)
vars <- seq(1, 20, 1)
vars <- paste("z", vars, sep="")  
Z <- data.matrix(x[,vars])
blp(theta_nlin=theta_nlin, x=x, brandFE=F, iv=T, supply=F, n=n, A=(t(Z)%*%Z))
summary(mean_u)    # regression results
theta_lin=coefficients(mean_u)
c(mean(b), quantile(b, 0.5), sd(b))    # markup statistics
c(mean(mc), quantile(mc, 0.5), sd(mc))    # marginal cost statistics

#Post-Nabisco merger
PNmerger <- x
PNmerger$xi <- residuals(mean_u)
PNmerger$mc <- mc
PNmerger[which(PNmerger$firm==3), "firm"] <- 6
PNmerger$firm <- as.factor(as.character(PNmerger$firm))
mer <- blpmerger(theta_nlin=theta_nlin, theta_lin=theta_lin, x=PNmerger, va=va, n=n, brandFE=F)
c(mean(mer[[1]]), quantile(mer[[1]], 0.5), sd(mer[[1]]))
c(mean(mer[[2]]), quantile(mer[[2]], 0.5), sd(mer[[2]]))

#GM-Quaker merger
GQmerger <- x
GQmerger$xi <- residuals(mean_u)
GQmerger$mc <- mc
GQmerger[which(GQmerger$firm==2), "firm"] <- 4
GQmerger$firm <- as.factor(as.character(GQmerger$firm))
mer <- blpmerger(theta_nlin=theta_nlin, theta_lin=theta_lin, x=GQmerger, va=va, n=n, brandFE=F)
c(mean(mer[[1]]), quantile(mer[[1]], 0.5), sd(mer[[1]]))
c(mean(mer[[2]]), quantile(mer[[2]], 0.5), sd(mer[[2]]))
###########################################################

###########################################################
# IV with brand FE (1-3)
vars <- seq(1, 20, 1)
vars <- paste("z", vars, sep="")  
Z <- data.matrix(cbind(x[,vars], data.matrix(dummy.data.frame(as.data.frame(x$brand)))))
blp(theta_nlin=theta_nlin, x=x, brandFE=T, iv=T, supply=F, n=n, A=(t(Z)%*%Z))
coefs    # regression results (coefficient estimates)
se    # regression results (standard errors)
theta_lin=coefs
c(mean(b), quantile(b, 0.5), sd(b))    # markup statistics
c(mean(mc), quantile(mc, 0.5), sd(mc))    # marginal cost statistics

nam <- names(coefficients(mean_u)[2:length(coefficients(mean_u))])
nam <- gsub("brand", "", nam)
resmd <- data.frame(nam, resmd)
colnames(resmd) <- c("brand", "resmd")

#Post-Nabisco merger
PNmerger <- x
PNmerger$xi <- residuals(mean_u)
PNmerger$mc <- mc
PNmerger[which(PNmerger$firm==3), "firm"] <- 6
PNmerger$firm <- as.factor(as.character(PNmerger$firm))
PNmerger <- join(PNmerger, resmd, by="brand")
mer <- blpmerger(theta_nlin=theta_nlin, theta_lin=theta_lin, x=PNmerger, va=va, n=n, brandFE=T)
c(mean(mer[[1]]), quantile(mer[[1]], 0.5), sd(mer[[1]]))
c(mean(mer[[2]]), quantile(mer[[2]], 0.5), sd(mer[[2]]))

#GM-Quaker merger
GQmerger <- x
GQmerger$xi <- residuals(mean_u)
GQmerger$mc <- mc
GQmerger[which(GQmerger$firm==2), "firm"] <- 4
GQmerger$firm <- as.factor(as.character(GQmerger$firm))
GQmerger <- join(GQmerger, resmd, by="brand")
mer <- blpmerger(theta_nlin=theta_nlin, theta_lin=theta_lin, x=GQmerger, va=va, n=n, brandFE=T)
c(mean(mer[[1]]), quantile(mer[[1]], 0.5), sd(mer[[1]]))
c(mean(mer[[2]]), quantile(mer[[2]], 0.5), sd(mer[[2]]))
###########################################################

###########################################################
# Full model (5-7)
theta_nlin <- c(0, 2, 0, 0, 0.3, 2.2, 0.01, 0.2, 5, 13, -0.2, 1.3, 0, -1, 0, 0, 0.2, 0, 0.3, -0.8)
vars <- seq(1, 20, 1)
vars <- paste("z", vars, sep="")  
Z <- data.matrix(cbind(x[,vars], dummy.data.frame(as.data.frame(x$brand))))
params <- optim(par=theta_nlin, fn=blp, x=x, n=n, brand=T, iv=T, supply=F, A=(t(Z)%*%Z),  
                method="Nelder-Mead", control=list(reltol=0.1, trace=T))
theta_nlin <- params$par
theta_nlin
blp(theta_nlin=theta_nlin, x=x, n=n, brand=T, iv=T, supply=F, A=(t(Z)%*%Z))
coefs
se
c(mean(b), quantile(b, 0.5), sd(b))
c(mean(mc), quantile(mc, 0.5), sd(mc))

nam <- names(coefficients(mean_u)[2:length(coefficients(mean_u))])
nam <- gsub("brand", "", nam)
resmd <- data.frame(nam, resmd)
colnames(resmd) <- c("brand", "resmd")

#Post-Nabisco merger
PNmerger <- x
PNmerger$xi <- residuals(mean_u)
PNmerger$mc <- mc
PNmerger[which(PNmerger$firm==3), "firm"] <- 6
PNmerger$firm <- as.factor(as.character(PNmerger$firm))
PNmerger <- join(PNmerger, resmd, by="brand")
mer <- blpmerger(theta_nlin=theta_nlin, theta_lin=theta_lin, x=PNmerger, va=va, n=n, brandFE=T)
c(mean(mer[[1]]), quantile(mer[[1]], 0.5), sd(mer[[1]]))
c(mean(mer[[2]], na.rm=T), quantile(mer[[2]], 0.5, na.rm=T), sd(mer[[2]], na.rm=T))

#GM-Quaker merger
GQmerger <- x
GQmerger$xi <- residuals(mean_u)
GQmerger$mc <- mc
GQmerger[which(GQmerger$firm==2), "firm"] <- 4
GQmerger$firm <- as.factor(as.character(GQmerger$firm))
GQmerger <- join(GQmerger, resmd, by="brand")
mer <- blpmerger(theta_nlin=theta_nlin, theta_lin=theta_lin, x=GQmerger, va=va, n=n, brandFE=T)
c(mean(mer[[1]]), quantile(mer[[1]], 0.5), sd(mer[[1]]))
c(mean(mer[[2]], na.rm=T), quantile(mer[[2]], 0.5, na.rm=T), sd(mer[[2]], na.rm=T))
