blp <- function(theta_nlin, x, brand=F, iv=T, supply=T, n, A=NULL){
  
  #Arrange the non linear coefficients
  sig <- theta_nlin[1:4]
  sig <- diag(sig)
  pi <- cbind(theta_nlin[5:8], theta_nlin[9:12], theta_nlin[13:16], theta_nlin[17:20])
  
  #Organize demographic variables
  dj1 <- seq(1, 20, 1)
  dj1 <- paste("v", dj1, sep="")
  dj2 <- seq(21, 40, 1)
  dj2 <- paste("v", dj2, sep="")
  dj3 <- seq(41, 60, 1)
  dj3 <- paste("v", dj3, sep="")
  dj4 <- seq(61, 80, 1)
  dj4 <- paste("v", dj4, sep="")
  
  xt <- split(x, as.factor(x$t))
  
  #Fixed point algorithm to compute mean utilities
  deltas <- NULL
  va <- NULL
  sij_a <- NULL
  
  for(i in 1:length(xt)){
    xtt <- xt[[i]]
    xi <- data.matrix(xtt[,c("constant", "price", "sugar", "mushy")])
    d <- data.matrix(xtt[,c(dj1, dj2, dj3, dj4)])
    sjo <- c(xtt$share, 1-sum(xtt$share)) 
    v <- t(rnorm(n))
    va <- rbind(va,v)
    v <- matrix(rep(t(v),ncol(xi)), ncol=n, byrow = TRUE)
    del <- matrix(1, nrow=nrow(xi), ncol=1)
    sij <- matrix(0, ncol=n, nrow=nrow(xi)+1)
    mu <- matrix(0, ncol=n, nrow=nrow(xi))
    for(j in 1:nrow(d)){
      dj <- rbind(d[j,dj1], d[j,dj2], d[j,dj3], d[j,dj4])
      muj <-  xi%*%sig%*%v + xi%*%pi%*%dj
      mu <- mu+(muj/(nrow(d)))
      u <- del%*%rep(1, n) + muj
      exp_u <- exp(rbind(u, rep(0, ncol(u))))
      sijt <- sweep(exp_u, 2, colSums(exp_u),`/`)
      sij <- sij + sijt
    }
    sj <- rowMeans(sij)
    delp <-  log(sjo) - log(sj) + c(del,1)
    dels <- rbind(c(c(del,1)), c(delp))
    tol <-  dist(dels)
    tol <- as.numeric(tol)
    
    repeat{
      del <- delp
      u <- t(t(del[1:(length(del)-1)]))%*%rep(1, n) + mu
      exp_u <- exp(rbind(u, rep(0, ncol(u))))
      sij <- sweep(exp_u, 2, colSums(exp_u),`/`)
      sj <- rowMeans(sij)
      delp <-  log(sjo) - log(sj) + del
      dels <- rbind(c(del), c(delp))
      tol <-  dist(dels)
      tol <- as.numeric(tol)
      if (tol<1e-14) break
    } 
    deltas <- rbind(deltas, t(t(delp[1:(length(delp)-1)])))
    u <- t(t(delp[1:(length(delp)-1)]))%*%rep(1, n) + mu
    exp_u <- exp(rbind(u, rep(0, ncol(u))))
    sij_a[[i]] <- sweep(exp_u, 2, colSums(exp_u),`/`)
  }
  va <<- va
  
  #Mean utility regression
  if(iv==T){
    vars <- seq(1, 20, 1)
    vars <- paste("z", vars, sep="")  
    vars <- paste(vars, collapse="+")
    vars <- as.name(vars)
    form <- paste("price", vars, sep="~")
    p_h <- predict(lm(as.formula(form), data=x))    
    if(brand==T) vars <- c("-1", "p_h", "brand") else vars <- c("p_h", "sugar", "mushy")
    vars <- paste(vars, collapse="+")
    vars <- as.name(vars)
    form <- paste("deltas", vars, sep="~")
  }else{
    if(brand==T) vars <- c("-1", "price", "brand") else vars <- c("price", "sugar", "mushy")
    vars <- paste(vars, collapse="+")
    vars <- as.name(vars)
    form <- paste("deltas", vars, sep="~")
  }
  mean_u <<- lm(as.formula(form), data=x)
  nam <- names(coefficients(mean_u))
  nam <- which(nam=="price"|nam=="p_h")
  
  #Minimum distance estimates for brand dummy f.e
  if(brand==T){
    ymd <- coefficients(mean_u)[2:length(coefficients(mean_u))]
    hvcov <- vcov(mean_u)[2:length(coefficients(mean_u)),2:length(coefficients(mean_u))]
    ymd <- matrix(c(as.numeric(na.omit(ymd))),nrow=nrow(hvcov), ncol=1) 
    xmd <- xt[[1]]
    xmd <- data.matrix(xmd[,c("constant", "sugar", "mushy")])
    hdmd <- solve(t(xmd)%*%solve(hvcov)%*%xmd)%*%t(xmd)%*%solve(hvcov)%*%matrix(c(ymd),nrow=nrow(hvcov), ncol=1)
    resmd <<- ymd-xmd%*%hdmd
    semd  <-  sqrt(diag(solve(t(xmd)%*%solve(hvcov)%*%xmd)))
    coefs <<- c(hdmd[1], coefficients(mean_u)[nam], hdmd[2:3])
    se <<- c(semd[1], sqrt(vcov(mean_u)[nam,nam]), semd[2:3]) 
    Rsq <<- 1-(t(resmd-mean(resmd))%*%(resmd-mean(resmd)))/(t(ymd-mean(ymd))%*%(ymd-mean(ymd)))
    Rsq_G <<- 1-(t(resmd)%*%solve(hvcov)%*%resmd)/(t(ymd-mean(ymd))%*%solve(hvcov)%*%(ymd-mean(ymd)))
  }
  
  #Marginal cost estimates
  elasticities <- NULL
  mc <- NULL
  b <- NULL
  for(i in 1:length(xt)){
    xtt <- xt[[i]]
    d <- data.matrix(xtt[,c(dj1, dj2, dj3, dj4)])   
    aij <- matrix(0, ncol=n, nrow=nrow(xtt))
    for(j in 1:nrow(d)){
      dj <- rbind(d[j,dj1], d[j,dj2], d[j,dj3], d[j,dj4])     
      a <-  coefficients(mean_u)[nam] + t(t(xtt$price))%*%sig[2,2]%*%va[i,] + t(t(xtt$price))%*%pi[2,]%*%dj
      aij <- aij+a
    }
    aij <- aij/nrow(d)
    sijt <- sij_a[[i]]
    sijt <- sijt[1:(nrow(sijt)-1),]
    sj <- rowMeans(sijt)
    ep <- rowMeans(aij*sijt*(1-sijt)) #own-price derivatives of demand
    ec <- ((-aij)*sijt)%*%t(sijt)/n #cross-price derivatives of demand
    diag(ec) <- ep
    elasticities[[i]] <- ec*(xtt$price/xtt$share) 
    su <- summary(xtt$firm)
    frm <- NULL
    for(k in 1:length(su)) {frm[[k]] <- matrix(1, nrow=su[k], ncol=su[k])}   
    om <-  data.matrix(bdiag(frm))
    om <- (-1)*om*ec #Omega matrix
    bt <- solve(om)%*%sj
    b <- rbind(b, bt)
    mct <- xtt$price-bt
    mc <- rbind(mc, mct)
  }
  mc <<- mc 
  b <<- b 
  
  #Estimation of pricing equation (this is not implemented)
  if(supply==T){
    if(brand==T) w <<- lm(mc ~ sugar + mushy + brand + t, data=x) else w <<- lm(mc ~ sugar + mushy + t, data=x)
    g <- matrix(c(residuals(mean_u),residuals(w)), nrow=length(c(residuals(mean_u),residuals(w))), ncol=1)
    ge <<- g
    if(iv==T){
      vars <- seq(1, 20, 1)
      vars <- paste("z", vars, sep="")  
      if(brand==T) Z <- data.matrix(rbind(cbind(x[,vars], dummy.data.frame(as.data.frame(x$brand))),cbind(x[,vars],dummy.data.frame(as.data.frame(x$brand))))) else Z <- data.matrix(rbind(x[,vars],x[,vars]))      
    }else{
      Z <- data.matrix(rbind(dummy.data.frame(as.data.frame(x$brand)),dummy.data.frame(as.data.frame(x$brand))))    
    }    
  }else{
    g <- residuals(mean_u)
    ge <<- residuals(mean_u)
    if(iv==T){
      vars <- seq(1, 20, 1)
      vars <- paste("z", vars, sep="")  
      if(brand==T) Z <- data.matrix(cbind(x[,vars], dummy.data.frame(as.data.frame(x$brand)))) else Z <- data.matrix(x[,vars])
    }else{
      if(brand==T) Z <- data.matrix(dummy.data.frame(as.data.frame(x$brand)))
    }
  }
  
  #Criterion function
  if(brand == T | iv == T) gmm <-  ((t(g)%*%Z)%*%solve(A)%*%(t(Z)%*%g))/nrow(xt[[1]])
  if(brand == F & iv == F) gmm <-  (t(g)%*%(g))/nrow(xt[[1]])
  if(length(gmm)==0) gmm <- 1e8
  return(gmm)    
}