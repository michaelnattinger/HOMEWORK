blpvcov <- function(theta_nlin, theta_lin, x, va, n, jacob, nmkt, nprod, deg, Z){
  
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

  S1 <- matrix(0, 44, 44)
  S2 <- matrix(0, 44, 44)
  S3 <- matrix(0, 44, 44)
  
  for(i in 1:length(xt)){
    
    xtt <- xt[[i]]
    xi <- data.matrix(xtt[,c("constant", "price", "sugar", "mushy")])
    d <- data.matrix(xtt[,c(dj1, dj2, dj3, dj4)])  
    vars <- seq(1, 20, 1)
    vars <- paste("z", vars, sep="")  
    Zt <- data.matrix(cbind(xtt[,vars], dummy.data.frame(as.data.frame(xtt$brand))))
    
    hbr <- xi[,c(1,3,4)]%*%theta_lin[c(1,3,4)] + xtt$resmd
    del <- data.matrix(cbind(xi[,2], dummy.data.frame(as.data.frame(xtt$brand))))%*%c(theta_lin[2], hbr) + xtt$xi
    
    sij <- matrix(0, ncol=n, nrow=nrow(xi)+1)
    for(j in 1:nrow(d)){
      dj <- rbind(d[j,dj1], d[j,dj2], d[j,dj3], d[j,dj4])
      muj <-  xi%*%sig%*%matrix(rep(va[i,], ncol(xi)), ncol=n, byrow = TRUE) + xi%*%pi%*%dj
      u <- del%*%rep(1, n) + muj
      exp_u <- exp(rbind(u, rep(0, ncol(u))))
      sijt <- sweep(exp_u, 2, colSums(exp_u),`/`)/nrow(d)
      sij <- sij + sijt
    }
    sij <- sij[1:(nrow(sij)-1),]
    sj <- rowMeans(sij)
    dsij_dxi <- sij*(1-sij)
    dsj_dxi <- rowMeans(dsij_dxi)
    H <- sij%*%t(sij)
    diag(H) <- dsj_dxi
    
    S1t <- t(Zt)%*%xtt$xi%*%t(xtt$xi)%*%Zt
    
    V2 <- diag(sj)-sj%*%t(sj)
    S2t <- t(Zt)%*%solve(H)%*%V2%*%t(solve(H))%*%Zt
    
    V3 <- (sj-xtt$share)%*%t(sj-xtt$share)
    S3t <- (t(Zt)%*%solve(H)%*%V3%*%t(solve(H))%*%Zt)/n
    
    S1 <- S1+S1t #Intrinsic error
    S2 <- S2+S2t #Error due to sampling
    S3 <- S3+S3t #Error due to simulation
  }
 
  S1 <- as.numeric(((t(x$xi)%*%t(t(x$xi)))/(nprod*nmkt-44)))*S1
  S2 <- (1/(n*nprod*nmkt))*S2
  S3 <- (1/(nprod*nmkt))*S3
  S <- S1+S2+S3
  V <- (ginv(t(jacob)%*%Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%jacob)%*%(t(jacob)%*%Z%*%solve(t(Z)%*%Z)%*%S%*%solve(t(Z)%*%Z)%*%t(Z)%*%jacob)%*%ginv(t(jacob)%*%Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%jacob))
  return(V)
}