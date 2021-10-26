jacob <- function(theta_nlin, theta_lin, x, va, n, Z){
  
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
  
  ddel_dsig_a <- NULL
  ddel_dpi_a <- NULL
  
  for(i in 1:length(xt)){
    
    xtt <- xt[[i]]
    xi <- data.matrix(xtt[,c("constant", "price", "sugar", "mushy")])
    d <- data.matrix(xtt[,c(dj1, dj2, dj3, dj4)])   
    
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
    
    dsij_ddelj <- rowMeans(sij*(1-sij))
    dsj_ddel <- sij%*%t(sij)
    diag(dsj_ddel) <- dsij_ddelj
    
    #partial share/partial sigma
    dsj_dsig <- NULL
    for(j in 1:ncol(xi)){
      xs <- colSums(sweep(sij, 1, t(t(xi[,j])), `*`))
      dsij_dsigk <- (sij*matrix(rep(va[i,], nrow(sij)), ncol=n, byrow = TRUE))*(t(t(xi[,j]))%*%rep(1,n)-matrix(rep(xs, nrow(sij)), ncol=n, byrow = TRUE))
      dsj_dsigk <- rowMeans(dsij_dsigk)
      dsj_dsig <- cbind(dsj_dsig, dsj_dsigk)
    }
    
    #partial share/partial pi
    dsj_dpi <- NULL
    for(j in 1:ncol(xi)){
      dsj_dpij <- NULL
      xs <- colSums(sweep(sij, 1, t(t(xi[,j])), `*`))
      for(k in 1:ncol(xi)){
        nam <- get(paste("dj", k, sep=""))
        dsij_dpik <- (sij*matrix(rep(d[k,nam], nrow(sij)), ncol=n, byrow = TRUE))*(t(t(xi[,k]))%*%rep(1,n)-matrix(rep(xs, nrow(sij)), ncol=n, byrow = TRUE))
        dsj_dpi <- cbind(dsj_dpi, rowMeans(dsij_dpik))
      }      
    }
    
    #partial delta/partial sigma
    ddel_dsig <- solve(dsj_ddel)%*%dsj_dsig
    ddel_dsig_a <- rbind(ddel_dsig_a, ddel_dsig)
    #partial delta/partial pi
    ddel_dpi <- solve(dsj_ddel)%*%dsj_dpi
    ddel_dpi_a <- rbind(ddel_dpi_a, ddel_dpi)
  }
  
  jacob <- cbind(Z, ddel_dsig_a, ddel_dpi_a)
  return(jacob)
}