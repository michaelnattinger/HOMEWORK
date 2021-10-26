blpmerger <- function(theta_nlin, theta_lin, x, va, n, brand=F){
  
  #Load function
  source("equilibrium.R")
  brand <<- brand 
  
  #Arrange the non linear coefficients
  sig <- theta_nlin[1:4]
  sig <- diag(sig)
  sig <<- sig
  hpi <<- cbind(theta_nlin[5:8], theta_nlin[9:12], theta_nlin[13:16], theta_nlin[17:20])
    
  #Organize demographic variables
  dj1 <- seq(1, 20, 1)
  dj1 <<- paste("v", dj1, sep="")
  dj2 <- seq(21, 40, 1)
  dj2 <<- paste("v", dj2, sep="")
  dj3 <- seq(41, 60, 1)
  dj3 <<- paste("v", dj3, sep="")
  dj4 <- seq(61, 80, 1)
  dj4 <<- paste("v", dj4, sep="")
  
  xt <- split(x, as.factor(x$t))

  #Post-merger equilibrium
  p_a <- NULL
  sj_a <- NULL
  
  for(i in 1:length(xt)){
    xtt <- xt[[i]]
    k <<- i 
    xtt <<- xtt[order(xtt$firm, xtt$id),]
    mer <- FixedPoint(Inputs=xtt$price, Function=equilibrium) #This line computes the fixed point in prices (price vector that solves p=mc+inv(om(p))*s(p))
    p_a <- rbind(p_a, t(t(mer$FixedPoint)))
    
    xtt$price <- t(t(mer$FixedPoint))
    xi <- data.matrix(xtt[,c("constant", "price", "sugar", "mushy")])
    d <- data.matrix(xtt[,c(dj1, dj2, dj3, dj4)])   
    
    #Predicted shares
    if(brand==F){
      del <- xi%*%theta_lin + xtt$xi
    } else {
      hbr <- xi[,c(1,3,4)]%*%theta_lin[c(1,3,4)] + xtt$resmd
      del <- data.matrix(cbind(xi[,2], dummy.data.frame(as.data.frame(xtt$brand))))%*%c(theta_lin[2], hbr) + xtt$xi
    }
    sij <- matrix(0, ncol=n, nrow=nrow(xi)+1)
    for(j in 1:nrow(d)){
      dj <- rbind(d[j,dj1], d[j,dj2], d[j,dj3], d[j,dj4])
      muj <-  xi%*%sig%*%matrix(rep(va[i,], ncol(xi)), ncol=n, byrow = TRUE) + xi%*%hpi%*%dj
      u <- del%*%rep(1, n) + muj
      exp_u <- exp(rbind(u, rep(0, ncol(u))))
      sijt <- sweep(exp_u, 2, colSums(exp_u),`/`)/nrow(d)
      sij <- sij + sijt
    }
    sj <- rowMeans(sij)
    sj_a <- rbind(sj_a, sj[1:(length(sj)-1)])
  }        
  return(list(p_a, sj_a))
  
  
}