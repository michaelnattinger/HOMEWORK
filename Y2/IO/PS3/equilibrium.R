equilibrium <- function(p){
  
  xtt$price <- p
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
    muj <-  xi%*%sig%*%matrix(rep(va[k,], ncol(xi)), ncol=n, byrow = TRUE) + xi%*%hpi%*%dj
    u <- del%*%rep(1, n) + muj
    exp_u <- exp(rbind(u, rep(0, ncol(u))))
    sijt <- sweep(exp_u, 2, colSums(exp_u),`/`)/nrow(d)
    sij <- sij + sijt
  }
  sj <- rowMeans(sij)
  
  #Predicted markups
  aij <- matrix(0, ncol=n, nrow=nrow(xtt))
  for(j in 1:nrow(d)){
    dj <- rbind(d[j,dj1], d[j,dj2], d[j,dj3], d[j,dj4])     
    a <-  theta_lin[2] + t(t(xtt$price))%*%sig[2,2]%*%va[k,] + t(t(xtt$price))%*%hpi[2,]%*%dj
    aij <- aij+a
  }
  aij <- aij/nrow(d)
  sijt <- sij[1:(nrow(sij)-1),]
  ep <- rowMeans(aij*sijt*(1-sijt))
  ec <- (-aij*sijt)%*%t(sijt)/n
  diag(ec) <- ep
  su <- summary(xtt$firm)
  frm <- NULL
  for(k in 1:length(su)) {frm[[k]] <- matrix(1, nrow=su[k], ncol=su[k])}   
  om <-  data.matrix(bdiag(frm))
  om <- (-1)*om*ec
  bt <- solve(om)%*%sj[1:(length(sj)-1)]
  hp <- xtt$mc+bt
  return(hp)
}