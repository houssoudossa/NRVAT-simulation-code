Estim.Gamm = function(X,Y)
{
  X1=X[,2]
  X2=X[,3]
  fit = glm(Y~X1+X2,family=binomial)
  gamm.hat=summary(fit)$coef[,1]
  return(gamm.hat)
}

###################################################
##################### Estime bien h^2 #############
###################################################
like_general<-function(h.S,kin2,mu.hat,Y){
  I=I
  n=n
  
  #X = Generate.X(n)
  #mu = Compute.mu(X,gamm)
  #Y = Generate.Y(I,n,mu,h.S,kin2)
  #gamm.hat = Estim.Gamm(X,Y)
  #mu.hat = Compute.mu(X,gamm.hat)
  
  llikelihood = rep(0,I)
  for(i in 1:I)
  {
    indices.fam=split(1:sum(n), rep(1:I, n))
    Y.fam = Y[indices.fam[[i]]]
    mu.fam = mu.hat[indices.fam[[i]]]
    
    Sigma=h.S*2*kin2[indices.fam[[i]],indices.fam[[i]]]+(1-h.S)*diag(rep(1,n[i]))
    
    for(j in 1:(n[i]-1))
    {
      for(k in (j+1):n[i])
      {
        ### Using MVTNORM ######
        ########################
        #Sigma1=Sigma[c(j,k),c(j,k)]
        
        #copCDF.u1.u2<- pmvnorm(upper=c(qnorm(mu.fam[j]),qnorm(mu.fam[k])),sigma=Sigma[c(j,k),c(j,k)])
        
        ### Using Copula #######
        ########################
        
        cop <- BiCop(family = 1, Sigma[j,k])
        copCDF.u1.u2<- BiCopCDF(mu.fam[j],mu.fam[k],cop)
        
        if ((Y.fam[j]==1)&&(Y.fam[k]==1)) lik.fam = log(copCDF.u1.u2) 
        if ((Y.fam[j]==0)&&(Y.fam[k]==1)) lik.fam = log(mu.fam[k]-copCDF.u1.u2)
        if ((Y.fam[j]==1)&&(Y.fam[k]==0)) lik.fam = log(mu.fam[j]-copCDF.u1.u2)
        if ((Y.fam[j]==0)&&(Y.fam[k]==0)) lik.fam = log(1-mu.fam[j]-mu.fam[k]+copCDF.u1.u2) 
        
        lik.fam=sum(lik.fam)
        
      }
      
    }
    
    llikelihood[[i]] = lik.fam
  }
  #return(sum(llikelihood))
  return(-sum(llikelihood))
}

#h.S_est=optim(h.S,like_general, kin2=kin2, mu.hat=mu.hat,Y=Y)$par