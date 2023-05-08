
############################
## Date : 09 - 10 - 2019  ##
## Auteur : Roland        ##
############################

#############################################
## Computation of Mu under Null hypothesis ##
#############################################

Compute.mu = function(X,gamm) {as.vector((exp(X%*%gamm)/(1+exp(X%*%gamm))))}

####################################################
## Computation of Mu under Alternative hypothesis ##
####################################################

Compute.mu.A = function(X,gamm,G,beta) {as.vector((exp(X%*%gamm+ G%*%beta)/(1+exp(X%*%gamm+ G%*%beta))))}

####################################################
## Computation of Sigma with the Heritability h^2 ##
####################################################

compute.Sigma = function(h.S,i) {h.S*rep(1,n[i])%*%t(rep(1,n[i]))+(1-h.S)*diag(rep(1,n[i]))}

################################
## Computation of The matrice ## 
################################

Compute.V = function(mu,I,n,h.S,kin2)
  
{
  V = NULL
  for(i in 1:I){
    indices.fam=split(1:sum(n), rep(1:I, n))
    
    mu.fam = mu[indices.fam[[i]]]
    V.fam = diag(mu.fam*(1-mu.fam))
    
    Sigma=h.S*2*kin2[indices.fam[[i]],indices.fam[[i]]]+(1-h.S)*diag(rep(1,n[i]))
    
    for(j in 1:(n[i]-1))
    {
      for(k in (j+1):n[i])
      {
        ############################
        ### Using MVTNORM ##########
        ############################
        #Sigma.2=Sigma[c(j,k),c(j,k)]
        #V.fam[j,k] = pmvnorm(upper=c(qnorm(mu.fam[j]),qnorm(mu.fam[k])),sigma=Sigma[c(j,k),c(j,k)])- mu.fam[j]*mu.fam[k]
        
        #V.fam[j,k]=pmvnorm(upper=c(qnorm(mu.fam[j]),qnorm(mu.fam[k])),corr=Sigma.2) - mu.fam[j]*mu.fam[k]
        
        ############################
        ######## Using Copula ######
        ############################
        
        cop <- BiCop(family = 1, Sigma[j,k])
        
        copCDF.u1.u2<- BiCopCDF(mu.fam[j],mu.fam[k],cop)
        
        V.fam[j,k]<- copCDF.u1.u2-mu.fam[j]*mu.fam[k]
        
        V.fam[k,j] <- V.fam[j,k]
      }
    }
    V[[i]] = V.fam
  }
  bdiag(V)
}

###########################################
## Computation of A (Fisher Information) ##
###########################################

Compute.A = function(X,I,mu) {(t(X)%*%diag(mu*(1-mu))%*%X)/(I)}

#######################
## Computation of B ###
#######################

Compute.B = function(X,V,I) {(t(X)%*%V%*%X)/(I)}

#########################################
## Computation of D (derivative of Mu) ##
#########################################

Compute.D = function(mu,X,n){
  D = NULL
  for(i in 1:I)
  {
    indices.fam=split(1:sum(n), rep(1:I, n))
    
    mu.fam = mu[indices.fam[[i]]]
    X.fam = X[indices.fam[[i]],]
    result = diag(mu.fam*(1-mu.fam))%*%X.fam
    D[[i]] = result
  }
  #As D is liste, we transforme it as matrice
  D=do.call(rbind,D)
  D
  #bdiag(D)
}
