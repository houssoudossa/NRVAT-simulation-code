############################
## Date : 03 - 12 - 2019  ##
## Auteur : Roland        ##
############################

####################################################################
## Generer deux covariables                                       ##
## Une covariable X1: loi uniforme sur [0,1]                      ##
## Une covariable discrete 0/1 avec une probabilite de succes 0.2 ##
####################################################################

Generate.X=function(n)
{
  X1 = runif(sum(n))
  X2 = sample(x=0:1,size=sum(n),replace=TRUE,prob=c(0.8,0.2))
  X  = cbind(rep(1,sum(n)),X1,X2)
  return(X)
}


###############################################################
##### Generate the binary Y using the kinship Matrix ##########  
###############################################################

Generate.Y=function(I,n,mu,h.S,kin2)
{
  Y=NULL
  
  for(i in 1:I){
    indices.fam=split(1:sum(n), rep(1:length(n), n))
    
    mu.fam = mu[indices.fam[[i]]]
    
    ######################################
    #### Compute the kinship matrix ######
    ######################################
    Sigma=h.S*2*as.matrix(kin2[indices.fam[[i]],indices.fam[[i]]])+(1-h.S)*diag(rep(1,n[i]))
    
    ############################
    ### Using MVTNORM ##########
    ############################
    
    #U = pnorm(t(rmvnorm(n=1,sigma=Sigma)))
    
    ############################
    ######## Using Copula ######
    ############################
    
    
    for(j in 1:(n[i]-1))
    {
      for(k in (j+1):n[i])
      {
        
        copspl="Gaussian"
        #rho.cop=Sigma[j,k]
        
        p.cop = normalCopula(param=Sigma[j,k], dim=n[i])
        U= t(rCopula(1,p.cop))
        
        Y.fam = as.numeric(U<=mu.fam)
        Y.fam=cbind(Y.fam)
        
      }
    }
    
    Y = rbind(Y,Y.fam)
  }
  return(Y)
}

###############################################################
##### Generate the binary Y Using the Random Vector ###########  
###############################################################

Generate.Y.b=function(I,n,h.S,X,kin2)   
{
  Y=NULL
  
  for(i in 1:I){
    
    indices.fam=split(1:sum(n), rep(1:length(n), n))
    
    X.fam= X[indices.fam[[i]],]
    
    ##### Using the Random Vector to compute Y #########
    ######################################################
    b=rmvnorm(n=1,sigma=h.S*2*kin2[indices.fam[[i]],indices.fam[[i]]])
    #b=rmvnorm(n=1,sigma=h.S*2*kin2[indices.fam[[i]],indices.fam[[i]]]+(1-h.S)*diag(rep(1,n[i])))
    
    mu.fam=as.vector((exp(X.fam%*%gamm+t(b))/(1+exp(X.fam%*%gamm+t(b)))))
    
    Y.fam = rbinom(n[i],1,mu.fam)
    
    Y.fam=cbind(Y.fam)
    
    Y = rbind(Y,Y.fam)
  }
  return(Y)
}

##############################################
# Function for getting the Weigth of the SNP #
##############################################
Generate.W=function(a1,b1,p)
{
  N=length(p)
  return(diag(dbeta(p,a1,b1),N,N))
}

#####################
## Kernel fonctions##
#####################

K1_Help= function(x,y){
  # Helper function for 2 way interaction kernel
  p = length(x)
  a = x*y
  b = cumsum(a)
  return(sum(a[-1]*b[-p]))
}


call_Kernel_IBS<-function(Z,n,p){
  
  #Kernel_IBS(double * Z, int * pn, int * pp, double * Kernel)
  K<- matrix(rep(0,n*n),nrow = n, ncol = n)	
  temp<-.C("Kernel_IBS",as.integer(as.vector(t(Z))),as.integer(n), as.integer(p),as.double(as.vector(K)))[[4]]
  matrix(temp,nrow=n)
}



call_Kernel_IBS_Weight<-function(Z,n,p,weights){
  
  #Kernel_IBS_Weight(int * Z, int * pn, int * pp, int *UseGivenWeight ,  double * weight, double * Kernel)
  given_weight = 1;
  if( is.null(weights)){
    weights = rep(0,p);
    given_weight = 0;
  } 
  else {
    # change!!
    weights<-weights^2;
  }
  K<- matrix(rep(0,n*n),nrow = n, ncol = n)	
  temp<-.C("Kernel_IBS_Weight",as.integer(as.vector(t(Z))),as.integer(n), as.integer(p),as.integer(given_weight),
           as.double(weights),as.double(as.vector(K)))[[6]]
  matrix(temp,nrow=n)
}


call_Kernel_2wayIX<-function(Z,n,p){
  
  #Kernel_IBS(double * Z, int * pn, int * pp, double * Kernel)
  K<- matrix(rep(0,n*n),nrow = n, ncol = n)	
  temp<-.C("Kernel_2wayIX",as.integer(as.vector(t(Z))),as.integer(n), as.integer(p),as.double(as.vector(K)))[[4]]
  matrix(temp,nrow=n)
}

################################################################################
# GAUSSIAN KERNEL rho > 0 
################################################################################

kernel.gaussian <- function(x, rho)
{
  out <- exp(-1*as.matrix(dist(x) ^ 2)/rho)
  return(out)
}
################################################################################
# POLYNOMIAL KERNEL
################################################################################

kernel.polynomial <- function(x, rho, gamma, d)
{
  return((rho * x %*% t(x) + gamma) ^ d)
}

################################################################################
# Integrating the differents KernelL functions
################################################################################


KN = function(Z, kernel,n,m,rho,gamma,d){
  ## Add linear, linear.W and quadratic.w
  if (kernel == "linear") {
    K = Z%*%t(Z)
  }
  
  
  if (kernel == "gauss") {
    K = exp(-1*as.matrix(dist(Z) ^ 2)/rho)
  }
  
  if (kernel == "poly") {
    K = (rho * Z %*% t(Z) + gamma) ^ d
  }
  
  if (kernel == "quadratic") {
    K = (Z%*%t(Z)+1)**2
  }
  
  
  if (kernel == "IBS") {
    K = call_Kernel_IBS(Z,n,m)
  }
  
  if (kernel == "2wayIX") {
    K = call_Kernel_2wayIX(Z,n,m)
  }  
  
  if (kernel == "IBS_OLD") {
    K1=matrix(nrow=n,ncol=n)
    for (i in 1:n) {
      K1[i,] = apply(abs(t(Z)-Z[i,]),2, sum)
    }
    K = (2*m-K1)/(2*m)
  }
  if (kernel == "2wayIX_OLD") {
    K = 1+Z%*%t(Z)
    N1=  matrix(nrow = n, ncol = n)
    for (i in 1:n){
      for (j in i:n){
        N1[j,i] = N1[i,j] = K1_Help(Z[i,], Z[j,])
      }
    }
    K = K+N1
  }
  return(K)
  
}

KNL = function(Z, kernel, weights,n,m,rho,gamma,d){
  ## Add linear, linear.W and quadratic.w
  if (kernel == "linear") {
    K = Z%*%t(Z)
  }
  
  if (kernel == "linear.w") {
    K = Z%*%weights%*%t(Z)
  }
  
  if (kernel == "gauss") {
    K = exp(-1*as.matrix(dist(Z%*%t(Z))^2)/rho)
  }
  if (kernel == "gauss.w") {
    K = exp(-1*as.matrix(dist(Z%*%weights%*%t(Z))^2)/rho)
  }
  
  if (kernel == "poly") {
    K = (rho * Z %*% t(Z) + gamma)^d
  }
  
  if (kernel == "poly.w") {
    K = (rho * Z%*%weights%*%t(Z) + gamma)^ d
  }
  
  if (kernel == "quadratic") {
    K = (Z%*%t(Z)+1)**2
  }
  
  if (kernel == "quadratic.w") {
    K = (Z%*%weights%*%t(Z)+1)**2
  }
  
  if (kernel == "IBS") {
    K = call_Kernel_IBS(Z,n,m)
  }
  if (kernel == "IBS.weighted") {
    
    K = call_Kernel_IBS_Weight(Z,n,m,weights)
  }
  if (kernel == "2wayIX") {
    K = call_Kernel_2wayIX(Z,n,m)
  }  
  if (kernel == "IBS.weighted_OLD") {
    #K = matrix(nrow = n, ncol = n)
    if (is.null(weights)) {
      qs = apply(Z, 2, mean)/(2)
      weights = 1/sqrt(qs)
    } else {
      weights<-weights^2
    }
    K1 = matrix(nrow =n, ncol = n)
    for (i in 1:n) {
      K1[i,] = apply(abs(t(Z)-Z[i,])*weights,2, sum)
    }
    K= 1-(K1)/(2*sum(weights))
  }
  
  if (kernel == "IBS_OLD") {
    K1=matrix(nrow=n,ncol=n)
    for (i in 1:n) {
      K1[i,] = apply(abs(t(Z)-Z[i,]),2, sum)
    }
    K = (2*m-K1)/(2*m)
  }
  if (kernel == "2wayIX_OLD") {
    K = 1+Z%*%t(Z)
    N1=  matrix(nrow = n, ncol = n)
    for (i in 1:n){
      for (j in i:n){
        N1[j,i] = N1[i,j] = K1_Help(Z[i,], Z[j,])
      }
    }
    K = K+N1
  }
  return(K)
  
}