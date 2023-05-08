################################################################################################
### When using parallele
################################################################################################

rm(list=ls())
setwd("directory_files")
#####################################################################################################
##########################################################################################
################# Fixing the Function which will be used inside the Parallele ############
##########################################################################################

#########################################################################################
##### My function ################
##################################
NRVAT.t1E <-function(n, I, gamm, h.S, a1, b1){
  setwd("directory_files")
  ####################
  ####################
  source("SimKin.R")      ### function used to simulate the response variable
  source("Functionkin.R") ### function used to compute the corrected variance-covariance matrix 
  source("Estimkin.R")    ### function used to estimates the nuissance parameters and the heritability parameter
  source("Get.G.R")       ### function used to get the genotypes obtained form simulate3.exe programm
  source("AFC.R")         ### function used to run the AFC model proposed by Saad et al.
  system(" simulate3.exe") # Here, one must be aware that under windows the simulate3.exe works with system(" simulate3.exe"). When using Linux, one must run  system("wine simulate3.exe")
  H=read.table(file ="directory_files/pedfile.dat") ### function used by simulate3.exe programm
  G=Get.G(H[,-c(1:9)])    ### Getting the genotypes obtained form simulate3.exe programm
  
  ####################################
  ### Prepare file for the kinship ##
  ###################################
  Data=H[,-c(5,6,7,10:49)]
  names(Data)=c("ped", "id", "father", "mother", "sex", "avail")
  ID <- data.frame(FID=Data$ped, IID=Data$id, PAT=Data$father, MAT=Data$mother)
  pedAll <- pedigree(id=Data$id, dadid=Data$father, momid=Data$mother, sex=Data$sex, famid=Data$ped)
  kinship2::kindepth
  kinship2::kinship
  kin2 <- as.matrix(kinship(pedAll))
  
  ############################## 
  #Getting Weight for the SNP 
  p=apply(G, 2, function(x) {sum(x, na.rm = T)/((length(x) - sum(is.na(x))) * 2)})
  
  W=Generate.W(a1,b1,p) # Weigth of the SNP
  ############################## 
  
  X = Generate.X(n)
  mu = Compute.mu(X,gamm)
  Y = Generate.Y(I,n,mu,h.S,kin2)
  gamm.hat = Estim.Gamm(X,Y)
  mu.hat = Compute.mu(X,gamm.hat)
  
  ###############################################
  ## Estimate of the dependance parameter (h^2)##
  ###############################################
  h.S_est=optim(par = h.S, fn = like_general, kin2=kin2, mu.hat=mu.hat,Y=Y, method = "Brent", lower = h.S, upper = 0.1)$par
  
  V.hat = Compute.V(mu.hat,I,n,h.S_est,kin2)
  #V.hat = Compute.V(mu.hat,I,n,h.S,kin2)
  D = Compute.D(mu.hat,X,n)
  A=Compute.A(X,I,mu.hat) 
  B=Compute.B(X,V.hat,I)
  I.A=solve(A)
  M1=V.hat%*%X%*%I.A%*%t(D)
  M2= D%*%solve(A)%*%B%*%I.A%*%t(D)
  
  ### Computation of the Corrected Matrix
  
  MAtt.Cov= V.hat - (1/I)*t(M1)-(1/I)*M1+(1/I)*M2
  
  ### Eigen Decomposition of the Corrected Matrix
  sigmaY = eigen(MAtt.Cov)
  U = sigmaY$vectors
  delta.1 = sigmaY$values
  delta.1[which(delta.1<0)]=0
  sqrt.delta.1 = sqrt(delta.1)
  sqrt.sigmaY = U%*%diag(sqrt.delta.1,sum(n),sum(n))%*%t(U)
  
  ### Computation of the Linear Kernel Functions
  AL=KNL(Z=G, kernel="linear.w",weights=W, n=sum(n),m=ncol(G),rho=ncol(G),gamma=1,d=2)
  #=KN(Z=G, kernel="linear.w",n=sum(n),m=ncol(G), rho=0.7)#m=p=ncol(G)
  KL=sqrt.sigmaY%*%AL%*%sqrt.sigmaY 
  eig=eigen(KL)
  DL=eig$values
  
  ### Computation of the P.Values
  S_obsL=t(Y-mu.hat) %*% AL %*%(Y-mu.hat)
  p.valueL=davies(S_obsL, DL)$Qq
  
  
  ### Computation of the Quadratic Kernel Functions
  AQ=KNL(Z=G, kernel="quadratic.w",weights=W, n=sum(n),m=ncol(G),rho=ncol(G),gamma=1,d=2)
  #KN(Z=G, kernel=methodQ,n=sum(n),m=ncol(G), rho=0.7)#m=p=ncol(G)
  KQ=sqrt.sigmaY%*%AQ%*%sqrt.sigmaY 
  eig=eigen(KQ)
  DQ=eig$values
  
  ### Computation of the P.Values
  S_obsQ=t(Y-mu.hat) %*% AQ %*%(Y-mu.hat)
  p.valueQ=davies(S_obsQ, DQ)$Qq
  
  
  ### Computation of the IBS Kernel Functions
  #AQ=KNL(Z=G, kernel="quadratic",weights=W, n=sum(n),m=ncol(G), rho=h.S)#m=p=ncol(G)
  AIB=KNL(Z=G, kernel="IBS.weighted",weights=W, n=sum(n),m=ncol(G),rho=ncol(G),gamma=1,d=2)#m=p=ncol(G)
  #KNL(Z=G, kernel=methodQ,weights=W, n=sum(n),m=ncol(G), rho=h.S_est)#m=p=ncol(G)
  KIB=sqrt.sigmaY%*%AIB%*%sqrt.sigmaY 
  eig=eigen(KIB)
  DIB=eig$values
  #DQ[which(DQ<0)]=0
  ### Computation of the P.Values
  S_obsIB=t(Y-mu.hat) %*% AIB %*%(Y-mu.hat)
  ###Get the Quadratic Statistic and P.value
  p.valueIB=davies(S_obsIB, DIB)$Qq
  
  
  ### Computation of the Gaussien Kernel Functions
  AG=KNL(Z=G, kernel="gauss",weights=W, n=sum(n),m=ncol(G),rho=ncol(G),gamma=1,d=2)#m=p=ncol(G)
  #KNL(Z=G, kernel=methodG,weights=W, n=sum(n),m=ncol(G), rho=h.S_est)#m=p=ncol(G)
  KG=sqrt.sigmaY%*%AG%*%sqrt.sigmaY 
  eig=eigen(KG)
  DG=eig$values
  #DG[which(DG<0)]=0
  ### Computation of the P.Values
  S_obsG=t(Y-mu.hat) %*% AG %*%(Y-mu.hat)
  ###Get the Gaussien Statistic and P.values
  p.valueG=davies(S_obsG, DG)$Qq
  
  
  ### Computation of the Polynomial Kernel Functions
  AP=KNL(Z=G, kernel="poly",weights=W, n=sum(n),m=ncol(G),rho=ncol(G),gamma=1,d=2)#m=p=ncol(G)
  #KNL(Z=G, kernel=methodG,weights=W, n=sum(n),m=ncol(G), rho=h.S_est)#m=p=ncol(G)
  KP=sqrt.sigmaY%*%AP%*%sqrt.sigmaY 
  eig=eigen(KP)
  DP=eig$values
  #DG[which(DG<0)]=0
  ### Computation of the P.Values
  S_obsP=t(Y-mu.hat) %*% AP %*%(Y-mu.hat)
  ###Get the Gaussien Statistic and P.values
  p.valueP=davies(S_obsP, DP)$Qq
  
  ###########################Using Plink Binary file to get the GDS file ######################
  ###############################
  ##### Using Plink in R ########
  ###############################
  
  ## Set your working directory as the path to the PLINK program files: ##
  
  ##### Getting a genotype in Ped file format ############
  GenoE=H[,-c(5:7)]
  #Geno$V1 <- 1:400
  write.table(GenoE, file = "GenoE.ped", sep = " ",
              row.names = FALSE, col.names = FALSE)
  
  ## Recode the allele in the example file ##
  system("plink --file GenoE --recode --alleleACGT --out GenoE") 
  
  ## Make a binary PED file ##
  # (provide the full path, not just the file name)
  system("plink --file GenoE --make-bed --out GenoE")
  
  # PLINK BED files
  bed.fn <- "directory_files/GenoE.bed"
  fam.fn <- "directory_files/GenoE.fam"
  bim.fn <- "directory_files/GenoE.bim"
  
  
  #################################################################################
  #### getting GDS file from binary ped file using "SeqArray" #####################
  #################################################################################
  
  # remove the temporary gds file
  unlink("genoE.gds", recursive = TRUE, force=TRUE)
  # clean up the fragments, reduce the file size
  #cleanup.gds("genoE.gds.tmp")
  #Delete file if it exists
  #file.remove("genoE.gds")
  
  # convert
  seqBED2GDS(bed.fn, fam.fn, bim.fn, "genoE.gds")
  #seqSummary("genoE.gds")
  
  # open a GDS file
  f <- seqOpen("genoE.gds")
  
  # display the contents of the GDS file in a hierarchical structure
  #f
  
  # look at reference and alternate alleles
  Ref.A=refChar(f)
  Alt.A=altChar(f)
  
  ### Allele
  Alell <- seqGetData(f, "allele")
  
  ##of selected variants:
  geno <- seqGetData(f, "genotype")
  #dim(geno)
  #head(geno)
  
  # take out sample id
  Sple <- seqGetData(f, "sample.id")
  #head(Sple)
  
  # a unique integer ID is assigned to each variant
  variant.id <- seqGetData(f, "variant.id")
  #length(variant.id)
  #head(variant.id)
  
  # Chromosome
  chr <- seqGetData(f, "chromosome")
  #head(chr)
  
  # Position
  pos<- seqGetData(f, "position")
  #head(pos)
  ####################################################
  ### code chunk for getting genotype in GDS format ##
  ####################################################
  geno.file = "directory_files/genoE.gds"
  ###########
  # close the file
  seqClose(f)
  
  ##################################################
  ## Creating The variables in order to use SMMAT ##
  ##################################################
  
  phen=data.frame(id=Sple,Y=Y, X=X) 
  colnames(phen)= c("id","disease","Int","age","sex")
  #head(phen)
  
  ###########################################
  ## Getting the GLMM under Null Hypothesis##
  ###########################################
  rownames(kin2)=colnames(kin2)=Sple
  model0 <- glmmkin(disease ~ age + sex, data = phen, kins = 2*kin2,id = "id", family = binomial(link = "logit"))
  
  
  ####################################################
  ### code chunk for getting group file ##
  ####################################################
  gf=data.frame(rep("Set1",20),chr,pos, Ref.A, Alt.A,rep(1,20))
  names(gf)=NULL
  write.table(gf, file = "SetID.txt", sep="\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
  
  group.file="directory_files/SetID.txt"
  
  #######################################################
  ### code chunk for getting SMMAT ("B", "S", "O", "E")##
  #######################################################
  
  BSOE=SMMAT(model0, group.file = group.file, geno.file = geno.file, 
             MAF.range = c(1e-7, 0.5), miss.cutoff = 1, method = "davies", 
             tests = c("O", "E"))
  
  SMMAT=cbind(BSOE[9],BSOE[11],BSOE[12],BSOE[13],BSOE[16])
  
  
  ############################################################
  ######### Getting AFC (Xcorrec, QLS and SKAT) statistic ####
  ############################################################
  COR = cor(G)
  COR[is.na(COR)] = 0
  Maf = matrix(NA , nr=ncol(G) , nc=2)
  Maf[,1] = colMeans(G[Y==1,])/2
  Maf[,2] = colMeans(G[Y==0,])/2
  WEIGHTS=matrix(diag(W),nc=1,nr=ncol(G))
  #WEIGHTS = matrix(1/(Maf[,2]+1),nc=1,nr=ncol(G))
  CORREC = Xcorrec(MAF = Maf , Pheno = Y, Kin = kin2 , Correlation=COR, Weights = WEIGHTS)
  QLS = Wqls(Genotypes=G, MAF=Maf, Pheno = Y, Kin = kin2 , Correlation=COR, Weights = WEIGHTS)
  SKAT = famSKAT(MAF = Maf , Pheno = Y, Kin = kin2 , Correlation=COR, Weights = WEIGHTS)
  
  ######################################################################################
  ######### Getting GEE_KM asymptoptic p-value and pertubation p-value from gskat ######
  ######################################################################################
  gskat::kindepth
  gskat::kinship
  KM=gskat_seq(Y,XC=X,Z=G,ID)
  KM1=gskat_seq(Y,XC=X,Z=G,ID,pw="Norm")
  
  #### Table of Results######
  
  results=cbind(gamm.hat[1], gamm.hat[2], gamm.hat[3], h.S_est, S_obsL, p.valueL,S_obsQ,p.valueQ, S_obsIB,p.valueIB,S_obsG,p.valueG,S_obsP,p.valueP,SMMAT[1],SMMAT[2],SMMAT[3],SMMAT[4],SMMAT[5],CORREC[3],CORREC[4],QLS[3],QLS[4],SKAT[1],SKAT[2], KM[1], KM[2], KM1[2])
} 

  #####################################################
  ######### Applying The Parallelisation ##############
  #####################################################
  library(parallel)
  library(mvtnorm)
  library(CompQuadForm)
  library(VineCopula)
  library(copula)
  library(SKAT)
  library(Matrix)
  library(GMMAT)
  library(SeqArray)
  library(SeqVarTools)
  library(gskat)
  library(kinship2)
  ##########################################################
  ############## Set the Parameter values  #################
  ##########################################################
  n=rep(c(3,4,8),40) # n (number on individuals) is different from one family to another.
  I=length(n)
  gamm = c(-2,1,1) # Covariates' Coefficients
  h.S = 0          # Here we run with three dirrent values of h^2 = {0, 0.2, 0.5} 
  a1=1             # parameter value for getting the Weigth of the SNP
  b1=25            # parameter value for getting the Weigth of the SNP
  ###########################
  ### Setting the Boostrap ##
  ###########################
  bstrp=10000
  seed = 1100
  write.table(seed, file = "seed.txt", row.names = FALSE, col.names=F)
  #############################################
  ####### Trying with lapply #################
  ############################################
  RNGkind("L'Ecuyer-CMRG") #for reproductible result
  set.seed(2700)
  system.time(result<-lapply(1:bstrp, function(i){NRVAT.t1E(n, I, gamm, h.S, a1, b1)}))
  
###############################################
### Get the Result from Parallele simulation ##
###############################################
 
result=matrix(unlist(result), nrow = bstrp, byrow = T)

##############################################################  
######################################################################################
### Set your working directory as the path to where you want the results file to be ##
######################################################################################  
setwd("directory_files")
write.table(results, file = "Result_h0.txt", sep=" \t", quote = FALSE, row.names = F, col.names = F)
#write.table(result, file = "Result_h02.txt", sep=" \t", quote = FALSE, row.names = F, col.names = F)
#write.table(result, file = "Result_h05.txt", sep=" \t", quote = FALSE, row.names = F, col.names = F)
###############################################
### Read the Results ##
###############################################
OUTPUT=try(read.table("directory_files/Result_h0.txt",header=FALSE, 
             col.names =c("Beta0", "gam1", "gam2", "h.S", "StatL", "p.v_obsbL","StatQ", "p.v_obsbQ","StatIB", "p.v_obsbIB", "StatG", "p.v_obsbG", "StatP", "p.v_obsbP","B.score", "B.pval","S.pval","O.pval","E.pval", "Afc.Xc.Stat","Afc.Xc.pval","Afc.Qls.Stat","Afc.Qls.pval","Afc.Skat.pval.Sat", "Afc.Skat.pval.Davi","gskat.pval.Asym","gskat.pval.PertRad","gskat.pval.PertNorm")))
