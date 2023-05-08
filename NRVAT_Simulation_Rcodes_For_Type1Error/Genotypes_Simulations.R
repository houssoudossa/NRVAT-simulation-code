rm(list=ls())
setwd("directory_files")
#####################################################################################################
##########################################################################################
################# Fixing the Function which will be used inside the Parallele ############
##########################################################################################

library(parallel)
library(mvtnorm)
library(CompQuadForm)
library(VineCopula)
library(copula)
library(SKAT)
library(Matrix)
library(kinship2)
################################################
bstrp=10000
############ Seed for getting the same replication of the Simulate3.exe
seed = 1100
write.table(seed, file = "seed.txt", row.names = FALSE, col.names=F)
#######################################################################
Geno = NULL

for (j in 1:bstrp){
  
  setwd("directory_files")
  
  ####################
  source("Get.G.R")
  system(" simulate3.exe") # Here, one must be aware that under windows the simulate3.exe works with system(" simulate3.exe"). When using Linux, one must run  system("wine simulate3.exe")
  H=read.table(file ="directory_files/pedfile.dat") ### function used by simulate3.exe programm
  G=Get.G(H[,-c(1:9)])    ### Getting the genotypes obtained form simulate3.exe programm
  Geno = cbind(Geno,G)
} 

Genotypes = Geno
##########################################################
  