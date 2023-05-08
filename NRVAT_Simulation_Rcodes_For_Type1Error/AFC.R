
###### Xcorrec
# MAF: matrix (#Snps * 2): First column contains Minor Allele Frequency (MAF) in cases; Second column contains MAF in controls
# Pheno: matrix (#subjects * 1): this one-column matrix contains 0's amd 1's: 1 for cases and 0 for controls. No missing values are allowed.
# Kin: The kinship matrix (#subjects * #subjects): the subjects must be ordered as the Pheno variable.
# Correlation: Correlation matrix between SNPs (#Snps * #Snps). The user should calculate this matrix beforehand. 
#   Either based on own genotype data (in cases, controls, or both) or based on public databases (e.g., 1000 Genomes Projects, ESP, etc.). 
#   NA values are not allowed. They have to be replaced by zeros.
######

#Correlation = cor(IND[pheno>=0,7:ncol(IND)])
#Correlation[is.na(Correlation)] = 0

Xcorrec = function(MAF, Pheno, Kin, Correlation, Weights)
{
	Na     = length(Pheno[Pheno[,1]==1,])
	Nu     = length(Pheno[Pheno[,1]==0,])
	N      = Na + Nu

	OneN  = matrix(1, nc=1, nr = N)						# The three following lines: prepare the phenotype variables
	Y  = Pheno
	OneHat = matrix( Na/N, nc=1 , nr=N)

	P  = (MAF[,1]*Na + MAF[,2]*Nu)/N					# Estimate MAF in all subjects
	if (is.null(Weights))
	{
		VarSnps = sqrt(P*(1-P))						# Variance of SNPs (2p(1-p))
	} else 
	{
		VarSnps = Weights*sqrt(P*(1-P))					# Variance of SNPs (2p(1-p)) accounting for the prespecified Snp weights
	}
	VarSnps = matrix(VarSnps,nc=1)
	cs      = 2*t(VarSnps) %*% Correlation %*% VarSnps			# This value will account for the correlation between Snps.

	if (is.null(Weights))
	{
	num   = 4*( sum (Na*MAF[,1] - Na*P) )^2					# Numerator of the Xcorrec test statistic
	}else{
	num   = 4*( sum (Na*Weights*MAF[,1] - Na*Weights*P) )^2					# Numerator of the Xcorrec test statistic
	}
	denom = 2*as.numeric(cs)* t(Y - OneHat) %*% Kin %*% (Y - OneHat)	# Denominator of the Xcorrec test statistic
	W   = num/denom								# Xcorrec test statistic
	Pvalue = 1-pchisq(W,1)							# Pvalue from a chi-square proba distribution
	out = t(data.frame(c(num, denom,W,Pvalue)))
	colnames(out)= c("Numerator","Denominator","Xcorrected","Pvalue")
	rownames(out) = "Statistics"
	return(out)
}


###### Wqls
# Genotypes: matrix(#Subjects * #Snps). The genotypes are coded as 0, 1, or 2 copies of the minor allele.
# MAF: matrix (#Snps * 2): First column contains Minor Allele Frequency (MAF) in cases; Second column contains MAF in controls
# Pheno: matrix (#subjects * 1): this one-column matrix contains 0's amd 1's: 1 for cases and 0 for controls. No missing values are allowed.
# Kin: The kinship matrix (#subjects * #subjects): the subjects must be ordered as the Pheno variable.
# Correlation: Correlation matrix between SNPs (#Snps * #Snps). The user should calculate this matrix beforehand. 
#   Either based on own genotype data (in cases, controls, or both) or based on public databases (e.g., 1000 Genomes Projects, ESP, etc.). 
#   NA values are not allowed. They have to be replaced by zeros.
######
Wqls = function(Genotypes, MAF, Pheno, Kin, Correlation, Weights)
{
	Na     = length(Pheno[Pheno[,1]==1,])
	Nu     = length(Pheno[Pheno[,1]==0,])
	N      = Na + Nu

        OneN  = matrix(1, nc=1, nr = N)                                         # The three following lines: prepare the phenotype variables
        Y  = Pheno
        OneHat = matrix( Na/N, nc=1 , nr=N)
	temp = Y - OneHat

	P  = (MAF[,1]*Na + MAF[,2]*Nu)/N		# Estimate MAF in all subjects
	Nsnps   = nrow(MAF)							    # Number of Snps
	if (is.null(Weights))
	{
		VarSnps = sqrt(P*(1-P))						# Variance of SNPs (2p(1-p))
	} else 
	{
		VarSnps = Weights*sqrt(P*(1-P))		# Variance of SNPs (2p(1-p)) accounting for the prespecified Snp weights
	}
	VarSnps = matrix(VarSnps,nc=1)
	cs      = 2*t(VarSnps) %*% Correlation %*% VarSnps			# This value will account for the correlation between Snps.

	if (Nsnps==1)
	{
		S=Genotypes
	}else 
	{ 
		if (is.null(Weights))
		{
		S = apply(Genotypes , 1 , sum)	# Rare variant score: sum of minor alleles accross all Snps.
		}else
		{
		S = Genotypes %*% Weights	# Rare variant score: sum of minor alleles accross all Snps.
		}
	}
	S = matrix(S,nc=1,nrow=length(S))
	Kin = 2*Kin;
	KinInv = solve(Kin) ;
	A = as.numeric(t(Y)%*% KinInv %*% OneN %*% solve(t(OneN) %*% KinInv %*% OneN)) ;
	V = KinInv%*%Y - A * KinInv %*% OneN
	if (is.null(Weights))
	{
	num = (sum((S[Y==1]) * rowSums(KinInv[Y==1,Y==1]))*(1-A) -2*sum(MAF[,2])*(A)*Nu)^2;
	}else
	{
	num = (sum((S[Y==1]) * rowSums(KinInv[Y==1,Y==1]))*(1-A) -2*sum(Weights*MAF[,2])*(A)*Nu)^2;
	}
	denom = as.numeric(cs) * (t(V)%*% Kin %*% V) ;
	W = num/denom ;
	Pvalue = 1-pchisq(W,1) ;
	out = t(data.frame(c(num, denom,W,Pvalue)))
	colnames(out)= c("Numerator","Denominator","WQLS","Pvalue")
	rownames(out) = "Statistics"
	return(out)
}

###### famSKAT
# MAF: matrix (#Snps * 2): First column contains Minor Allele Frequency (MAF) in cases; Second column contains MAF in controls
# Pheno: matrix (#subjects * 1): this one-column matrix contains 0's amd 1's: 1 for cases and 0 for controls. No missing values are allowed.
# Kin: The kinship matrix (#subjects * #subjects): the subjects must be ordered as the Pheno variable.
# Correlation: Correlation matrix between SNPs (#Snps * #Snps). The user should calculate this matrix beforehand. 
#    Either based on own genotype data (in cases, controls, or both) or based on public databases (e.g., 1000 Genomes Projects, ESP, etc.). 
#    NA values are not allowed. They have to be replaced by zeros.
######
famSKAT = function( MAF, Pheno, Kin, Correlation, Weights)
{
	Na     = length(Pheno[Pheno[,1]==1,])
	Nu     = length(Pheno[Pheno[,1]==0,])
	N      = Na + Nu

        OneN  = matrix(1, nc=1, nr = N)                                         # The three following lines: prepare the phenotype variables
        Y  = Pheno
        OneHat = matrix( Na/N, nc=1 , nr=N)

	P  = (MAF[,1]*Na + MAF[,2]*Nu)/N		# Estimate MAF in all subjects
	if (is.null(Weights))
	{
		VarSnps = sqrt(P*(1-P))						# Variance of SNPs (2p(1-p))
	} else 
	{
		VarSnps = Weights*sqrt(P*(1-P))		# Variance of SNPs (2p(1-p)) accounting for the prespecified Snp weights
	}
	VarSnps = matrix(VarSnps,nc=1)
	cz = 2* sum((Y - OneHat) %*% t(Y - OneHat) * Kin )
	Vz = cz * VarSnps %*% t(VarSnps)*Correlation
	if (is.null(Weights))
	{
		Q =  4*Na^2*(sum ( (MAF[,1] - P)^2 ))		# Quadratic form without Weights
	} else
	{
		Q =  4*Na^2*(sum ( Weights*(MAF[,1] - P)^2 ))			# Quadratic form with Weights
	}
	E_Q = sum(diag(Vz))							# Satterwaite approximation (1)
	V_Q = 2* sum( diag (Vz %*% Vz))	# Satterwaite approximation (2)
	Delta = V_Q/(2*E_Q) ;						# Satterwaite approximation (3)
	df= 2*E_Q^2/V_Q;							  # Satterwaite approximation (4)
	Qscaled = Q / Delta ;						# Satterwaite approximation (5)
	Pvalue_Sat = 1-pchisq(Qscaled, df) ;					# Satterwaite approximation (6)

	eig = eigen(Vz, symmetric = T, only.values = T)				  # Davies approximation (1)
	evals = eig$values[eig$values > 1e-06 * eig$values[1]]	# Davies approximation (2)
	Pvalue_Dav =davies(Q, evals, acc = 1e-5)$Qq				      # Davies approximation (3)
	out = t(data.frame(c(Pvalue_Sat,Pvalue_Dav)))
        colnames(out)= c("Satterwaite","Davies")
        rownames(out) = "Pvalue"
        return(out)
}

