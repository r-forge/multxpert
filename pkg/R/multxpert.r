# PDunnett is a secondary function which computes the cumulative 
# distribution function of the Dunnett distribution in one-sided 
# hypothesis testing problems with a balanced one-way layout and 
# equally weighted null hypotheses
#library(mvtnorm)
pdunnett<-function(x,df,m)
# X, Argument
# DF, Number of degrees of freedom
# M, Number of comparisons
{
	# Correlation matrix
	corr<-matrix(0.5,m,m)
	for (i in 1:m) corr[i,i]<-1
	p<-pmvt(lower=rep(-Inf,m), upper=rep(x,m), delta=rep(0,m), df=df, corr=corr, algorithm=GenzBretz(maxpts=25000, abseps=0.00001, releps=0))[1]    
	return(p)
}
# End of pdunnett

# QDunnett is a secondary function which computes a quantile of the 
# Dunnett distribution in one-sided hypothesis testing problems 
# with a balanced one-way layout and equally weighted null hypotheses
#library(mvtnorm)
qdunnett<-function(x,df,m)
# X, Argument
# DF, Number of degrees of freedom
# M, Number of comparisons
{
	# Correlation matrix
	corr<-matrix(0.5,m,m)
	for (i in 1:m) corr[i,i]<-1
	temp<-qmvt(x,interval=c(0,4),tail="lower.tail",df=df, delta=rep(0,m),corr=corr, algorithm=GenzBretz(maxpts=25000, abseps=0.00001, releps=0))[1]
	return(temp$quantile)
}
# End of qdunnett

# LevelProc is a secondary function which computes the significance level for 
# the next family for multistage parallel gatekeeping procedures in hypothesis 
# testing problems with equally weighted null hypotheses
levelproc<-function(reject,proc,gamma,level)
# REJECT, Vector of rejection decisions (1, if rejected; 0, if accepted)
# PROC, Procedure name
# GAMMA, Truncation parameter
# LEVEL, Significance level in the current family
{
	# Total number of null hypotheses in the current family
	m<-length(reject)
	# Number of null hypotheses accepted in the current family
	a<-m-sum(reject)
	# Ratio of accepted hypotheses in the current family
	ratio<-a/m
	
	if (ratio==0) error<-0
	if (ratio>0)
	{
		if (proc=="Holm" | proc=="Hommel" | proc=="Hochberg")
		{
			error<-gamma+(1-gamma)*ratio
		}
		else if (proc=="Fallback")
		{
			accept<-1-reject
			sum<-0
			smallest<-1
			for (i in 1:length(accept))
			{
				# Loop over accepted null hypotheses
				if (accept[i]==1)
				{
					# i is not the smallest index
					if (smallest==0)
					{
						largest<-1
						for (j in 1:(i-1))
						{
							if (accept[j]==1) largest<-j
						}
						l<-largest
					}
					
					# i is the smallest index
					if (smallest==1)
					{
						l<-0
						smallest<-0
					}
					sum<-sum+(i-l)
				}
				error<-gamma*(sum/m)+(1-gamma)*ratio
			}
		}
	}
	
	# Significance level in the next family
	nextlevel<-level-error*level
	return(nextlevel)
}
# End of levelproc

# ParGateEval is a secondary function which evaluates decision rules for 
# multistage parallel gatekeeping procedures in hypothesis testing problems 
# with equally weighted null hypotheses
pargateeval<-function(gateproc,alpha,independence)
# GATEPROC, List of gatekeeping procedure parameters
# ALPHA, Global familywise error rate
# INDEPENDENCE, Boolean indicator (TRUE, Independence condition is imposed; FALSE, 
# Independence condition is not imposed)
{
	
	# Number of families
	nfams<-length(gateproc)
	
	# Significance level in the first family
	level<-alpha
	
	# Test families from first to last
	for(i in 1:nfams)
	{
		# Raw p-values in the current family
		p<-gateproc[[i]]$rawp
		# Null hypotheses are equally weighted in the current family
		w<-rep(1/length(p),length(p))
		# Procedure and truncation parameter in the current family
		if (gateproc[[i]]$proc=="Bonferroni") 
		{
			proc<-"Holm"
			gamma<-0
			gateproc[[i]]$proc<-"Holm"
			gateproc[[i]]$procpar<-0
		}
		else
		{
			proc<-gateproc[[i]]$proc
			gamma<-gateproc[[i]]$procpar
		}
		# Compute the local adjusted p-values in the current family
		adjp<-pvaltrunc(p,w,proc,gamma)
		# Rejection decisions in the current family (0 if rejected and 1 if accepted)
		rejection<-(adjp<=level)
		# Save the local adjusted p-values in the current family
		gateproc[[i]][5]<-list(adjp=adjp)
		# Save the significance level in the current family
		gateproc[[i]][6]<-list(level=level)
		# Save the rejection decisions in the current family
		gateproc[[i]][7]<-list(rejection=rejection)
		# Compute the significance level in the next family
		level<-levelproc(rejection,proc,gamma,level)
	}
	
	# Retest families from last to first if the independence condition is not imposed    
	if (independence==FALSE)    
	{        
		for(i in (nfams+1):(2*nfams-1))
		{
			# Family index
			k<-2*nfams-i     
			prevrej<-gateproc[[i-1]][[7]]   
			# Full alpha level if all null hypotheses are rejected in the previous family            
			if (sum(prevrej)==length(prevrej)) level<-alpha else level<-0            
			# Label in the current family
			label<-gateproc[[k]]$label
			# Raw p-values in the current family
			p<-gateproc[[k]]$rawp
			# Null hypotheses are equally weighted in the current family
			w<-rep(1/length(p),length(p))
			# Procedure in the current family
			if (gateproc[[k]]$proc=="Bonferroni") proc<-"Holm" else proc<-gateproc[[k]]$proc
			# Truncation parameter in the current family (regular procedure is used)
			gamma<-1
			# Compute the local adjusted p-values in the current family
			adjp<-pvaltrunc(p,w,proc,gamma)
			# Rejection decisions in the current family (0 if rejected and 1 if accepted)
			rejection<-(adjp<=level)
			# Create a new list for the current family
			gateproc[[i]]<-list(label=label, rawp=p, proc=proc, procpar=gamma, adjp=adjp, level=level, rejection=rejection)            
		}
	}
	
	return(gateproc)
	
}
# End of pargateeval# Bonf function is a secondary function which computes the weighted Bonferroni
# p-value for an intersection hypothesis
bonf<-function(p,w)
# P, Vector of raw p-values
# W, Vector of hypothesis weights
{
	# Number of null hypotheses included in the intersection
	k<-length(w[w!=0])
	if (k>0) bonf<-min(p[w!=0]/w[w!=0])
	if (k==0) bonf<-1
	return(bonf)
}
# End of bonf




# Simes function is a secondary function which computes a truncated version of the
# weighted Simes p-value for an intersection hypothesis
simes<-function(p,w,gamma)
# P, Vector of raw p-values
# W, Vector of hypothesis weights
# GAMMA, Truncation parameter
{
	# Number of null hypotheses included in the intersection
	k<-length(w[w!=0])
	if (k>1)
	{
		temp<-matrix(0,3,k)# TODO: Add comment
		# P-values
		temp[1,]<-p[w!=0]
		# Weights
		temp[2,]<-w[w!=0]
		# Normalized weights
		temp[3,]<-w[w!=0]/sum(w[w!=0])
		# Sort by p-values
		sorted<-temp[,order(temp[1,])]
		numer<-sorted[1,]
		denom<-gamma*cumsum(sorted[3,])+(1-gamma)*sorted[2,]
		simes<-min(numer/denom)
	}
	if (k==1)
	{
		numer<-p[w!=0]
		denom<-gamma+(1-gamma)*w[w!=0]
		simes<-numer/denom
	}
	if (k==0) simes<-1
	return(simes)
}
# End of simes

# IncSimes function is a secondary function which computes a truncated version of the weighted
# incomplete Simes p-value for an intersection hypothesis
incsimes<-function(p,w,gamma)
# P, Vector of raw p-values
# W, Vector of hypothesis weights
# GAMMA, Truncation parameter
{
	# Number of null hypotheses included in the intersection
	k<-length(w[w!=0])
	if (k>1)
	{
		temp<-matrix(0,3,k)
		# P-values
		temp[1,]<-p[w!=0]
		# Weights
		temp[2,]<-w[w!=0]
		# Normalized weights
		temp[3,]<-w[w!=0]/sum(w[w!=0])
		# Sort by p-values
		sorted<-temp[,order(temp[1,])]
		modw<-w[w!=0]
		modw[1]<-0
		modw[2:k]<-sorted[3,1:k-1]
		numer<-sorted[1,]
		denom<-gamma*sorted[3,]/(1-cumsum(modw))+(1-gamma)*sorted[2,]
		incsimes<-min(numer/denom)
	}
	if (k==1)
	{
		numer<-p[w!=0]
		denom<-gamma+(1-gamma)*w[w!=0]
		incsimes<-numer/denom
	}
	if (k==0) incsimes<-1
	return(incsimes)
}
# End of incsimes

# PValTrunc function is a secondary function which computes adjusted p-values for
# truncated versions of the Holm, Hommel, Hochberg, fixed-sequence and fallback procedures
# in general hypothesis testing problems (equally or unequally weighted null hypotheses)
pvaltrunc<-function(rawp,weight,proc,gamma)
# RAWP, Vector of raw p-values
# WEIGHT, Vector of hypothesis weights
# PROC, Procedure name
# GAMMA, Truncation parameter
{
	# Number of null hypotheses
	nhyps<-length(rawp)

	# Number of weights
	nweis<-length(weight)

	# Number of intersection hypotheses
	nints<-2**nhyps-1

	# Decision matrix
	h<-matrix(0,nints,nhyps)
	for(i in 1:nhyps)
	{
		for(j in 0:(nints-1))
		{
			k<-floor(j/2**(nhyps-i))
			if (k/2==floor(k/2)) h[j+1,i]<-1
		}
	}

	# Temporary vector of weights
	tempw<-matrix(0,1,nhyps)

	# Holm procedure
	if (proc=="Holm")
	{
		# Local p-values for Holm procedure
		holmp<-matrix(0,nints,nhyps)

		for(i in 1:nints)
		{
			tempw<-(weight*h[i,]/sum(weight*h[i,]))*gamma+(1-gamma)*weight*h[i,]
			holmp[i,]<-h[i,]*bonf(rawp,tempw)
		}

		# Compute adjusted p-values
		holm<-rep(0,nhyps)
		for(i in 1:nhyps) holm[i]<-pmin(1, max(holmp[,i]))

		# Delete temporary objects
		rm(holmp)

		# List of adjusted p-values
		res<-holm
	}

	# Hommel procedure
	else if (proc=="Hommel")
	{
		# Local p-values for Hommel procedure
		hommp<-matrix(0,nints,nhyps)

		for(i in 1:nints) hommp[i,]<-h[i,]*simes(rawp,weight*h[i,],gamma)

		# Compute adjusted p-values
		hommel<-rep(0,nhyps)
		for(i in 1:nhyps) hommel[i]<-max(hommp[,i])

		# Delete temporary objects
		rm(hommp)

		# List of adjusted p-values
		res<-hommel
	}

	# Hochberg procedure
	else if (proc=="Hochberg")
	{
		# Local p-values for Hochberg procedure
		hochp<-matrix(0,nints,nhyps)

		for(i in 1:nints) hochp[i,]<-h[i,]*incsimes(rawp,weight*h[i,],gamma)

		# Compute adjusted p-values
		hochberg<-rep(0,nhyps)
		for(i in 1:nhyps) hochberg[i]<-max(hochp[,i])

		# Delete temporary objects
		rm(hochp)

		# List of adjusted p-values
		res<-hochberg
	}

	# Fallback procedure
	else if (proc=="Fallback")
	{
		# Local p-values for fallback procedure
		fallp<-matrix(0,nints,nhyps)
		fallwin<-matrix(0,nhyps,nhyps)

		# Compute window variables
		for(i in 1:nhyps)
		{
			for(j in 1:nhyps)
			{
				if (i>=j) fallwin[i,j]<-1
			}
		}

		for(i in 1:nints)
		{
			tempw[]<-0
			for(j in 1:nhyps)
			{
				if (h[i,j]==1) tempw[j]<-(sum(weight*fallwin[j,])-sum(tempw))
				if (h[i,j]==0) tempw[j]<-0
			}
			tempw<-tempw*gamma+(1-gamma)*weight*h[i,]
			fallp[i,]<-h[i,]*bonf(rawp,tempw)
		}

		# Compute adjusted p-values
		fallback<-rep(0,nhyps)
		for(i in 1:nhyps) fallback[i]<-pmin(1,max(fallp[,i]))

		# Delete temporary objects
		rm(fallp)
		rm(fallwin)

		# List of adjusted p-values
		res<-fallback
	}

	# Fixed-sequence procedure
	else if (proc=="Fixed-sequence")
	{
		# Compute adjusted p-values using a recursive algorithm
		fixedseq<-rep(0,nhyps)
		fixedseq[1]<-rawp[1]
		if (nhyps>1)
		{
			for(i in 2:nhyps) fixedseq[i]<-max(fixedseq[i-1],rawp[i])
		}

		# List of adjusted p-values
		res<-fixedseq
	}

	rm(tempw)
	return(res)
}
# End of pvaltrunc
# ParAdjP function computes adjusted p-values for the single-step Dunnett
# procedure and step-down Dunnett procedure in one-sided hypothesis testing
# problems with a balanced one-way layout and equally weighted null hypotheses
paradjp<-function(stat,n,proc=c("Single-step Dunnett", "Step-down Dunnett"))
# STAT, Vector of test statistics
# N, Common sample size in each treatment group
# PROC, Procedure name
{

	if (n<=0) stop("Sample size must be positive")

	# Number of null hypotheses
	m<-length(stat)

        if(m==0) stop("No test statistics are specified")

        if(length(n)==0) stop("No sample size specified")

	# Number of degrees of freedom
	nu<-(m+1)*(n-1)
	# Raw p-values
	rawp<-1-pt(stat,2*(n-1))
	# Adjusted p-values
	adjp<-rep(0,m)

        if(!all(proc %in% c("Single-step Dunnett", "Step-down Dunnett"))) stop("Procedure name is not recognized. ParAdjP function supports only the single-step Dunnett and step-down Dunnett procedures")

        # Number of procedures specified
        nproc <- length(proc)

        # Set up matrix to contain adjusted p-values
        adjp <- matrix(0,m,nproc)
        dimnames(adjp) <- list(NULL, paste(proc, ".adj.pvalue", sep=""))

	if (is.element("Single-step Dunnett", proc))
        {
        #for (i in 1:m) adjp[i]<-1-pdunnett(stat[i],nu,m)
        adjp[,"Single-step Dunnett.adj.pvalue"] <- sapply(stat, function(x) {1-pdunnett(x,nu,m)})
        }
	if (is.element("Step-down Dunnett", proc))
        {
        adjptmp <- rep(NA, length(stat))
		# Sort test statistics from largest to smallest
		or<-order(stat,decreasing=TRUE)
		sorted<-stat[or]
		for (i in 1:m)
		{
			if (i==1)
			{
				#adjp[1,"Step-down Dunnett.adj.pvalue"] <- 1-pdunnett(sorted[1],nu,m)
				adjptmp[1] <- 1-pdunnett(sorted[1],nu,m)
                                #print(1-pdunnett(sorted[1],nu,m))
				#maxp<-adjp[1, "Step-down Dunnett.adj.pvalue"]
				maxp<-adjptmp[1]
			}
			if (i>1 & i<m)
			{
				#adjp[i, "Step-down Dunnett.adj.pvalue"] <- max(maxp,1-pdunnett(sorted[i],nu,m-i+1))
				adjptmp[i] <- max(maxp,1-pdunnett(sorted[i],nu,m-i+1))

				#maxp<-max(maxp,adjp[i, "Step-down Dunnett.adj.pvalue"])
				maxp<-max(maxp,adjptmp[i])
			}
			if (i==m) {
                            #adjp[m, "Step-down Dunnett.adj.pvalue"] <- max(maxp,1-pt(sorted[m],nu))
                            adjptmp[m] <- max(maxp,1-pt(sorted[m],nu))
                        }
		}
		# Return to original ordering
		temp<-adjptmp
		adjptmp[or]<-temp
                adjp[,"Step-down Dunnett.adj.pvalue"] <- adjptmp
	}

	# Data frame returned by the function
	result<-data.frame(stat, round(rawp, 4), adjp)
	names(result)[1]<-"Test.statistic"
	names(result)[2]<-"Raw.pvalue"
	#names(result)[3]<-"Adj.pvalue"

	return(result=result)
}
# End of paradjp
# The ParCI function computes one-sided multiplicity-adjusted confidence
# intervals (simultaneous confidence intervals) for the single-step and
# step-down Dunnett procedures in one-sided testing problems with
# a balanced one-way layout and equally weighted null hypotheses
parci<-function(stat,n,est,stderror,covprob=0.975,proc)
# STAT, Vector of test statistics
# N, Common sample size in each treatment group
# EST, Vector of point estimates
# STDERROR, Vector of standard errors associated with the point estimates
# COVPROB, Simultaneous coverage probability
# PROC, Procedure name
    {
    # Number of null hypotheses
    m<-length(stat)

    if (m==0) stop("No test statistics are specified")

    if (m!=length(est)) stop("STAT and EST vectors have different lengths")
    if (m!=length(stderror)) stop("STAT and STDERROR vectors have different lengths")

    if (covprob>=1) stop("Simultaneous coverage probability must be <1")
    if (covprob<=0) stop("Simultaneous coverage probability must be >0")

    if(!all(proc %in% c("Single-step Dunnett", "Step-down Dunnett"))) stop("Procedure name is not recognized. ParCI function supports only the single-step Dunnett and step-down Dunnett procedures")
    
    if (n<=0) stop("Sample size must be positive")

    # number of procedures specified
    nproc <- length(proc)

    # set up matrix to contain confidence limits
    cimat <- matrix(0,m,nproc)
    dimnames(cimat) <- list(NULL, paste(proc, ".conf.limit", sep=""))

    # set up matrix to contain adjusted p-values
    adjpmat <- matrix(0,m,nproc)
    dimnames(adjpmat) <- list(NULL, paste(proc, ".adj.pvalue", sep=""))

    # Degrees of freedon
    nu<-(m+1)*(n-1)

    # Compute adjusted p-values
    result <- paradjp(stat,n,proc)

    #adjpmat <- result[, grep(".adj.pvalue", names(result), value=TRUE)]

    # One-sided familywise error rate
    alpha<-1-covprob

    # Rejection/acceptance of null hypotheses
    #reject<-(adjp<=alpha)

    # Vectors of confidence limits
    ci<-rep(0,m)

    zero<-rep(0,m)

	if (is.element("Single-step Dunnett", proc)) {
            adjpmat[, "Single-step Dunnett.adj.pvalue"] <- round(result[, "Single.step.Dunnett.adj.pvalue"], 4)

            reject <- (result[, "Single.step.Dunnett.adj.pvalue"] <= alpha)
           # Critical value

           c<-qdunnett(1-alpha,nu,m)
           cimat[, "Single-step Dunnett.conf.limit"] <-round(est-c*stderror, 4)
        }

	if (is.element("Step-down Dunnett", proc)) {
            adjpmat[, "Step-down Dunnett.adj.pvalue"] <- round(result[, "Step.down.Dunnett.adj.pvalue"], 4)

            reject <- (result[, "Step.down.Dunnett.adj.pvalue"] <= alpha)

	    # All null hypotheses are rejected
  	    if (sum(reject)==m) {
               # Critical value
               c<-qt(1-alpha,nu)
               cimat[, "Step-down Dunnett.conf.limit"] <- round(pmax(zero,est-c*stderror), 4)
            }

            # Some null hypotheses are accepted
  	    if (sum(reject)<m) {
                for (i in 1:m) {
                   if (reject[i]==1) cimat[i, "Step-down Dunnett.conf.limit"]<-0
                   if (reject[i]==0) {
                      # Critical value
                      c<-qdunnett(1-alpha,nu,m-sum(reject))
                      cimat[i, "Step-down Dunnett.conf.limit"] <- round(est[i]-c*stderror[i],4)
                   }
                }
            }
        }


    # Data frame returned by the function
    result<-data.frame(stat, est, stderror, adjpmat, cimat)
    names(result)[1]<-"Test.statistic"
    names(result)[2]<-"Estimate"
    names(result)[3]<-"Std.error"
    #names(result)[4]<-"Adj.pvalue"
    #names(result)[5]<-"Conf.limit"

    return(result=result)

    }
# End of parci

# ParGateAdjP function computes adjusted p-values and generates decision rules
# for multistage parallel gatekeeping procedures in hypothesis testing problems
# with multiple families of null hypotheses (null hypotheses are assumed
# to be equally weighted within each family)
pargateadjp<-function(gateproc, independence, alpha=0.05, printDecisionRules=FALSE)
# GATEPROC, List of gatekeeping procedure parameters
# INDEPENDENCE, Boolean indicator (TRUE, Independence condition is imposed; FALSE,
# Independence condition is not imposed)
# ALPHA: Global familywise error rate
# PRINTDECISIONRULES: Boolean indicator which controls printing of decision rules
{
	# Number of families
	nfams<-length(gateproc)
	if (nfams<=1) stop("Function requires more than one family of null hypotheses")

	for (i in 1:nfams)
	{
		pr<-gateproc[[i]]$proc
		if (pr!="Bonferroni" & pr!="Holm" & pr!="Hommel" & pr!="Hochberg" & pr!="Fallback")
			stop("Procedure name is not recognized. ParGateAdjP function supports only the Bonferroni, Holm, Hommel, Hochberg and fallback procedures")
	}

	if (alpha <= 0) stop("Alpha must be positive")
	if (alpha >= 1) stop("Alpha must be less than 1")

	temp<-gateproc

	for (i in 1:nfams)
	{
		# Number of null hypotheses
		nhyps<-length(temp[[i]]$rawp)
		adjp<-rep(0,nhyps)
		# Placeholder for adjusted p-values
		gateproc[[i]][5]<-list(adjp=adjp)

		for (j in 1:nhyps)
		{

			# Find the lowest alpha level at which the current null hypothesis is rejected
			upper<-1
			lower<-0
			for (k in 1:20)
			{
				current<-(lower+upper)/2
				# Evaluate decision rules for the multistage parallel gatekeeping procedure
				res<-pargateeval(temp,current,independence)
				# Rejection decision for the current null hypothesis
				if (independence==TRUE | i==nfams) reject<-res[[i]][[7]][j]
				# Rejection decisions after retesting if the independence condition is not imposed
				if (independence==FALSE & i<nfams)
				{
					# If the current null hypothesis was retested
					if (res[[2*nfams-i]][[6]]>0) reject<-res[[2*nfams-i]][[7]][j] else reject<-res[[i]][[7]][j]
				}
				# Update the interval
				if (reject==TRUE) upper<-current
				if (reject==FALSE) lower<-current
			}

			# Global adjusted p-value
			gateproc[[i]][[5]][j]<-(lower+upper)/2
		}
	}

	# Build a data frame with the raw and global adjusted p-values
	count<-0
	for (i in 1:nfams) {
		count<-count+length(gateproc[[i]]$rawp)
	}
        result <- data.frame()
	k<-1
	for (i in 1:nfams)
	{
		# Number of null hypotheses
		nhyps<-length(gateproc[[i]]$rawp)
		for (j in 1:nhyps)
		{
			result[k,1]<-gateproc[[i]]$label
			result[k,2]<-gateproc[[i]]$proc
			result[k,3]<-gateproc[[i]]$procpar
			result[k,4]<-round(gateproc[[i]]$rawp[j], 4)
			result[k,5]<-round(gateproc[[i]][[5]][j], 4)
			k<-k+1
		}
	}
	names(result)[1]<-"Family"
	names(result)[2]<-"Procedure"
	names(result)[3]<-"Parameter"
	names(result)[4]<-"Raw.pvalue"
	names(result)[5]<-"Adj.pvalue"

	if(printDecisionRules==TRUE) { pargaterule(gateproc,alpha,independence)}

	return(result=result)
}
# End of pargateadjp
# ParGateRule function is a secondary function which generates decision rules for multistage parallel
# gatekeeping procedures in hypothesis testing problems with multiple families
# of null hypotheses (null hypotheses are assumed to be equally weighted within
# each family)
pargaterule<-function(gateproc,alpha,independence)
# GATEPROC, List of gatekeeping procedure parameters
# ALPHA, Global familywise error rate
# INDEPENDENCE, Boolean indicator (TRUE, Independence condition is imposed; FALSE,
# Independence condition is not imposed)
{

    # Number of families
    nfams<-length(gateproc)

	# Evaluate decision rules for the multistage parallel gatekeeping procedure
	gateproc<-pargateeval(gateproc,alpha,independence)

	cat("Hypothesis testing problem\n\n")
	cat("Global familywise error rate=", alpha, "\n",sep="")
	if (independence==TRUE) cat("Independence condition is imposed (the families are tested from first to last) \n",sep="")
	if (independence==FALSE) cat("Independence condition is not imposed (the families are tested from first to last and then re-tested from last to first) \n",sep="")
	cat("\n\n")

	# Retesting indicator
	retest<-FALSE

    # Sequentially numbered hypotheses labels
    hyplabel<-vector("list",nfams)
    cumtotal<-0

	# Test families from first to last
	for (i in 1:nfams)
	{
		# Local alpha level
		la<-round(gateproc[[i]][[6]],4)
		# Procedure parameter
		pp<-round(gateproc[[i]]$procpar,4)

		# Family and procedure
		cat("Family ", i, " (", gateproc[[i]]$label, ") is tested using", sep="")
		cat(" ", gateproc[[i]]$proc, " procedure (truncation parameter=", pp, ") at alpha", i, "=", la, ".\n\n", sep="")

		# Number of null hypotheses
		nhyps<-length(gateproc[[i]]$rawp)
        # Hypotheses labels
        hyplabel[[i]]<-seq(1:nhyps)+cumtotal
        cumtotal<-cumtotal+nhyps
		# Number of rejected null hypotheses
		rejcount<-sum(gateproc[[i]][[7]])

		for (j in 1:nhyps)
		{
			# Raw p-value
			rp<-round(gateproc[[i]][[2]][j],4)
			# Adjusted p-value
			ap<-round(gateproc[[i]][[5]][j],4)

			cat("Null hypothesis ", hyplabel[[i]][j], " (raw p-value=", rp, ")", sep="")
			if (gateproc[[i]][[7]][j]==TRUE)
			{
				cat(" is rejected.\n\n", sep="")
			}
			if (la==0) cat(" is automatically accepted. \n\n", sep="")
			if (la>0 & gateproc[[i]][[7]][j]==FALSE) cat(" is accepted.\n\n", sep="")

		}

        # Details
        cat("Details on the decision rule for this family can be obtained by running the PValAdjP function for ",gateproc[[i]]$proc," procedure with gamma=",pp," and alpha=",la,".\n\n", sep="")

		# Consclusion
		if (rejcount==0 & i<nfams)
		{
			cat("No null hypotheses are rejected in Family ", i, " and the parallel gatekeeping procedure cannot pass this family.",sep="")
			cat(" Testing stops and all remaining null hypotheses are automatically accepted.\n\n\n")
		}

		if (rejcount>0 & i<nfams)
		{
			cat("One or more null hypotheses are rejected in Family ", i, " and the parallel gatekeeping procedure passes this family.",sep="")
			cat(" Based on the error rate function of ", gateproc[[i]]$proc, " procedure (truncation parameter=", pp, "),", sep="")
			cat(" alpha",i+1,"=",round(gateproc[[i+1]][[6]],4)," is carried over to Family ", i+1, ".\n\n\n",sep="")
		}

		if (i==nfams & independence==FALSE)
		{
			if (rejcount==nhyps)
			{
				retest<-TRUE
				cat("All null hypotheses are rejected in Family ", i, " and the parallel gatekeeping procedure passes this family.", sep="")
				cat(" Retesting begins and alpha",i+1,"=",round(alpha,4)," is carried over to Family ", i-1, ".\n\n\n",sep="")
			}
			if (rejcount<nhyps)
			{
				cat("Some null hypotheses are accepted in Family ", i, " and the parallel gatekeeping procedure cannot pass this family.", sep="")
				cat(" Retesting will not be performed.\n\n\n",sep="")
			}
		}
	}

	# Retest families from last to first if the independence condition is not imposed
	if (independence==FALSE & retest==TRUE)
	{
		for(i in (nfams+1):(2*nfams-1))
		{
			# Family index
			k<-2*nfams-i
			# Local alpha level
			la<-round(gateproc[[i]][[6]],4)
			# Procedure parameter
			pp<-round(gateproc[[i]]$procpar,4)

			cat("Family ", k, " (", gateproc[[k]]$label, ") is retested using", sep="")
			cat(" ", gateproc[[i]]$proc, " procedure (truncation parameter=", pp, ") at alpha", i, "=", la, ".\n\n", sep="")

			# Number of null hypotheses
			nhyps<-length(gateproc[[i]]$rawp)
			# Number of rejected null hypotheses
			rejcount<-sum(gateproc[[i]][[7]])

			for (j in 1:nhyps)
			{
				# Raw p-value
				rp<-round(gateproc[[i]][[2]][j],4)
				# Adjusted p-value
				ap<-round(gateproc[[i]][[5]][j],4)

				cat("Null hypothesis ", hyplabel[[k]][j], " (raw p-value=", rp, ")", sep="")
				if (gateproc[[i]][[7]][j]==TRUE)
				{
					cat(" is rejected.\n\n", sep="")
				}
				if (la==0) cat(" is automatically accepted.\n\n", sep="")
				if (la>0 & gateproc[[i]][[7]][j]==FALSE) cat(" is accepted.\n\n", sep="")
			}

            # Details
            cat("Details on the decision rule for this family can be obtained by running the PValAdjP function for ",gateproc[[i]]$proc," procedure with gamma=",pp," and alpha=",la,".\n\n", sep="")

            # Conclusions
			if (k>1)
			{
				if (rejcount==nhyps)
				{
					retest<-TRUE
					cat("All null hypotheses are rejected in Family ", k, " and the parallel gatekeeping procedure passes this family and", sep="")
					cat(" alpha",i+1,"=",round(alpha,4)," is carried over to Family ", k-1, ".\n\n\n",sep="")
				}
				if (rejcount<nhyps)
				{
					cat("Some null hypotheses are accepted in Family ", k, " and the parallel gatekeeping procedure cannot pass this family.", sep="")
					cat(" Retesting stops.\n\n\n",sep="")
				}
			}
		}
	}

}
# End of pargaterule

# PValAdjP function computes adjusted p-values and generates decision rules
# for the Bonferroni, Holm, Hommel, Hochberg, fixed-sequence and fallback procedures
pvaladjp<-function(rawp,weight=rep(1/length(rawp), length(rawp)),alpha=0.05,
		           proc=c("Bonferroni", "Holm", "Hommel", "Hochberg", "Fixed-sequence", "Fallback"),
				   printDecisionRules=FALSE)
# RAWP, Vector of raw p-values
# WEIGHT, Vector of hypothesis weights
# ALPHA, Familywise error rate
# PROC, Procedure name
# PRINTDECISIONRULES: Boolean indicator which controls printing of decision rules

{

	# Number of null hypotheses
	m<-length(rawp)

	if (m==0) stop("No p-values are specified")

	for (i in 1:m)
	{
		if (rawp[i]<0) stop("P-values must be positive")
		if (rawp[i]>1) stop("P-values must be less than 1")
	}

	index <-order(rawp)

	if (alpha <= 0) stop("Alpha must be positive")
	if (alpha >= 1) stop("Alpha must be less than 1")

	# Number of weights
	nweis<-length(weight)

	if (m!=nweis) stop("RAWP and WEIGHT vectors have different lengths")

	if (sum(weight)>1) stop("Sum of hypothesis weights must be <=1")

	for (i in 1:nweis)
	{
		if (weight[i]<0) stop("Hypothesis weights must be >=0")
	}

	if(!all(proc %in% c("Bonferroni", "Holm", "Hommel", "Hochberg", "Fixed-sequence", "Fallback")))
		stop("Procedure name is not recognized. PValAdjP function supports only the Bonferroni, Holm, Hommel, Hochberg, Fixed-sequence, and Fallback procedures")

	# number of procedures specified
	nproc <- length(proc)

	# set up matrix to contain adjusted p-values
	adjp <- matrix(0,m,nproc)
	dimnames(adjp) <- list(NULL, paste(proc, ".adj.pvalue", sep=""))

	if (is.element("Bonferroni", proc)) {
		adjp[,"Bonferroni.adj.pvalue"]<-pvaltrunc(rawp,weight,"Holm",0)
	}
	if (is.element("Holm", proc)) {
		adjp[,"Holm.adj.pvalue"]<-pvaltrunc(rawp,weight,"Holm",1)
	}
	if (is.element("Hommel", proc)) {
		adjp[,"Hommel.adj.pvalue"]<-pvaltrunc(rawp,weight,"Hommel",1)
	}
	if (is.element("Hochberg", proc)) {
		adjp[,"Hochberg.adj.pvalue"]<-pvaltrunc(rawp,weight,"Hochberg",1)
	}
	if (is.element("Fixed-sequence", proc)) {
		adjp[,"Fixed-sequence.adj.pvalue"]<-pvaltrunc(rawp,weight,"Fixed-sequence",1)
	}
	if (is.element("Fallback", proc)) {
		adjp[,"Fallback.adj.pvalue"]<-pvaltrunc(rawp,weight,"Fallback",1)
	}

	# Data frame returned by the function
	result<-data.frame(round(rawp,4),round(weight,4), round(adjp,4))
	names(result)[1]<-"Raw.pvalue"
	names(result)[2]<-"Weight"

	if(printDecisionRules==TRUE) {
            if(length(unique(round(weight,3))) > 1) stop("Weights must be equal for decision rule calculations to be valid")
            pvalrule(rawp, alpha, proc)
        }

	return(result=result)
}
# End of pvaladjp

# The PvalCI function computes one-sided multiplicity-adjusted confidence
# intervals (simultaneous confidence intervals) for the Bonferroni, Holm and fixed-sequence
# procedures in general hypothesis testing problems with equally or
# unequally weighted null hypotheses
pvalci<-function(rawp,est,stderror,weight=rep(1/length(rawp),length(rawp)),covprob=0.975,proc=c("Bonferroni", "Holm", "Fixed-sequence"))
# RAWP, Vector of raw p-values
# EST, Vector of point estimates
# STDERROR, Vector of standard errors associated with the point estimates
# WEIGHT, Vector of hypothesis weights
# COVPROB, Simultaneous coverage probability
# PROC, Procedure name
    {

    # Number of null hypotheses
    m<-length(rawp)

    if (m==0) stop("No p-values are specified")

    for (i in 1:m)
        {
        if (rawp[i]<0) stop("P-values must be positive")
        if (rawp[i]>1) stop("P-values must be less than 1")
        }

    if (m!=length(weight)) stop("RAWP and WEIGHT vectors have different lengths")
    if (m!=length(est)) stop("RAWP and EST vectors have different lengths")
    if (m!=length(stderror)) stop("RAWP and STDERROR vectors have different lengths")

    if (sum(weight)>1) stop("Sum of hypothesis weights must be <=1")

    for (i in 1:length(weight))
        {
        if (weight[i]<0) stop("Hypothesis weights must be >=0")
        }

    if (covprob>=1) stop("Simultaneous coverage probability must be <1")
    if (covprob<=0) stop("Simultaneous coverage probability must be >0")

    if (!all(proc %in% c("Bonferroni", "Holm", "Fixed-sequence")))
        stop("Procedure name is not recognized")
    #if (proc!="Bonferroni" & proc!="Holm" & proc!="Fixed-sequence") stop("Procedure name is not recognized")

    # number of procedures specified
    nproc <- length(proc)

        # set up matrix to contain confidence limits
    cimat <- matrix(0,m,nproc)
    dimnames(cimat) <- list(NULL, paste(proc, ".conf.limit", sep=""))

 	# One-sided familywise error rate
 	alpha<-1-covprob

    # Compute adjusted p-values
    result <- pvaladjp(rawp=rawp,weight=weight,alpha=alpha,proc=proc)
    adjpmat <- result[, grep(".adj.pvalue", names(result), value=TRUE)]
    #adjp<-out[c(-1,-2)]
    #print(out)

        # Rejection/acceptance of null hypotheses
	#reject<-(adjp<=alpha)

    # Vectors of confidence limits
    ci<-rep(0,m)

 	zero<-rep(0,m)


   # adjp<-out[,3]

    # Bonferroni procedure
	if (is.element("Bonferroni", proc)) {
           reject <- (result[,"Bonferroni.adj.pvalue"] <= alpha)
           cimat[,"Bonferroni.conf.limit"] <-est-(stderror*qnorm(1-(alpha*weight)))
        }

    # Holm procedure
	if (is.element("Holm", proc)) {
            reject <- (result[,"Holm.adj.pvalue"] <= alpha)
        bonfci<-est-(stderror*qnorm(1-(alpha*weight)))
        # All null hypotheses are rejected
 		if (sum(reject)==m) cimat[,"Holm.conf.limit"] <-pmax(zero,bonfci)
        # Some null hypotheses are accepted
        if (sum(reject)<m)
            {
            for(i in 1:m)
                {
                if (reject[i]==1) cimat[i, "Holm.conf.limit"]<-0
                if (reject[i]==0)
                    {
                    adjalpha<-(alpha*weight[i])/sum(weight[reject==0])
                    cimat[i,"Holm.conf.limit"]<-est[i]-(stderror[i]*qnorm(1-adjalpha))
                    }
                }
            }
        }

   	# Fixed-sequence procedure
	if (is.element("Fixed-sequence", proc)) {
            reject <- (result[,"Fixed.sequence.adj.pvalue"] <= alpha)
		# All null hypotheses are accepted
  		if (sum(reject)==0)
            {
            cimat[1,"Fixed-sequence.conf.limit"] <-est[1]-stderror[1]*qnorm(1-alpha)
            for(i in 2:m) cimat[i, "Fixed-sequence.conf.limit"]<-NA
            }
        # All null hypotheses are rejected
		if (sum(reject)==m)
            {
            temp1<-est-stderror*qnorm(1-alpha)
            cimat[,"Fixed-sequence.conf.limit"] <-min(temp1)
            }
        # Some null hypotheses are accepted and some are rejected
        if (sum(reject)>0 & sum(reject)<m)
            {
            cimat[1, "Fixed-sequence.conf.limit"]<-0
            for(i in 2:m)
                {
                if (reject[i]==1) cimat[i, "Fixed-sequence.conf.limit"]<-0
                if (reject[i]==0 & reject[i-1]==1) cimat[i, "Fixed-sequence.conf.limit"]<-est[i]-stderror[i]*qnorm(1-alpha)
                if (reject[i]==0 & reject[i-1]==0) cimat[i, "Fixed-sequence.conf.limit"]<-NA
                }
            }
        }

	# Data frame returned by the function
    #result<-data.frame(rawp,est,stderror,weight,adjp,ci)
    result<-data.frame(rawp,est,stderror,weight,adjpmat,cimat)
    names(result)[1]<-"Raw.pvalue"
    names(result)[2]<-"Estimate"
    names(result)[3]<-"Std.error"
    names(result)[4]<-"Weight"
    #names(result)[5]<-"Adj.pvalue"
    #names(result)[6]<-"Conf.limit"

    return(result=result)
}
# End of pvalci
# PValRule function  is a secondary function which generates decision rules for the Bonferroni, Holm, 
# Hommel, Hochberg, fixed-sequence and fallback procedures in hypothesis testing problems with 
# equally weighted null hypotheses 
pvalrule<-function(rawp,alpha,proc)
# RAWP, Vector of raw p-values
# ALPHA, Familywise error rate
# PROC, Procedure name
{
	# Number of p-values
	m<-length(rawp)
	
	cat("Hypothesis testing problem\n\n")
	cat("Familywise error rate: alpha=",alpha,"\n",sep="")
	
	# List of original null hypotheses and raw p-values
	cat("Original null hypotheses: ",sep="")
	for (i in 1:m) cat("H",i," ",sep="")
	cat("(equally weighted)\n")    
	cat("Original raw p-values: ")
	for (i in 1:m) cat("p",i,"=",round(rawp[i],4)," ",sep="")        
	cat("\n")
	
	# List of ordered null hypotheses and raw p-values    
	o<-order(rawp)
	orderp<-rawp[o]        
	cat("Ordered null hypotheses: ",sep="")
	for (i in 1:m) cat("H(",i,") ",sep="")
	cat("\n")
	cat("Ordered raw p-values: ")
	for (i in 1:m) cat("p(",i,")=",round(orderp[i],4)," ",sep="")        
	cat("\n\n")
	
	# Print decision rules for each specified procedure
	proc.to.print <- proc
	for(p in proc.to.print) {
	
		# Bonferroni procedure (parameters needed: m, rawp, alpha)
		if (p=="Bonferroni")
		{
			cat("\n\n")
			cat("Decision rules for the Bonferroni procedure:\n\n",sep="")
			for (i in 1:m)
			{
				cat("Original null hypothesis H",i,sep="")
				level<-alpha/m
				if (rawp[i]<=level) cat(" is rejected since p",i,"<=alpha/",m,"=",round(level,4),"\n\n",sep="")
				else cat(" is accepted since p",i,">alpha/",m,"=",round(level,4),"\n\n",sep="")
			}
		}    
        
        # Fixed-sequence procedure    
        if (p=="Fixed-sequence")
            {
            cat("\n\n")
            cat("Decision rules for the fixed-sequence procedure:\n\n",sep="")        
            for (i in 1:m)
                {
                cat("Step ",i,"\n",sep="")
                cat("Original null hypothesis H",i,sep="")
                if (i==1)
                    {
                    # Current null hypothesis is rejected                
                    if (rawp[1]<=alpha) 
                        {
                        current<-1
                        cat(" is rejected since p1<=alpha=",round(alpha,4),"\n\n",sep="")
                        }
                    # Current null hypothesis is accepted
                    if (rawp[1]>alpha) 
                        {
                        current<-0                 
                        cat(" is accepted since p1>alpha=",round(alpha,4),"\n\n",sep="")
                        }
                    }
                if (i>1)
                    {
                    # Preceding null hypothesis
                    prec<-current
                    # Current null hypothesis is rejected                
                    if (rawp[i]<=alpha & prec==1) current<-1
                    # Current null hypothesis is accepted
                    if (rawp[i]>alpha | prec==0) current<-0                
                    if (rawp[i]<=alpha & prec==1) cat(" is rejected since p",i,"<=alpha=",round(alpha,4)," and the preceding original null hypothesis is rejected\n\n",sep="")
                    if (rawp[i]>alpha & prec==1) cat(" is accepted since p",i,">alpha=",round(alpha,4),"\n\n",sep="")
                    if (prec==0) cat(" is accepted since the preceding original null hypothesis is accepted\n\n",sep="")                		       
                }            
                }
            }
            
        # Fallback procedure    
        if (p=="Fallback")
            {
            cat("\n\n")
            cat("Decision rules for the fallback procedure:\n\n",sep="")
            numer<-1
            for (i in 1:m)
                {
                cat("Step ",i,"\n",sep="")
                cat("Original null hypothesis H",i,sep="")            
                level<-numer*alpha/m
                if (i<m)
                    {                               
                    # Current null hypothesis is rejected                
                    if (rawp[i]<=level) 
                        {
                        cat(" is rejected since p",i,"<=",numer,"*alpha/",m,"=",round(level,4)," and the significance level used in this test is carried over to the next original null hypothesis\n\n",sep="")
                        numer<-numer+1
                        }
                    # Current null hypothesis is accepted
                    if (rawp[i]>level) 
                        {
                        current<-0 
                        cat(" is accepted since p",i,">",numer,"*alpha/",m,"=",round(level,4),"\n\n",sep="")
                        numer<-1
                        }
                    }
                if (i==m)
                    {
                    if (rawp[i]<=level) cat(" is rejected since p",i,"<=",numer,"*alpha/",m,"=",round(level,4),"\n\n",sep="")
                    else cat(" is accepted since p",i,">",numer,"*alpha/",m,"=",round(level,4),"\n\n",sep="")                    
                    }
                }        
            }                    
        
		# Holm procedure (parameters needed: m, rawp, orderp, alpha)
		if (p=="Holm")
		{
			cat("Decision rules for the Holm procedure:\n\n",sep="")
			for (i in 1:m)
			{
				cat("Step ",i,"\n",sep="")
				cat("Ordered null hypothesis H(",i,") [Original null hypothesis H",o[i],"]",sep="")
				level<-alpha/(m-i+1)
				if (i==1)
				{                
					# Current null hypothesis is rejected                
					if (orderp[1]<=level) current<-1
					# Current null hypothesis is accepted
					if (orderp[1]>level) current<-0                
					if (orderp[1]<=level) cat(" is rejected since p(1)<=alpha/",m,"=",round(level,4),"\n\n",sep="")
					else cat(" is accepted since p(1)>alpha/",m,"=",round(level,4),"\n\n",sep="")
				}
				if (i>1)
				{
					# Preceding null hypothesis
					prec<-current
					# Current null hypothesis is rejected                
					if (orderp[i]<=level & prec==1) current<-1
					# Current null hypothesis is accepted
					if (orderp[i]>level | prec==0) current<-0                
					if (orderp[i]<=level & prec==1) cat(" is rejected since p(",i,")<=alpha/",m-i+1,"=",round(level,4)," and the preceding ordered null hypothesis is rejected\n\n",sep="")
					if (orderp[i]>level & prec==1) cat(" is accepted since p(",i,")>alpha/",m-i+1,"=",round(level,4),"\n\n",sep="")
					if (prec==0) cat(" is accepted since the preceding ordered null hypothesis is accepted\n\n",sep="")                		       
				}
			}
		}  
		
		# Hochberg procedure (parameters needed: m, rawp, orderp, alpha)
		if (p=="Hochberg")
		{
			cat("Decision rules for the Hochberg procedure:\n\n",sep="")
			for (i in 1:m)
			{
				cat("Step ",i,"\n",sep="")            
				level<-alpha/i
				if (i==1)
				{                
					if (orderp[m]<=level) 
					{
						cat("Ordered null hypothesis H(",m,") [Original null hypothesis H",o[m],"]",sep="")
						cat(" is rejected since p(",m,")<=alpha/1=",round(level,4),"\n\n",sep="")  
						if (m>1)    
						{
							for (j in (m-1):1)
							{
								cat("Ordered null hypothesis H(",j,") [Original null hypothesis H",o[j],"]",sep="")
								cat(" is rejected since the preceding ordered null hypothesis is rejected\n\n",sep="")                        
							}
						}    
						break
					}
					if (orderp[m]>level)
					{
						cat("Ordered null hypothesis H(",m,") [Original null hypothesis H",o[m],"]",sep="")
						cat(" is accepted since p(",m,")>alpha/1=",round(level,4),"\n\n",sep="")
					}
				}
				if (i>1)
				{
					if (orderp[m-i+1]<=level) 
					{
						cat("Ordered null hypothesis H(",m-i+1,") [Original null hypothesis H",o[m-i+1],"]",sep="")
						cat(" is rejected since p(",m-i+1,")<=alpha/",i,"=",round(level,4),"\n\n",sep="")
						if (m-i+1>1)
						{
							for (j in (m-i):1)
							{
								cat("Ordered null hypothesis H(",j,") [Original null hypothesis H",o[j],"]",sep="")
								cat(" is rejected since the preceding ordered null hypothesis is rejected\n\n",sep="")                        
							}
						}    
						break    
					}
					if (orderp[m-i+1]>level)
					{
						cat("Ordered null hypothesis H(",m-i+1,") [Original null hypothesis H",o[m-i+1],"] ",sep="")
						cat("is accepted since p(",m-i+1,")>alpha/",i,"=",round(level,4),"\n\n",sep="")
					}    
				}
			}
		}         
		
		# Hommel procedure (parameters needed: m, rawp, orderp, alpha)
		if (p=="Hommel")
		{
            cat("Decision rules for the Hommel procedure:\n\n",sep="")
            for (i in 1:m)
                {
                cat("Step ",i,"\n",sep="")            
                if (i==1)
                    {                
                    if (orderp[m]<=alpha) 
                        {
                        cat("Ordered null hypothesis H(",m,") [Original null hypothesis H",o[m],"]",sep="")
                        cat(" is rejected since p(",m,")<=alpha/1=",round(alpha,4),"\n\n",sep="")  
                        if (m>1)    
                            {
                            for (j in (m-1):1)
                                {
                                cat("Ordered null hypothesis H(",j,") [Original null hypothesis H",o[j],"]",sep="")
                                cat(" is rejected since the preceding ordered null hypothesis is rejected\n\n",sep="")                            }
                            }    
                        break   
                        }
                    if (orderp[m]>alpha)
                        {
                        cat("Ordered null hypothesis H(",m,") [Original null hypothesis H",o[m],"]",sep="")
                        cat(" is accepted since p(",m,")>alpha/1=",round(alpha,4),"\n\n",sep="")
                        }
                    }
                if (i>1)
                    {
                    # All p-values are greater than corresponding significance levels
                    all<-1
                    hyprej<-0
                    jrej<-0
                    irej<-0
                    for (j in 1:i)
                        {
                        if (orderp[m-i+j]<=j*alpha/i) 
                            {
                            if (hyprej==0)
                                {
                                all<-0 
                                hyprej<-m-i+j 
                                jrej<-j 
                                irej<-i
                                }
                            }
                        }
                    if (all==1)
                        {
                        cat("Ordered null hypothesis H(",m-i+1,") [Original null hypothesis H",o[m-i+1],"] ",sep="")
                        cat("is accepted since",sep="")
                        for (j in 1:(i-1)) cat(" p(",m-i+j,")>",j,"*alpha/",i,"=",round(j*alpha/i,4)," and",sep="")
                        cat(" p(",m,")>alpha/1=",round(alpha,4),"\n\n",sep="")
                        }
                    if (all==0)
                        {
                        cat("p(",hyprej,")<=",jrej,"*alpha/",irej,"=",round(jrej*alpha/irej,4),sep="")
                        cat(" and thus all remaining hypotheses are tested with alpha/",i-1,"=",round(alpha/(i-1),4),": \n\n",sep="")                      
                        allsub<-1
                        for (j in (m-i+1):1)
                            {                      
                            cat("Ordered null hypothesis H(",j,") [Original null hypothesis H",o[j],"]",sep="")
                            if (orderp[j]<=alpha/(i-1)) 
                                {
                                if (allsub==1)
                                    {
                                    allsub<-0
                                    cat(" is rejected since p(",j,")<=alpha/",i-1,"=",round(alpha/(i-1),4),"\n\n",sep="") 
                                    }
                                else cat(" is rejected since the preceding ordered null hypothesis is rejected\n\n",sep="")
                                }
                            if (orderp[j]>alpha/(i-1)) cat(" is accepted since p(",j,")>alpha/",i-1,"=",round(alpha/(i-1),4),"\n\n",sep="")                              
                            }
                        break    
                        }
                    }
                }
			
		}         
		
		cat("\n\n")
	}
}
# End of pvalrule# Multxpert package implements commonly used p-value-based and parametric
# multiple testing procedures and parallel gatekeeping procedures

# For more information about the Multxpert package, visit the Multiplicity Expert web site
# http://multxpert.com/wiki/MultXpert_package

.onLoad <- function(lib, pkg) {
	if (interactive()) {
		cat('multxpert: Implementation of commonly used p-value based and parametric\n',
				'multiple testing procedures and parallel gatekeeping procedures.\n',
				'For more information visit http://multxpert.com/wiki/MultXpert_package\n', sep='')
	}
}