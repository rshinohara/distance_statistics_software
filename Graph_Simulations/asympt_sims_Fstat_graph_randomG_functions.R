
library(graph)
library(gRbase)
library(R.utils)
library(vegan)
V<-LETTERS[1:5]

# setwd('~/Dropbox/Documents/Distance_Statistics/') # laptop
# setwd('/Volumes/RAID Array/Dropbox/Documents/Distance_Statistics/') # desktop
setwd('~/Distance_Statistics/') # cluster

resetPar <- function() {
    dev.new()
    op <- par(no.readonly = TRUE)
    dev.off()
    op
}

#################################################
# FUNCTIONS FOR GRAPH CASES
#################################################

# This function is a Monte Carlo version of combn.
combn.random<-function(vec,n,B.mc) {
	if (B.mc > choose(length(vec),n)) {
		print('ERROR - TOO MANY COMBINATIONS REQUESTED ...')
		return(NULL)
	} else {
		num.combns<-0 # number of combinations found
		combn.mat<-c()

		# until enough unique combinations have been found
		while (num.combns<B.mc) {
			combn.mat<-rbind(combn.mat,t(replicate(B.mc-num.combns,sort(sample(vec,n,replace=FALSE)))))
			duplicate.rows<-c(which(duplicated(data.frame(combn.mat))))
			if (length(duplicate.rows)>0) {
				combn.mat<-combn.mat[-duplicate.rows,]
			}
			num.combns<-dim(combn.mat)[1]
		}
		return(t(combn.mat[1:B.mc,]))
	}
}

#This function simulates a graph with a given set of nodes and with a symmetric probability matrix showing 
random.graph<-function(V,pmat) {
	if ((!isSymmetric(pmat,tol=1E-5))|(length(V)!=(dim(pmat)[1]))) { print("ERROR - invalid inputs... ") 
	return(NA)} else {
		edList <- vector("list", length=length(V))
		names(edList) <- V

		#Sample edges
		get.binom<-function(p) rbinom(1,1,p)
		pvec<-pmat[upper.tri(pmat)]; edVec<-sapply(pvec,get.binom)
		edMat<-pmat; ind<-lower.tri(pmat); edMat[upper.tri(pmat)]<-edVec;  edMat[ind]<-t(edMat)[ind]
	
		
		#Find list of edges
		for (i in 1:length(V))
			edList[[i]] <- list(edges=V[-i][which(edMat[i,-i]==1)], weights=rep(1,sum(edMat[i,-i]==1)))

		#Return the graph
		gR <- graphNEL(nodes=V, edgeL=edList)
		return(gR)
		}
}




### THIS DOES NOT ALLOW DUPLICATES
# This function estimates the distance-based sample variance in a list of graphs based on B pairs of graphs.
var.graph.mc<-function(sample.g,d,B.mc=1000,n.cores.=n.cores,exact=FALSE) {
	n<-length(sample.g)
	if ((exact)||(B.mc > choose(length(sample.g),2))) {
		xmat<-combn(1:length(sample.g),2)
	} else {
		xmat<-combn.random(1:length(sample.g),2,B.mc)
	}
	sample.dist<-mapply(d,sample.g[xmat[1,]],sample.g[xmat[2,]])
	mean(sample.dist)
}

anova.graph.SS<-function(sample.g1,sample.g2,d=edge.diff,B.mc=1E4,n.cores.=n.cores,exact=FALSE) {
	return(var.graph.mc(c(sample.g1,sample.g2),d,B.mc,n.cores.,exact)*(length(sample.g1)+length(sample.g2)-1))
}

anova.graph.SSE<-function(sample.g1,sample.g2,d=edge.diff,B.mc=1E4,n.cores.=n.cores,exact=FALSE) {
	SSE<-((length(sample.g1)-1)*var.graph.mc(sample.g1,d,B.mc,n.cores.,exact)+(length(sample.g2)-1)*var.graph.mc(sample.g2,d,B.mc,n.cores.,exact))
	return(SSE)
}

# This function estimates the distance between a subject and a list of other subjects.
get.d.list<-function(x,x.list,d) lapply(x.list,d,x)

# Cuts the time in half.
get.d.mat.list<-function(x.list1,x.list2,d) {
	x.list<-c(x.list1,x.list2);n.list<-length(x.list)
	mat.out<-matrix(0,n.list,n.list)
	for (jj in 1:(n.list-1)) {
		mat.out[jj,(jj+1):n.list]<-unlist(get.d.list(x.list[[jj]],x.list[(jj+1):n.list],d))
	}
	mat.out<-mat.out+t(mat.out)
	return(mat.out)
} 

get.A.mat<-function(dists.mat,d,eta,sig) {
	exp.dist<- apply(dists.mat,1,mean) # This is the expected distance from each subject.
	A.term1<-(1-1/eta*(d%o%d)-1/(1-eta)*((1-d)%o%(1-d)))*dists.mat
	A.term2<-(1-(1-d)/eta)%o%((1-d)*exp.dist-sig^2*eta);A.term2<-A.term2+t(A.term2)
	A.term3<-(1-d/(1-eta))%o%(d*exp.dist-sig^2*(1-eta));A.term3<-A.term3+t(A.term3)
	A.term4<-(1-d-eta)%o%(1-d-eta)*sig^2*eta^(-1)*(1-eta)^(-1)
	return(A.term1-A.term2-A.term3-A.term4)
}


#This function does the graph simulation and returns the simulated data.
simulate.data<-function(n,p,tau,d) {

	dd<-rbinom(n,1,1-p)
	n.2<-sum(dd);n.1<-n-n.2

 	pmat<-cbind(c(1,1,1-tau,tau, tau),c(1,1, tau, tau, tau),c(1-tau, tau,1,1,1),c(tau, tau,1,1,1),c(tau, tau,1,1,1))
 	sample.g1<-replicate(n.1,random.graph(V,pmat),simplify=FALSE)
 	sample.g2<-replicate(n.2,random.graph(V,pmat),simplify=FALSE)
	
 	return(list(dd=dd,n.1=n.1,n.2=n.2,sample.g1=sample.g1,sample.g2=sample.g2))

}

#This function does the proposed testing.
do.proposed.test<-function(n,d,dd,n.1,n.2,sample.g1,sample.g2,d.mat,B.mc,exact=FALSE,return.details=FALSE,return.sample=TRUE,return.sample.size=5000) {

	SS<-anova.graph.SS(sample.g1,sample.g2,d=d,B.mc=B.mc,n.cores=1,exact)
	SSE<-anova.graph.SSE(sample.g1,sample.g2,d=d,B.mc=B.mc,n.cores=1,exact)

	numer<-SS-SSE
	denom<-SSE/(n-2)
	ganova.F<-numer/denom

	A.mat<-get.A.mat(d.mat,c(rep(0,n.1),rep(1,n.2)),n.1/n,sqrt(SS/n))
	e.A.mat<-eigen(A.mat)
	evals<-e.A.mat$values[order(abs(e.A.mat$values),decreasing=TRUE)][1:10]
	numer.approx<-get.ganova.approx(evals,K=10,N=B.mc)/n+(SS/n)
	p.val<-mean(ganova.F<numer.approx/denom)

### RETURN RESULTS
	if (return.details==TRUE) {
		return(list(SS=SS,SSE=SSE,numer=numer,denom=denom,ganova.F=ganova.F,evals=evals,numer.approx=numer.approx,p.val=p.val))
	} else {
		if (return.sample==TRUE) {
			return(list(ganova.F=ganova.F,p.val=p.val,asympt.est.sample=numer.approx[1:return.sample.size]/denom))
		} else {
			return(list(ganova.F=ganova.F,p.val=p.val))
		}
	}

}

get.ganova.approx<-function(evals,K,N=1E4){
	tmp<-rep(0,N)
	for (j in 1:K) tmp<-tmp+(evals[j]*(rchisq(N,1)-1))
	return(tmp)
}

#This function takes a sample from the null and compares the observed data to it.
is.rejected<-function(x,vec) return(x>quantile(unlist(vec),probs=0.95))

pos.part<-function(x) x*(x>0)