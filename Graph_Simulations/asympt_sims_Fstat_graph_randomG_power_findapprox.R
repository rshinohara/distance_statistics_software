#################################################
#	Confirmatory analyses for GANOVA aymptotics - marginal distribution of U-statistics
#	October 18, 2017
#	bsub < asympt_sims_Fstat_graph_randomG_power_findapprox.sh
#################################################

library(ggplot2)
library(graph)
library(gRbase)
library(R.utils)
library(parallel)
library(vegan)
library(RColorBrewer)

cols<-brewer.pal(3,'Dark2')

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


#This function calculates the distance between two graphs as the total number of edges differing between them.
edge.diff<-function(g1,g2) {

	edges.1<-edgeList(g1,matrix=TRUE)
	edges.2<-edgeList(g2,matrix=TRUE)
	edges.1.string<-apply(edges.1,1,paste,collapse='_')
	edges.2.string<-apply(edges.2,1,paste,collapse='_')

	return(length(setdiff(edges.1.string, edges.2.string))+length(setdiff(edges.2.string, edges.1.string)))
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
do.proposed.test<-function(n,d,dd,n.1,n.2,sample.g1,sample.g2,d.mat,exact=FALSE,return.details=FALSE,return.sample=TRUE,return.sample.size=5000) {

	p.vals.grid<-rep(NA,length(mc.grid))
	test.times<-rep(NA,length(mc.grid))
	
	base.time.start<-proc.time()[3]
	SS<-mean(d.mat)*n
	SSE<-mean(d.mat[dd%o%dd==1])*n.2+mean(d.mat[(1-dd)%o%(1-dd)==1])*n.1
	numer<-SS-SSE
	denom<-SSE/(n-2)
	ganova.F<-numer/denom

### CONDUCT ASYMPTOTIC TEST
	A.mat<-get.A.mat(d.mat,c(rep(0,n.1),rep(1,n.2)),n.1/n,sqrt(SS/n))
	e.A.mat<-eigen(A.mat)
	evals<-e.A.mat$values[order(abs(e.A.mat$values),decreasing=TRUE)][1:10]
	
	base.time.stop<-proc.time()[3]

	for (jj in 1:(length(mc.grid))) {
		approx.time.start<-proc.time()[3]
		numer.approx<-get.ganova.approx(evals,K=10,N=mc.grid[jj])/n+(SS/n)
		p.vals.grid[jj]<-mean(ganova.F<numer.approx/denom)
		approx.time.stop<-proc.time()[3]
		test.times[jj]<-base.time.stop-base.time.start+approx.time.stop-approx.time.start
	}
	p.vals.errs<-(p.vals.grid-p.vals.grid[length(p.vals.grid)])^2

### RETURN RESULTS
	if (return.details==TRUE) {
		return(list(test.times=test.times,SS=SS,SSE=SSE,numer=numer,denom=denom,ganova.F=ganova.F,evals=evals,numer.approx=numer.approx,p.vals.grid=p.vals.grid,p.vals.errs=p.vals.errs))
	} else {
		if (return.sample==TRUE) {
			return(list(test.times=test.times,ganova.F=ganova.F,p.vals.grid=p.vals.grid,p.vals.errs=p.vals.errs,asympt.est.sample=numer.approx[1:return.sample.size]/denom))
		} else {
			return(list(test.times=test.times,ganova.F=ganova.F,p.vals.grid=p.vals.grid,p.vals.errs=p.vals.errs))
		}
	}

}


#This function does the graph simulation and returns the ganova statistic.
do.simulation<-function(n,p,tau1,tau2,d,exact=FALSE) {

	dd<-rbinom(n,1,1-p)
	n.2<-sum(dd);n.1<-n-n.2

 	pmat<-function(tau) cbind(c(1,1,1-tau,tau, tau),c(1,1, tau, tau, tau),c(1-tau, tau,1,1,1),c(tau, tau,1,1,1),c(tau, tau,1,1,1))
 	sample.g1<-replicate(n.1,random.graph(V,pmat(tau1)),simplify=FALSE)
 	sample.g2<-replicate(n.2,random.graph(V,pmat(tau2)),simplify=FALSE)
	
	group.var<-as.factor(c(rep(0,n.1),rep(1,n.2)))
	proposed.dmat.time<-system.time(d.mat<-get.d.mat.list(sample.g1,sample.g2,d=d))[[3]]
	proposed.test<-do.proposed.test(n=n,d=d,dd=dd,n.1=n.1,n.2=n.2,sample.g1=sample.g1,sample.g2=sample.g2,d.mat=d.mat,exact=exact)
	proposed.test.times<-proposed.test$test.times+proposed.dmat.time
	
	adonis.p<-rep(NA,length(perm.grid))
	adonis.times<-rep(NA,length(perm.grid))
	for (j in 1:length(perm.grid)) adonis.times[j]<-system.time(adonis.p[j]<-adonis(as.dist(d.mat)~group.var,permutations=perm.grid[j])$aov[6][[1]][1])[[3]]+proposed.dmat.time
	adonis.p.errs<-(adonis.p-adonis.p[length(perm.grid)])^2

	return(list(proposed.test.T=proposed.test$ganova.F,proposed.test.p=proposed.test$p.vals.grid,proposed.test.p.err=proposed.test$p.vals.errs,asympt.est.sample=proposed.test$asympt.est.sample,adonis.p=adonis.p,adonis.p.errs=adonis.p.errs,adonis.times=adonis.times,proposed.test.times=proposed.test.times))

}

get.ganova.approx<-function(evals,K,N=1E4){
	tmp<-rep(0,N)
	for (j in 1:K) tmp<-tmp+(evals[j]*(rchisq(N,1)-1))
	return(tmp)
}

#This function takes a sample from the null and compares the observed data to it.
is.rejected<-function(x,vec) return(x>quantile(unlist(vec),probs=0.95))

pos.part<-function(x) x*(x>0)

## NOW, compare the distribution of the F-statistic under the null and alternative
p<-0.5; B<-5E2
V<-LETTERS[1:5]
set.seed(976741)
n<-250

perm.grid<-c(2.5E4,5E4,7.5E4,1E5,2.5E5)
mc.grid<-c(1E5,2.5E5,5E5,7.5E5,1E6)
n.cores<-as.numeric(Sys.getenv("LSB_DJOB_NUMPROC"))

### DO SIMULATIONS
system.time(sims<-mclapply(rep(n,B),do.simulation,p,0.05,0.1,edge.diff,mc.cores=n.cores))

### ASSESS ERROR AND COMPUTING TIMES
power.proposed<-apply(matrix(unlist(lapply(sims,function(x) x$proposed.test.p)),nrow=B,byrow=TRUE)<0.05,2,mean,na.rm=TRUE)
power.competitor<-apply(matrix(unlist(lapply(sims,function(x) x$adonis.p)),nrow=B,byrow=TRUE)<0.05,2,mean,na.rm=TRUE)
p.err.proposed<-apply(matrix(unlist(lapply(sims,function(x) x$proposed.test.p.err)),nrow=B,byrow=TRUE),2,mean,na.rm=TRUE)
p.err.competitor<-apply(matrix(unlist(lapply(sims,function(x) x$adonis.p.errs)),nrow=B,byrow=TRUE),2,mean,na.rm=TRUE)
time.proposed<-apply(matrix(unlist(lapply(sims,function(x) x$proposed.test.times)),nrow=B,byrow=TRUE),2,mean,na.rm=TRUE)
time.competitor<-apply(matrix(unlist(lapply(sims,function(x) x$adonis.times)),nrow=B,byrow=TRUE),2,mean,na.rm=TRUE)

save(sims,power.proposed,power.proposed,p.err.proposed,p.err.competitor,time.proposed,time.competitor,file='ganova_graph_approx_times_smalleffect_2017-10-18.RData')
# load('ganova_graph_2017-10-05_approx_times_smalleffect_2017-10-06.RData')


pdf('approx_time_graph_n250_smalleffect_2017-10-18.pdf',height=7,width=14)
par(mfrow=c(1,2))
x.lims<-c(0.99*min(c(time.proposed)),1.01*max(c(time.proposed)))
y.lims<-c(min(c(p.err.competitor,p.err.proposed)),max(c(p.err.competitor,p.err.proposed)))
eps.x<-(x.lims[2]-x.lims[1])/30;eps.y<-(y.lims[2]-y.lims[1])/30
x.lims[2]<-x.lims[2]+eps.x;y.lims[2]<-y.lims[2]+eps.y
plot(time.proposed,p.err.proposed,type='l',xlim=x.lims,ylim=y.lims,lwd=3,col=cols[2],xlab='Average Test Time (s)',ylab='p-value MSE',main='Proposed Test')
points(time.proposed,p.err.proposed,pch=2)
text(x=time.proposed+eps.x,y=p.err.proposed+eps.y,labels=mc.grid)
x.lims<-c(0.5*min(c(time.competitor)),1.1*max(c(time.competitor)))
plot(time.competitor,p.err.competitor,type='l',xlim=x.lims,ylim=y.lims,lwd=3,col=cols[2],xlab='Average Test Time (s)',ylab='p-value MSE',main='Permutation Test')
points(time.competitor,p.err.competitor,pch=2)
text(x=time.competitor+eps.x,y=p.err.competitor+eps.y,labels=perm.grid)
dev.off()

set.seed(843)

### DO SIMULATIONS
system.time(sims<-mclapply(rep(n,B),do.simulation,p,0.05,0.2,edge.diff,mc.cores=50))

### ASSESS ERROR AND COMPUTING TIMES
power.proposed<-apply(matrix(unlist(lapply(sims,function(x) x$proposed.test.p)),nrow=B,byrow=TRUE)<0.05,2,mean,na.rm=TRUE)
power.competitor<-apply(matrix(unlist(lapply(sims,function(x) x$adonis.p)),nrow=B,byrow=TRUE)<0.05,2,mean,na.rm=TRUE)
p.err.proposed<-apply(matrix(unlist(lapply(sims,function(x) x$proposed.test.p.err)),nrow=B,byrow=TRUE),2,mean,na.rm=TRUE)
p.err.competitor<-apply(matrix(unlist(lapply(sims,function(x) x$adonis.p.errs)),nrow=B,byrow=TRUE),2,mean,na.rm=TRUE)
time.proposed<-apply(matrix(unlist(lapply(sims,function(x) x$proposed.test.times)),nrow=B,byrow=TRUE),2,mean,na.rm=TRUE)
time.competitor<-apply(matrix(unlist(lapply(sims,function(x) x$adonis.times)),nrow=B,byrow=TRUE),2,mean,na.rm=TRUE)

save(sims,power.proposed,power.proposed,p.err.proposed,p.err.competitor,time.proposed,time.competitor,file='ganova_graph_approx_times_largeeffect_2017-10-18.RData')

pdf('approx_time_graph_n250_largeeffect_2017-10-18.pdf',height=7,width=14)
par(mfrow=c(1,2))
x.lims<-c(0.99*min(c(time.proposed)),1.01*max(c(time.proposed)))
y.lims<-c(min(c(p.err.competitor,p.err.proposed)),max(c(p.err.competitor,p.err.proposed)))
eps.x<-(x.lims[2]-x.lims[1])/30;eps.y<-(y.lims[2]-y.lims[1])/30
x.lims[2]<-x.lims[2]+eps.x;y.lims[2]<-y.lims[2]+eps.y
plot(time.proposed,p.err.proposed,type='l',xlim=x.lims,ylim=y.lims,lwd=3,col=cols[2],xlab='Average Test Time (s)',ylab='p-value MSE',main='Proposed Test')
points(time.proposed,p.err.proposed,pch=2)
text(x=time.proposed+eps.x,y=p.err.proposed+eps.y,labels=mc.grid)
x.lims<-c(0.5*min(c(time.competitor)),1.1*max(c(time.competitor)))
plot(time.competitor,p.err.competitor,type='l',xlim=x.lims,ylim=y.lims,lwd=3,col=cols[2],xlab='Average Test Time (s)',ylab='p-value MSE',main='Permutation Test')
points(time.competitor,p.err.competitor,pch=2)
text(x=time.competitor+eps.x,y=p.err.competitor+eps.y,labels=perm.grid)
dev.off()
