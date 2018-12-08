#################################################
#	Power analysis for GANOVA with functional outcomes.
#	November 28, 2018
#	bsub < asympt_sims_Fstat_function_randomG_power.sh
#################################################

library(ggplot2)
library(R.utils)
library(parallel)
library(vegan)

# setwd('~/Dropbox/Documents/Distance_Statistics/') # laptop
# setwd('/Volumes/RAID Array/Dropbox/Documents/Distance_Statistics/') # desktop
setwd('/project/taki3/Distance_Statistics/') # cluster

resetPar <- function() {
    dev.new()
    op <- par(no.readonly = TRUE)
    dev.off()
    op
}

#################################################
# FUNCTIONS FOR FUNCTIONAL CASES
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

#This function calculates the distance between two graphs as the total number of edges differing between them.
ise<-function(f1,f2,grid=seq(0,1,length.out=100)) sum((f1(grid)-f2(grid))^2)
	
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

	exp.dist<- apply(dists.mat,1,mean) 
	h.1<-dists.mat*(1-(1-d)%o%(1-d)/(eta) - (d)%o%(d)/(1-eta))
	phi.21<-2*eta*((1-d)*exp.dist-eta*sig^2)
	phi.22<-(d-(1-eta))/eta^2
	h.2<-(phi.21%o%phi.22 + phi.22%o%phi.21)/2
	phi.31<-2*(1-eta)*(d*exp.dist-(1-eta)*sig^2) 
	phi.32<--(1-eta)^(-2)*(d-(1-eta))
	h.3<-(phi.31%o%phi.32 + phi.32%o%phi.31)/2
	h.4<-sig^2/((1-eta)*eta)*(d-(1-eta))%o%(d-(1-eta))
	out.mat<-h.1-h.2-h.3-h.4
	return(out.mat)
}

#This function does the graph simulation and returns the simulated data.
simulate.data<-function(n,p,tau1,tau2,sigma,d) {

	dd<-rbinom(n,1,1-p); if (sum(dd)<2) { dd[1]<-dd[2]<-1 }; if (sum(1-dd)<2) { dd[1]<-dd[2]<-0 }; 
	n.2<-sum(dd);n.1<-n-n.2

 	sample.g1<-replicate(n.1,random.f(tau1,sigma),simplify=FALSE)
 	sample.g2<-replicate(n.2,random.f(tau2,sigma),simplify=FALSE)
	
 	return(list(dd=dd,n.1=n.1,n.2=n.2,sample.g1=sample.g1,sample.g2=sample.g2))

}

#This function does the proposed testing.
do.proposed.test<-function(n,d,dd,n.1,n.2,sample.g1,sample.g2,d.mat,B.mc,exact=FALSE,return.details=FALSE,return.sample=TRUE,return.sample.size=5000) {
	
	SS<-mean(d.mat)*n
	SSE<-mean(d.mat[dd%o%dd==1])*n.2+mean(d.mat[(1-dd)%o%(1-dd)==1])*n.1

	numer<-SS-SSE
	denom<-SSE/(n-2)
	ganova.F<-numer/denom

	A.mat<-get.A.mat(d.mat,c(rep(0,n.1),rep(1,n.2)),n.1/n,sqrt(SS/n))
	e.A.mat<-eigen(A.mat)
	#evals<-e.A.mat$values[order(abs(e.A.mat$values),decreasing=TRUE)][1:10]
	evals.all<-e.A.mat$values[order(abs(e.A.mat$values),decreasing=TRUE)]
	n.evals<-sum(cumsum(evals.all)/sum(evals.all)<0.95)
	evals<-e.A.mat$values[order(abs(e.A.mat$values),decreasing=TRUE)][1:n.evals]
	numer.approx<-get.ganova.approx(evals,K=length(evals),N=B.mc)/n+(SS/n)
	p.val<-mean(ganova.F<numer.approx/denom)

### RETURN RESULTS
	if (return.details==TRUE) {
		return(list(ganova.F=ganova.F,evals=evals,evals.all=evals.all,p.val=p.val))
	} else {
		if (return.sample==TRUE) {
			return(list(ganova.F=ganova.F,p.val=p.val,asympt.est.sample=numer.approx[1:return.sample.size]/denom))
		} else {
			return(list(ganova.F=ganova.F,p.val=p.val))
		}
	}

}

#This function does the graph simulation and returns the proposed and standard test.
do.simulation<-function(n,p,tau1=0,tau2=0,sigma,d=ise,B.mc=2.5E5,exact=FALSE) {
	
	simd.data<-simulate.data(n,p,tau1,tau2,sigma,d) 
	group.var<-as.factor(c(rep(0,simd.data$n.1),rep(1,simd.data$n.2)))
	proposed.dmat.time<-system.time(d.mat<-get.d.mat.list(simd.data$sample.g1,simd.data$sample.g2,d=d))[[3]]
	proposed.test.time<-system.time(proposed.test<-do.proposed.test(n=n,d=d,dd=c(rep(0,simd.data$n.1),rep(1,simd.data$n.2)),n.1=simd.data$n.1,n.2=simd.data$n.2,sample.g1=simd.data$sample.g1,sample.g2=simd.data$sample.g2,d.mat=d.mat,B.mc=B.mc,exact=exact,return.details = TRUE))[[3]]+proposed.dmat.time
	competitor.test.time<-system.time(adonis.p<-adonis(as.dist(d.mat)~group.var,permutations=2.5E5)$aov[6][[1]][1])[[3]]+proposed.dmat.time

	return(list(proposed.test.T=proposed.test$ganova.F,proposed.test.p=proposed.test$p.val,asympt.est.sample=proposed.test$asympt.est.sample,adonis.p=adonis.p,proposed.dmat.time=proposed.dmat.time,proposed.test.time=proposed.test.time,competitor.test.time=competitor.test.time,evals=proposed.test$evals,evals.all=proposed.test$evals.all))

}

get.ganova.approx<-function(evals,K,N=1E4){
	tmp<-rep(0,N)
	for (j in 1:K) tmp<-tmp+(evals[j]*(rchisq(N,1)-1))
	return(tmp)
}

#This function takes a sample from the null and compares the observed data to it.
is.rejected<-function(x,vec) return(x>quantile(unlist(vec),probs=0.95))

# Data-generating function
random.f<-function(tau,sigma){
	Z1<-rnorm(1,0,sigma);Z2<-rnorm(1,0,sigma)
	return(function(t) tau*exp(10*(t-0.5))/(1+exp(10*(t-0.5))) + t*Z1 + Z2)
}

# ### PLOT DATA-GENERATING FUNCTIONS
# library(RColorBrewer)
# t<-seq(0,1,length.out=1E3)
# other.cols<-brewer.pal(4,'Dark2')
# plot(t,t,ylim=c(-1,5),type='n',ylab='f(t)')
# for (j in 1:4) {
# 	for (b in 1:20) {
# 		lines(t,random.f(c(0,0.1,0.3,0.5)[j],3E-1)(t),col=other.cols[j],lwd=1)
# 	}
# }
# legend(0,4,legend=paste('tau =',c(0,0.1,0.3,0.5)),col=other.cols,lwd=3)
 
### PLOT DATA-GENERATING FUNCTIONS
library(RColorBrewer)
pdf('example_fns.pdf')
t<-seq(0,1,length.out=1E3)
other.cols<-brewer.pal(4,'Dark2')
plot(t,t,ylim=c(-1,2),type='n',ylab='f(t)')
for (j in 1:2) {
	for (b in 1:100) {
		lines(t,random.f(c(0,0.5)[j],3E-1)(t),col=other.cols[j],lwd=0.2)
	}
}
for (j in 1:2) {
	lines(t,random.f(c(0,0.5)[j],3E-10)(t),col=other.cols[j],lwd=3)
}
legend(0,2,legend=paste('tau =',c(0,0.5)),col=other.cols,lwd=3)
dev.off() 


# ## NOW, compare the distribution of the F-statistic under the null and alternative
B<-1000
V<-LETTERS[1:5]
tau1<-0.05
n.cores<-as.numeric(Sys.getenv("LSB_DJOB_NUMPROC"))
sigma<-3E-1
tau1<-0

for (p in c(1/2,1/3)) {

	print(paste('For p=',p,'...'))

	power<-array(NA,dim=c(3,4))
	adonis.power<-array(NA,dim=c(3,4))
	time.proposed<-array(NA,dim=c(3,4))
	time.competitor<-array(NA,dim=c(3,4))

	dimnames(power)[[1]]<-paste0('tau=2',c(0.1,0.3,0.5))
	dimnames(power)[[2]]<-paste0('n=',c(30,50,100,250))

	dimnames(adonis.power)[[1]]<-paste0('tau2=',c(0.1,0.3,0.5))
	dimnames(adonis.power)[[2]]<-paste0('n=',c(30,50,100,250))

	dimnames(time.proposed)[[1]]<-paste0('tau2=',c(0.1,0.3,0.5))
	dimnames(time.proposed)[[2]]<-paste0('n=',c(30,50,100,250))

	dimnames(time.competitor)[[1]]<-paste0('tau2=',c(0.1,0.3,0.5))
	dimnames(time.competitor)[[2]]<-paste0('n=',c(30,50,100,250))

	i1<-1
	for (tau2 in c(0.1,0.3,0.5)) {
	### ASSESS SIMULATIONS
		i2<-1
		for (n in c(30,50,100,250)) {

			set.seed(976741)
			### DO SIMULATIONS

			### DO SIMULATIONS
			existing.files<-list.files(pattern=glob2rx(paste0('ganova_power_function_2018-11-25_n',n,'_sigma',sigma,'_tau1',tau1,'_tau2',tau2,'_p',round(p,2),'.RData')))
			if (length(existing.files)==0) {	
				system.time(sims<-mclapply(rep(n,B),do.simulation,p,tau1,tau2,sigma,ise,mc.cores=n.cores))
				save(sims,file=paste('ganova_power_function_2018-11-25_n',n,'_sigma',sigma,'_tau1',tau1,'_tau2',tau2,'_p',round(p,2),'.RData',sep=''))
			} else {
				load(existing.files[length(existing.files)])
			}

			try(power[i1,i2]<-mean(unlist(lapply(sims,function(x) x$proposed.test.p))<0.05,na.rm=TRUE))
			try(adonis.power[i1,i2]<-mean(unlist(lapply(sims,function(x) x$adonis.p))<0.05,na.rm=TRUE))
			try(time.proposed[i1,i2]<-mean(unlist(lapply(sims,function(x) x$proposed.test.time)),na.rm=TRUE))
			try(time.competitor[i1,i2]<-mean(unlist(lapply(sims,function(x) x$competitor.test.time)),na.rm=TRUE))
			

			print(paste('For n=',n,' the power is ', power[i1,i2] ,'...'))
			print(paste('For n=',n,' the adonis power is ', adonis.power[i1,i2] ,'...'))

			print(paste('For n=',n,' the average computation time is ', time.proposed[i1,i2] ,'...'))
			print(paste('For n=',n,' the average adonis computation time is ', time.competitor[i1,i2] ,'...'))

			print(paste('Finished with n=',n,'...'))
			i2<-i2+1
		}
		i1<-i1+1
	}
	save(power,adonis.power,time.proposed,time.competitor,file=paste0('ganova_power_2018-11-28_power_sims_function_p',round(p,2),'_sigma',sigma,'.RData'))
}

