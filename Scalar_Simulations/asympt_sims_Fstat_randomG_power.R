#################################################
#	Power analysis for GANOVA with scalar outcomes.
#	December 19, 2018
#	bsub < asympt_sims_Fstat_randomG_power.sh
#################################################

library(parallel)
library(vegan)

# setwd('~/Dropbox/Documents/Distance_Statistics/') # laptop
# setwd('/Volumes/RAID Array/Dropbox/Documents/Distance_Statistics/') # desktop
setwd('/project/taki3/Distance_Statistics/') # cluster

#################################################
# FUNCTIONS FOR SCALAR CASES
#################################################

# Distance functions
h.euclid<-function(u) (u[1]-u[2])^2/2
h.abs<-function(u) abs(u[1]-u[2])/2

# These functions calculate the sums of squares.
get.SS<-function(x,h) exp(log(length(x)-1)-log(choose(length(x),2))+log(sum(apply(combn(x,2),2,h))))
get.SSE<-function(x1,x2,h) get.SS(x1,h)+get.SS(x2,h)

# This function estimates the distance between a subject and a list of other subjects.
get.d.vec<-function(x,x.vec) sapply(x.vec,function(u) h(c(x,u)))

get.A.mat<-function(x,d,h,eta,sig) {

	n<-length(x)
	mat.out<-matrix(0,n,n)
	for (jj in 1:(n-1)) mat.out[jj,(jj+1):n]<-get.d.vec(x[jj],x[(jj+1):n])
	dists.mat<-mat.out+t(mat.out)
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

get.ganova.approx<-function(evals,K,N=1E4){
	tmp<-rep(0,N)
	for (j in 1:K) tmp<-tmp+(evals[j]*(rchisq(N,1)-1))
	return(tmp)
}


#This function does the simulation and returns the ganova statistic.
do.simulation<-function(n,p,mu,mu2,sig,h,B.mc=2.5E5,return.details=FALSE,return.sample=TRUE,return.sample.size=5000) {

### SIMULATE DATA
	d<-rbinom(n,1,1-p); if (sum(d)<2) { d[1]<-d[2]<-1 }; if (sum(1-d)<2) { d[1]<-d[2]<-0 }; 
	n2<-sum(d);n1<-n-n2
	x1<-rnorm(n1,mu1,sig)
	x2<-rnorm(n2,mu2,sig)
	x<-rep(NA,n)
	x[d==0]<-x1; x[d==1]<-x2 

### CALCULATE TEST STATISTIC
	start.time<-proc.time()[3]
	SS<-get.SS(c(x1,x2),h)
	SSE<-get.SSE(x1,x2,h)
	numer<-((SS-SSE)/1)
	denom<-(SSE/(length(x1)+length(x2)-2))
	ganova.F<-numer/denom

### CONDUCT ASYMPTOTIC TEST
	A.mat<-get.A.mat(c(x1,x2),c(rep(0,n1),rep(1,n2)),h,n1/n,sqrt(SS/n))
	e.A.mat<-eigen(A.mat)
	#evals<-e.A.mat$values[order(abs(e.A.mat$values),decreasing=TRUE)][1:10]
	evals.all<-e.A.mat$values[order(abs(e.A.mat$values),decreasing=TRUE)]
	n.evals<-sum(cumsum(evals.all)/sum(evals.all)<0.95)
	evals<-e.A.mat$values[order(abs(e.A.mat$values),decreasing=TRUE)][1:n.evals]
	
	#numer.approx<-get.ganova.approx(evals,K=length(evals),N=B.mc)/n+(SS/n)
	#p.val<-mean(ganova.F<numer.approx/denom)
	
	#numer.approx<-get.ganova.approx(evals,K=length(evals),N=B.mc)/n+(SS/n)
	#p.val<-mean(ganova.F<numer.approx/denom)

	Q.n<-numer/denom
	# mean(Q.n<1+get.ganova.approx(evals,K=length(evals),N=B.mc)/SS)
	p.val<-mean(Q.n<1+get.ganova.approx(evals,K=length(evals),N=B.mc)/SS)

	stop.time<-proc.time()[3]

### CONDUCT CLASSICAL F-TEST
	F.p.val<-anova(lm(x~d),lm(x~1),test='F')$Pr[2]

### RETURN RESULTS
	if (return.details==TRUE) {
		return(list(test.time=stop.time-start.time,evals=evals,evals.all=evals.all,ganova.F=ganova.F,p.val=p.val,F.p.val=F.p.val))
	} else {
		if (return.sample==TRUE) {
			numer.approx<-get.ganova.approx(evals,K=length(evals),N=B.mc)/n+(SS/n)
			return(list(test.time=stop.time-start.time,ganova.F=ganova.F,p.val=p.val,F.p.val=F.p.val,asympt.est.sample=numer.approx[1:return.sample.size]/denom))
		} else {
			return(list(test.time=stop.time-start.time,ganova.F=ganova.F,p.val=p.val,F.p.val=F.p.val))
		}
	}
}

#Simulate data under the null
#################################################
# SIMULATE AND COMPUTE FOR SCALAR CASES
#################################################

mu1<-0
mu2<-1
sig<-2
B<-1E4

is.rejected<-function(x,vec) return(x>quantile(unlist(vec),probs=0.95))
n.cores<-as.numeric(Sys.getenv("LSB_DJOB_NUMPROC"))

### DO SIMULATIONS

for (p in c(1/2,1/3)) {

	print(paste('For p=',p,'...'))

	power<-array(NA,dim=c(2,3,4))
	anova.power<-array(NA,dim=c(3,4))
	test.times<-array(NA,dim=c(2,3,4))

	dimnames(power)[[1]]<-c('euclid','abs')
	dimnames(power)[[2]]<-paste0('sig=',c(1,2,3))
	dimnames(power)[[3]]<-paste0('n=',c(30,50,100,500))

	dimnames(anova.power)[[1]]<-paste0('sig=',c(1,2,3))
	dimnames(anova.power)[[2]]<-paste0('n=',c(30,50,100,500))

	dimnames(test.times)[[1]]<-c('euclid','abs')
	dimnames(test.times)[[2]]<-paste0('sig=',c(1,2,3))
	dimnames(test.times)[[3]]<-paste0('n=',c(30,50,100,500))

	i1<-1
	for (sig in c(1,2,3)) {
		i2<-1
		for (n in c(30,50,100,500)) {
			
			set.seed(61343)
			h<-h.euclid;dist.name<-"euclid"
			filename<-paste('ganova_power_2018-12-19_n',n,'_sig',sig,dist.name,'_p',round(p,2),'.RData',sep='')
			existing.files<-list.files(pattern=glob2rx(paste0('ganova_power_2018-12-19_n',n,'_sig',sig,dist.name,'_p',round(p,2),'.RData')))
			if (length(existing.files)>0) {
				load(existing.files[length(existing.files)])
			} else {
				system.time(sims<-mclapply(rep(n,B),do.simulation,p,mu1,mu2,sig,h,mc.cores=n.cores,return.details=TRUE))
				save(sims,file=filename)
			}

			try(power[1,i1,i2]<-mean(unlist(lapply(sims,function(x) x$p.val))<0.05,na.rm=TRUE))
			try(test.times[1,i1,i2]<-mean(unlist(lapply(sims,function(x) x$test.time)),na.rm=TRUE))
			try(anova.power[i1,i2]<-mean(unlist(lapply(sims,function(x) x$F.p.val))<0.05,na.rm=TRUE))

			print(paste('For ',dist.name,', sig=',sig,'n, =',n,' the power is ', power[1,i1,i2]))
			print(paste('For ',dist.name,', sig=',sig,'n, =',n,' the ANOVA power is ', anova.power[i1,i2]))

			print(paste('Finished with n=',n,', sig=',sig,' using Euclidean distance ...'))
			print(paste('Finished simulations with n=',n,', sig=',sig,' using Euclidean distance ...'))

	#######

			h<-h.abs;dist.name<-"abs"
			filename<-paste('ganova_power_2018-12-19_n',n,'_sig',sig,dist.name,'_p',round(p,2),'.RData',sep='')
			existing.files<-list.files(pattern=glob2rx(paste0('ganova_power_2018-12-19_n',n,'_sig',sig,dist.name,'_p',round(p,2),'.RData')))
			if (length(existing.files)>0) {
				load(existing.files[length(existing.files)])
			} else {
				system.time(sims<-mclapply(rep(n,B),do.simulation,p,mu1,mu2,sig,h,mc.cores=n.cores,return.details=TRUE))
				save(sims,file=filename)
			}

			try(test.times[2,i1,i2]<-mean(unlist(lapply(sims,function(x) x$test.time)),na.rm=TRUE))
			try(power[2,i1,i2]<-mean(unlist(lapply(sims,function(x) x$p.val))<0.05,na.rm=TRUE))
			
			print(paste('For ',dist.name,', sig=',sig,'n, =',n,' the power is ', power[2,i1,i2]))
			print(paste('Finished simulations with n=',n,', sig=',sig,' using absolute distance ...'))
			print(paste('Finished with n=',n,', sig=',sig,' using absolute distance ...'))
			i2<-i2+1
		}
		i1<-i1+1
	}

	print(anova.power)
	print(power)
	print(test.times)
	save(test.times,power,anova.power,file=paste0('ganova_2018-12-19_power_sims','_p',round(p,2),'.RData'))
}
