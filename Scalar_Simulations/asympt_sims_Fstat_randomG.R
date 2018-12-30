#################################################
#	Confirmatory analyses for GANOVA aymptotics - marginal distribution of U-statistics
#	December 19, 2018
#	bsub < asympt_sims_Fstat_randomG.sh
#################################################

library(parallel)
library(RColorBrewer)

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
do.simulation<-function(n,p,mu,sig,h,B.mc=2.5E5,return.details=FALSE,return.sample=TRUE,return.sample.size=5000) {

### SIMULATE DATA
	d<-rbinom(n,1,1-p); if (sum(d)<2) { d[1]<-d[2]<-1 }; if (sum(1-d)<2) { d[1]<-d[2]<-0 }; 
	x<-rnorm(n,mu,sig)
	n2<-sum(d);n1<-n-n2
	x1<-x[d==0]
	x2<-x[d==1]
 
### CALCULATE TEST STATISTIC
	start.time<-proc.time()[3]
	SS<-get.SS(c(x1,x2),h)
	SSE<-get.SSE(x1,x2,h)
	numer<-((SS-SSE)/1)
	denom<-(SSE/(length(x1)+length(x2)-2))
	ganova.F<-numer/denom

### CONDUCT ASYMPTOTIC TEST
	A.mat<-get.A.mat(c(x1,x2),c(rep(0,n1),rep(1,n2)),h,n1/n,sqrt(SS/n)) ## FIXED HERE (\eta = P(d=0))
	e.A.mat<-eigen(A.mat)
	#evals<-e.A.mat$values[order(abs(e.A.mat$values),decreasing=TRUE)][1:10]
	evals.all<-e.A.mat$values[order(abs(e.A.mat$values),decreasing=TRUE)]
	n.evals<-sum(cumsum(evals.all)/sum(evals.all)<0.95)
	evals<-e.A.mat$values[order(abs(e.A.mat$values),decreasing=TRUE)][1:n.evals]

#	numer.approx<-get.ganova.approx(evals,K=length(evals),N=B.mc)/n+(SS/n)
#	p.val<-mean(ganova.F<numer.approx/denom)

	# sum.approx<-get.ganova.approx(evals,K=length(evals),N=B.mc)
	# numer.approx<-SS/n + sum.approx/n
	# p.val<-mean(ganova.F<numer.approx/denom)
	# mean(ganova.F<numer.approx/denom)

	Q.n<-numer/denom
	# mean(Q.n<1+get.ganova.approx(evals,K=length(evals),N=B.mc)/SS)
	p.val<-mean(Q.n<1+get.ganova.approx(evals,K=length(evals),N=B.mc)/SS)

	stop.time<-proc.time()[3]

### DO CLASSICAL F-TEST
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

is.rejected<-function(x,vec) return(x>quantile(unlist(vec),probs=0.95))

pos.part<-function(x) x*(x>0)

#################################################
# SIMULATE AND COMPUTE FOR SCALAR CASES
#################################################

#Simulate data under the null
B<-1E4
mu<-0
sig<-1

### JUST FOR TESTING
# n<-5E2;p<-0.5;B<-100
# p<-0.5; mu<-0; sig<-2; B<-1000
# h<-h.euclid;dist.name<-"euclid"
# h<-h.abs;dist.name<-"abs"

### DO SIMULATIONS
cols<-brewer.pal(3,'Pastel1')
other.cols<-brewer.pal(9,'YlGnBu')
n.cores<-as.numeric(Sys.getenv("LSB_DJOB_NUMPROC"))

for (p in c(1/2,1/3)) {
#p<-1/3

	print(paste('For p=',p,'...'))

	type1err<-array(NA,dim=c(2,4))
	F.type1err<-array(NA,dim=c(2,4))
	test.times<-array(NA,dim=c(2,4))

	dimnames(type1err)[[1]]<-c('euclid','abs')
	dimnames(type1err)[[2]]<-paste0('n=',c(30,50,100,500))

	dimnames(F.type1err)[[1]]<-c('euclid','abs')
	dimnames(F.type1err)[[2]]<-paste0('n=',c(30,50,100,500))

	dimnames(test.times)[[1]]<-c('euclid','abs')
	dimnames(test.times)[[2]]<-paste0('n=',c(30,50,100,500))

	i1<-1 
	for (dist.name in c('euclid','abs')) {
		i2<-1
		for (n in c(30,50,100,500)) {

			set.seed(61343)
		### RUN SIMULATIONS
			if (!file.exists(paste('ganova_type1err_2018-12-19_dist_',dist.name,'_n',n,'_p',round(p,2),'.RData',sep=''))) {
				if (dist.name=='euclid') h<-h.euclid
				if (dist.name=='abs') h<-h.abs
				system.time(sims<-mclapply(rep(n,B),do.simulation,p,mu,sig,h,mc.cores=n.cores,return.details=TRUE))
				save(sims,file=paste('ganova_type1err_2018-12-19_dist_',dist.name,'_n',n,'_p',round(p,2),'.RData',sep=''))
			} else {
				load(paste('ganova_type1err_2018-12-19_dist_',dist.name,'_n',n,'_p',round(p,2),'.RData',sep=''))
			}

		### MAKE PLOTS FOR COMPARISON TO EMPIRICAL NULL DISTRIBUTION
		#	jpeg(paste0('ganova_2017-10-19_qq_dist_',dist.name,'_n',n,'_p',round(p,2),'.jpg'))
		#	qqplot(unlist(lapply(sims,function(x) x$ganova.F)),pos.part(lapply(sims,function(x) x$asympt.est.sample)[[1]]),type='l',col=sample(other.cols),ylab='Asymptotic Approximation',xlab='Empirical Null Distribution',xlim=c(0,10),ylim=c(0,10))
		#	for (b in 2:B) {
		#		qq.pts<-qqplot(unlist(lapply(sims,function(x) x$ganova.F)),pos.part(lapply(sims,function(x) x$asympt.est.sample)[[b]]),type='l',xlab='Asymptotic Approximation',ylab='Empirical Null Distribution',plot.it=FALSE)
		#		lines(qq.pts,col=sample(other.cols))
		#	}
		#	abline(a=0,b=1,col='#fec44f')
		#	dev.off()

		### RECORD TYPE I ERROR RATES
			try(type1err[i1,i2]<-mean(unlist(lapply(sims,function(x) x$p.val))<0.05,na.rm=TRUE))
			try(F.type1err[i1,i2]<-mean(unlist(lapply(sims,function(x) x$F.p.val))<0.05,na.rm=TRUE))
			try(test.times[i1,i2]<-mean(unlist(lapply(sims,function(x) x$test.time)),na.rm=TRUE))
			
			print(paste('Finished with n=',n,'...'))

			i2<-i2+1
		}
		i1<-i1+1
	}

	print(type1err)
	print(F.type1err)
	print(test.times)
	save(type1err,F.type1err,test.times,file=paste0('ganova_2018-12-19_type1err_sims','_p',round(p,2),'.RData'))
}
