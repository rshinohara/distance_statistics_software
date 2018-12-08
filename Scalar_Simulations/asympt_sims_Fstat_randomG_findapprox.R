#################################################
#	Assess the necessary number of samples/permutations for our methods to be able to compare computation times.
#	October 1, 2017
#	bsub < asympt_sims_Fstat_randomG.sh
#################################################

library(parallel)
library(RColorBrewer)

# setwd('~/Dropbox/Documents/Distance_Statistics/') # laptop
# setwd('/Volumes/RAID Array/Dropbox/Documents/Distance_Statistics/') # desktop
setwd('~/Distance_Statistics/') # cluster

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
	# xmat<-cbind(rep(x,each=length(x)),rep(x,length.out=length(x)^2))
	# dists<-apply(xmat,1,h)
	# dists.mat<-matrix(dists,nrow=length(x))
	
	exp.dist<- apply(dists.mat,1,mean) # This is the expected distance from each subject.
	A.term1<-(1-1/eta*(d%o%d)-1/(1-eta)*((1-d)%o%(1-d)))*dists.mat
	A.term2<-(1-(1-d)/eta)%o%((1-d)*exp.dist-sig^2*eta);A.term2<-A.term2+t(A.term2)
	A.term3<-(1-d/(1-eta))%o%(d*exp.dist-sig^2*(1-eta));A.term3<-A.term3+t(A.term3)
	A.term4<-(1-d-eta)%o%(1-d-eta)*sig^2*eta^(-1)*(1-eta)^(-1)
	return(A.term1-A.term2-A.term3-A.term4)
}

get.ganova.approx<-function(evals,K,N=1E4){
	tmp<-rep(0,N)
	for (j in 1:K) tmp<-tmp+(evals[j]*(rchisq(N,1)-1))
	return(tmp)
}

#This function does the simulation and returns the ganova statistic.
do.simulation<-function(n,p,mu,sig,h,B.mc.grid=c(1E3,2E3,5E3,1E4,2E4,5E4),B.mc.max=1E5,return.details=FALSE,return.sample=TRUE,return.sample.size=5000) {

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
	A.mat<-get.A.mat(c(x1,x2),c(rep(0,n1),rep(1,n2)),h,n1/n,sqrt(SS/n))
	e.A.mat<-eigen(A.mat)
	evals<-e.A.mat$values[order(abs(e.A.mat$values),decreasing=TRUE)][1:10]
	mid.time<-proc.time()[3]
	p.vals.grid<-rep(NA,length(B.mc.grid)+1)
	approx.times<-rep(NA,length(B.mc.grid)+1)
	for (jj in 1:(length(B.mc.grid)+1)) {
		approx.time.start<-proc.time()[3]
		numer.approx<-get.ganova.approx(evals,K=10,N=c(B.mc.grid,B.mc.max)[jj])/n+(SS/n)
		p.vals.grid[jj]<-mean(ganova.F<numer.approx/denom)
		approx.time.stop<-proc.time()[3]
		approx.times[jj]<-approx.time.stop-approx.time.start
	}
	test.times<-mid.time-start.time+approx.times
	p.vals.errs<-(p.vals.grid-p.vals.grid[length(p.vals.grid)])^2

### DO CLASSICAL F-TEST
	F.p.val<-anova(lm(x~d),lm(x~1),test='F')$Pr[2]

### RETURN RESULTS
	if (return.details==TRUE) {
		return(list(test.times=test.times,SS=SS,SSE=SSE,numer=numer,denom=denom,ganova.F=ganova.F,evals=evals,numer.approx=numer.approx,p.vals.grid=p.vals.grid,p.vals.errs=p.vals.errs,F.p.val=F.p.val))
	} else {
		if (return.sample==TRUE) {
			return(list(test.times=test.times,ganova.F=ganova.F,p.vals.grid=p.vals.grid,p.vals.errs=p.vals.errs,F.p.val=F.p.val,asympt.est.sample=numer.approx[1:return.sample.size]/denom))
		} else {
			return(list(test.times=test.times,ganova.F=ganova.F,p.vals.grid=p.vals.grid,p.vals.errs=p.vals.errs,F.p.val=F.p.val))
		}
	}
}

is.rejected<-function(x,vec) return(x>quantile(unlist(vec),probs=0.95))

pos.part<-function(x) x*(x>0)

#################################################
# SIMULATE AND COMPUTE FOR SCALAR CASES
#################################################

#Simulate data under the null
p<-0.5; B<-1E4; mu<-0

F.type1err<-array(NA,dim=c(2,5))
type1err<-array(NA,dim=c(2,5,7))
test.times<-array(NA,dim=c(2,5,7))
pval.err.rates<-array(NA,dim=c(2,5,7))

### JUST FOR TESTING
# n<-5E2;p<-0.5;B<-100
# p<-0.5; mu<-0; sig<-2; B<-1000
# h<-h.euclid;dist.name<-"euclid"
# h<-h.abs;dist.name<-"abs"

### DO SIMULATIONS
set.seed(61343); sig<-1
cols<-brewer.pal(5,'Pastel1')
other.cols<-brewer.pal(9,'YlGnBu')
i1<-1 
for (dist.name in c('euclid','abs')) {
	i2<-1
	for (n in c(30,50,100,250,500)) {

	### RUN SIMULATIONS
		if (dist.name=='euclid') h<-h.euclid
		if (dist.name=='abs') h<-h.abs
		system.time(sims<-mclapply(rep(n,B),do.simulation,p,mu,sig,h,mc.cores=20))
		save(sims,file=paste('ganova_type1err_times_2017-10-01_dist_',dist.name,'_n',n,'.RData',sep=''))

	### MAKE PLOTS FOR COMPARISON TO EMPIRICAL NULL DISTRIBUTION
		jpeg(paste0('ganova_2017-10-01_qq_dist_',dist.name,'_n',n,'.jpg'))
		qqplot(unlist(lapply(sims,function(x) x$ganova.F)),pos.part(lapply(sims,function(x) x$asympt.est.sample)[[1]]),type='l',col=sample(other.cols),ylab='Asymptotic Approximation',xlab='Empirical Null Distribution',xlim=c(0,10),ylim=c(0,10))
		for (b in 2:B) {
			qq.pts<-qqplot(unlist(lapply(sims,function(x) x$ganova.F)),pos.part(lapply(sims,function(x) x$asympt.est.sample)[[b]]),type='l',xlab='Asymptotic Approximation',ylab='Empirical Null Distribution',plot.it=FALSE)
			lines(qq.pts,col=sample(other.cols))
		}
		abline(a=0,b=1,col='#fec44f')
		dev.off()

	### ESTIMATE ERROR RATES
		pval.err.rates[i1,i2,]<-apply(matrix(unlist(lapply(sims,function(x) x$p.vals.errs)),nrow=B,byrow=TRUE),2,mean)

	### RECORD TYPE I ERROR RATES
		type1err[i1,i2,]<-apply(matrix(unlist(lapply(sims,function(x) x$p.vals.grid)),nrow=B,byrow=TRUE)<0.05,2,mean,na.rm=TRUE)
		F.type1err[i1,i2]<-mean(unlist(lapply(sims,function(x) x$F.p.val))<0.05,na.rm=TRUE)
		test.times[i1,i2,]<-apply(matrix(unlist(lapply(sims,function(x) x$test.times)),nrow=B,byrow=TRUE),2,mean)
		
		print(paste('Finished with n=',n,'...'))

		i2<-i2+1
	}
	i1<-i1+1
}
save(type1err,F.type1err,test.times,file='ganova_2017-10-01_type1err_times_sims.RData')

pdf('test.pdf')
plot(test.times[1,1,],pval.err.rates[1,1,],type='n',xlim=c(min(test.times),max(test.times)),ylim=c(min(pval.err.rates),max(pval.err.rates)))
for (j in 1:5)
	lines(test.times[1,j,],pval.err.rates[1,j,],type='l',col=cols[j],lwd=3)
legend(max(test.times)-(max(test.times)-min(test.times))/3,max(pval.err.rates)-(max(pval.err.rates)-min(pval.err.rates))/3,paste0('n=',c(30,50,100,250,500)),col=cols,lwd=3)
dev.off()


pdf('approx_time_scalar_n250.pdf')
eps.x<-(max(test.times[1,4,])-min(test.times[1,4,]))/30
eps.y<-(max(pval.err.rates[1,4,])-min(pval.err.rates[1,4,]))/30
plot(test.times[1,4,],pval.err.rates[1,4,],type='l',xlim=c(min(test.times[1,4,]),max(test.times[1,4,])+eps.x),ylim=c(min(pval.err.rates[1,4,]),max(pval.err.rates[1,4,])+eps.y),lwd=3,col=cols[2],xlab='Average Test Time (s)',ylab='p-value MSE')
points(test.times[1,4,],pval.err.rates[1,4,],pch=2)
eps.x<-(max(test.times[1,4,])-min(test.times[1,4,]))/30
eps.y<-(max(pval.err.rates[1,4,])-min(pval.err.rates[1,4,]))/30
text(x=test.times[1,4,]+eps.x,y=pval.err.rates[1,4,]+eps.y,labels=c('1E3','2E3','5E3','1E4','2E4','5E4','1E5'))
dev.off()

pdf('approx_time_scalar_abs_n250.pdf')
eps.x<-(max(test.times[2,4,])-min(test.times[2,4,]))/30
eps.y<-(max(pval.err.rates[2,4,])-min(pval.err.rates[2,4,]))/30
plot(test.times[2,4,],pval.err.rates[2,4,],type='l',xlim=c(min(test.times[2,4,]),max(test.times[2,4,])+eps.x),ylim=c(min(pval.err.rates[2,4,]),max(pval.err.rates[2,4,])+eps.y),lwd=3,col=cols[2],xlab='Average Test Time (s)',ylab='p-value MSE')
points(test.times[2,4,],pval.err.rates[2,4,],pch=2)
eps.x<-(max(test.times[2,4,])-min(test.times[2,4,]))/30
eps.y<-(max(pval.err.rates[2,4,])-min(pval.err.rates[2,4,]))/30
text(x=test.times[2,4,]+eps.x,y=pval.err.rates[2,4,]+eps.y,labels=c('1E3','2E3','5E3','1E4','2E4','5E4','1E5'))
dev.off()


pdf('test.pdf')
plot(test.times[1,1,],type1err[1,1,],type='n',xlim=c(min(test.times),max(test.times)),ylim=c(min(type1err),max(type1err)))
for (j in 1:5)
	lines(test.times[1,j,],type1err[1,j,],type='l',col=cols[j],lwd=3)
legend(max(test.times)-(max(test.times)-min(test.times))/3,max(type1err)-(max(type1err)-min(type1err))/3,paste0('n=',c(30,50,100,250,500)),col=cols,lwd=3)
dev.off()



print(type1err)
print(F.type1err)
print(test.times)
save(type1err,F.type1err,test.times,file='ganova_2017-10-01_type1err_times_sims.RData')
