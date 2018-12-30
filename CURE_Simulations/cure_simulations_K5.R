#################################################
#	CURE Simulations
#	December 29, 2018
#	bsub < cure_simulations_K5.sh
#################################################

set.seed(821349)
library(vegan)
library(ggplot2)
library(parallel)
library(psych)
eps<-1E-10

#################################################
#################################################

setwd('/project/taki3/Distance_Statistics/') # cluster
data.path<-('/project/taki2/CURE_Data/') # cluster

#################################################
#################################################
#READ IN DATA
#################################################
#################################################

if (!file.exists('CURE_DTI_2018-12-19.RData')) {

	dti.dem<-read.csv(paste(data.path,'DTI/Demographics.csv',sep=''))
	labels<-read.csv(paste(data.path,'DTI/301_labels_coordinates_and_names.csv',sep=''))

	#SORT OUT RIDs (DTI)
	dti.dem$RID<-as.numeric(substr(levels(dti.dem$ID)[dti.dem$ID],3,6))

	#READ IN DATA - use original DTI data (301x301)
	dti.data<-list() ; i<-1
	dti.ids<-c()
	for (id in dti.dem$RID) {
		dti.data[[i]]<-read.csv(paste(data.path,'DTI/DTI_Connectivity/R',formatC(id, width = 4, format = "d", flag = "0"),'_DTI.csv',sep=''),header=FALSE)
		dti.ids<-c(dti.ids,id)
		i<-i+1
	}

	### GET GROUP LABELS
	levels(dti.dem$DX)<-c('ASD','TDC')
	dti.group<-rep(NA,length(dti.ids)); dti.group<-factor(dti.group); levels(dti.group)<-c('ASD','TDC')
	dti.age<-rep(NA,length(dti.ids))
	for (i in 1:length(dti.ids)) {
		dti.group[i]<-dti.dem$DX[which(dti.dem$RID==dti.ids[i])]
		dti.age[i]<-dti.dem$Age[which(dti.dem$RID==dti.ids[i])]
	}

	## DROP SUBJECTS MISSING AGE

	dti.data<-dti.data[!is.na(dti.age)]
	dti.group<-dti.group[!is.na(dti.age)]
	dti.age<-dti.age[!is.na(dti.age)]
	med.age<-median(dti.age)

	n<-length(dti.group);n.1<-sum(dti.group=='TDC');n.2<-sum(dti.group=='ASD')
	dd<-c(rep(0,n.1),rep(1,n.2))

	save.image('CURE_DTI_2018-12-19.RData')

} else {
	load('CURE_DTI_2018-12-19.RData')
}

#################################################
#################################################
#CREATE ALTERNATIVES
#################################################
#################################################

get.upper.tri<-function(mat) mat[upper.tri(mat)]

list.vecs<-lapply(dti.data,get.upper.tri)
dti.mat<-matrix(unlist(list.vecs),ncol=length(list.vecs))

if (!file.exists('CURE_DTI_2018-12-19_alternatives.RData')) {
	regress.on.col<-function(y,x) summary(lm(y~x))$coef[2,4]
	p.vals<-apply(dti.mat,1,regress.on.col,x=dti.group)

	which.edges.alternatives<-list(which(order(p.vals)<=20),which(order(p.vals)<=100),which(order(p.vals)<=200))
	edges.alternatives<-list((order(p.vals)<=20),(order(p.vals)<=100),(order(p.vals)<=200))
	save.image('CURE_DTI_2018-12-19_alternatives.RData')
} else {
	load('CURE_DTI_2018-12-19_alternatives.RData')
}

#################################################
#################################################

# This function estimates the distance between a subject and a list of other subjects.
get.d.list<-function(x,x.list,d) lapply(x.list,d,x)

# Cuts the time in half.
get.d.mat.list<-function(x.list,d) {
	n.list<-length(x.list)
	mat.out<-matrix(0,n.list,n.list)
	for (jj in 1:(n.list-1)) {
		mat.out[jj,(jj+1):n.list]<-unlist(get.d.list(x.list[[jj]],x.list[(jj+1):n.list],d))
	}
	mat.out<-mat.out+t(mat.out)
	return(mat.out)
} 

get.A.mat<-function(dists.mat,d,sig) {
	exp.dist<- apply(dists.mat,1,mean) 
	  pi.s.hat<-table(d)/length(d)
	  d.levels<-sort(unique(d))
	  summand<-matrix(1,length(d),length(d))
	  for (ii in 1:length(d.levels)) {
	    summand<-summand-(d==d.levels[ii])%o%(d==d.levels[ii])/pi.s.hat[ii]
	  }
  	h.1<-dists.mat*summand
  	summand<-matrix(0,length(d),length(d))
  	for (ii in 1:length(d.levels)) {
  	    phi.21<-2*pi.s.hat[ii]*((d==d.levels[ii])*exp.dist-pi.s.hat[ii]*sig^2)
  	    phi.22<-((d!=d.levels[ii])-(1-pi.s.hat[ii]))/pi.s.hat[ii]^2
  	    summand<-summand+(phi.21%o%phi.22 + phi.22%o%phi.21)/2
  	}
  	h.2<-summand
  	summand<-matrix(0,length(d),length(d))
  	for (ii in 1:length(d.levels)) {
  	  summand<-summand+sig^2/pi.s.hat[ii]*(1*(d==d.levels[ii])-pi.s.hat[ii])%o%(1*(d==d.levels[ii])-pi.s.hat[ii])
  	}
  	h.3<-summand
  	out.mat<-h.1-h.2-h.3
  	return(out.mat)
}

get.ganova.approx<-function(evals,K,N=1E4){
	tmp<-rep(0,N)
	for (j in 1:K) tmp<-tmp+(evals[j]*(rchisq(N,1)-1))
	return(tmp)
}

simulate.data<-function(n,mu,k) {

### SAMPLE GROUP SIZES
	dd<-sample(1:5,n,replace=TRUE); if (length(unique(dd))==1) { dd[1]<-1; dd[2]<-2 }
	
### SAMPLE DATA
	simulated.data.mat<-dti.mat[,sample(1:(dim(dti.mat)[2]),n,replace=TRUE)]+mu*(t((dd-1)%x%t(edges.alternatives[[k]])))

	return(list(data.mat=simulated.data.mat,dd=dd))
}

#This function does the proposed testing.
do.proposed.test<-function(n,d,ddd,d.mat,B.mc=2.5E5,exact=FALSE,return.details=FALSE,return.sample=FALSE,return.sample.size=5000) {

	SS<-mean(d.mat)*n
	SSE<-0; ddd.levels<-sort(unique(ddd))
	for (ii in 1:length(ddd.levels)) SSE<-SSE+mean(d.mat[(ddd==ddd.levels[ii])%o%(ddd==ddd.levels[ii])==1])*sum(1*ddd==1*ddd.levels[ii]) 
	
	numer<-(SS-SSE)/(length(ddd.levels)-1)
	denom<-SSE/(n-length(ddd.levels))
	ganova.F<-numer/denom

	A.mat<-get.A.mat(d.mat,ddd,sqrt(SS/n))
	e.A.mat<-eigen(A.mat)
	evals.all<-e.A.mat$values[order(abs(e.A.mat$values),decreasing=TRUE)]
	n.evals<-sum(cumsum(evals.all)/sum(evals.all)<0.95)
	evals<-e.A.mat$values[order(abs(e.A.mat$values),decreasing=TRUE)][1:n.evals]
	
	ganova.F<-numer/denom

	Q.n<-numer/denom
	Q.approx<-1+get.ganova.approx(evals,K=length(evals),N=B.mc)/SS/(length(ddd.levels)-1)
	p.val<-mean(Q.n<Q.approx)

### RETURN RESULTS
	if (return.details==TRUE) {
		return(list(ganova.F=ganova.F,evals=evals,evals.all=evals.all,p.val=p.val))
	} else {
		if (return.sample==TRUE) {
			return(list(ganova.F=ganova.F,p.val=p.val,asympt.est.sample=Q.approx[1:return.sample.size]))
		} else {
			return(list(ganova.F=ganova.F,p.val=p.val))
		}
	}

}

mat.to.list<-function(mat) lapply(seq_len(ncol(mat)), function(i) mat[,i])
vec.d<-function(x1,x2) sum(x1-x2)^2

#This function does the graph simulation and returns the proposed and standard test.
do.simulation<-function(n,mu,k,vec.d,B.mc=2.5E5,exact=FALSE) {
	
	simd.data<-simulate.data(n,mu,k) 

	proposed.dmat.time<-system.time(d.mat<-get.d.mat.list(mat.to.list(simd.data$data.mat),d=vec.d))[[3]]
	proposed.test.time<-system.time(proposed.test<-do.proposed.test(n=n,d=vec.d,dd=simd.data$dd,d.mat=d.mat,B.mc=B.mc,exact=exact,return.details = TRUE))[[3]]+proposed.dmat.time
	competitor.test.time<-system.time(adonis.p<-adonis(as.dist(d.mat)~factor(simd.data$dd),permutations=2.5E5)$aov[6][[1]][1])[[3]]+proposed.dmat.time

	return(list(proposed.test.T=proposed.test$ganova.F,proposed.test.p=proposed.test$p.val,proposed.test.p.old=proposed.test$p.val.old,asympt.est.sample=proposed.test$asympt.est.sample,adonis.p=adonis.p,proposed.dmat.time=proposed.dmat.time,proposed.test.time=proposed.test.time,competitor.test.time=competitor.test.time,evals=proposed.test$evals,evals.all=proposed.test$evals.all))

}


#################################################
#################################################
#DO SIMULATIONS
#################################################
#################################################

B<-1E3
n.cores<-as.numeric(Sys.getenv("LSB_DJOB_NUMPROC"))
job.index<-as.numeric(Sys.getenv("LSB_JOBINDEX"))
p<-1/2

mu.vals<-c(0,7.5,15,30)
ks<-c(20,100,200)
k.vals<-1:3
n.vals<-c(30,50,100,264)
par.matrix<-expand.grid(mu.vals,k.vals,n.vals)
mu<-par.matrix[job.index,1]
k<-par.matrix[job.index,2]
n<-par.matrix[job.index,3]

### DO SIMULATIONS
print(paste('Running simulations for n=',n,', mu=',mu,', k=',ks[k],'...'))
existing.files<-list.files(pattern=glob2rx(paste0('ganova_simulations_dti_K5_2018-12-29_test_n',n,'_mu',mu,'_k',ks[k],'.RData')))
if (length(existing.files)==0) {	
	system.time(sims<-mclapply(rep(n,B),do.simulation,mu,k,vec.d,mc.cores=n.cores))
	save(sims,file=paste('ganova_simulations_dti_K5_2018-12-29_test_n',n,'_mu',mu,'_k',ks[k],'.RData',sep=''))
} else {
	print('Simulations previously run!')
	load(existing.files[length(existing.files)])
}

try(rejection.rate.proposed<-mean(unlist(lapply(sims,function(x) x$proposed.test.p))<0.05,na.rm=TRUE))
try(rejection.rate.competitor<-mean(unlist(lapply(sims,function(x) x$adonis.p))<0.05,na.rm=TRUE))
try(time.proposed<-mean(unlist(lapply(sims,function(x) x$proposed.test.time)),na.rm=TRUE))
try(time.competitor<-mean(unlist(lapply(sims,function(x) x$competitor.test.time)),na.rm=TRUE))

if (mu==0) {
	print(paste('The type 1 error rate is ', rejection.rate.proposed ,'...'))
	print(paste('The adonis type 1 error rate is ', rejection.rate.competitor ,'...'))

} else {
	print(paste('The power is ', rejection.rate.proposed ,'...'))
	print(paste('The adonis power is ', rejection.rate.competitor ,'...'))
}

print(paste('The average computation time is ', time.proposed,'...'))
print(paste('The average adonis computation time is ', time.competitor,'...'))
