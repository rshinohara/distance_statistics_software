#################################################
#	CURE Analysis - Updated for resubmission
#	November 25, 2018
#	bsub < cure_analysis_resubmission.sh
#################################################

set.seed(821349)
library(vegan)
library(ggplot2)
library(parallel)
library(psych)

#################################################
#################################################

setwd('~/Distance_Statistics/') # cluster
data.path<-('/project/taki2/CURE_Data/') # cluster

#################################################
#################################################
#READ IN DATA
#################################################
#################################################

if (!file.exists('CURE_DTI_2018-11-12.RData')) {

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

	save.image('CURE_DTI_2018-11-12.RData')

} else {
	load('CURE_DTI_2018-11-12.RData')
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
#	if (length(levels(d))>2) {
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
#	} else {
#	  eta<-sum(d==0)/length(d)
#	  h.1.2<-dists.mat*(1-(1-d)%o%(1-d)/(eta) - (d)%o%(d)/(1-eta))
#	  phi.21<-2*eta*((1-d)*exp.dist-eta*sig^2)
#	  phi.22<-(d-(1-eta))/eta^2
#	  h.2.2<-(phi.21%o%phi.22 + phi.22%o%phi.21)/2
#	  phi.31<-2*(1-eta)*(d*exp.dist-(1-eta)*sig^2) 
#	  phi.32<--(1-eta)^(-2)*(d-(1-eta))
#	  h.3.2<-(phi.31%o%phi.32 + phi.32%o%phi.31)/2
#	  h.4.2<-sig^2/((1-eta)*eta)*(d-(1-eta))%o%(d-(1-eta))
#	  out.mat<-h.1.2-h.2.2-h.3.2-h.4.2
#	  return(out.mat)
#	}
}

get.ganova.approx<-function(evals,K,N=1E4){
	tmp<-rep(0,N)
	for (j in 1:K) tmp<-tmp+(evals[j]*(rchisq(N,1)-1))
	return(tmp)
}

#This function does the proposed testing.
do.proposed.test<-function(n,dd,d.mat,B.mc=2.5E5,return.details=FALSE,return.sample=TRUE,return.sample.size=5000) {

	SS<-mean(d.mat)*n
	SSE<-0; dd.levels<-sort(unique(dd))
	for (ii in 1:length(dd.levels)) SSE<-SSE+mean(d.mat[(dd==dd.levels[ii])%o%(dd==dd.levels[ii])==1])*sum(1*dd==1*dd.levels[ii]) 
	
	numer<-SS-SSE
	denom<-SSE/(n-2)
	ganova.F<-numer/denom

	A.mat<-get.A.mat(d.mat,dd,sqrt(SS/n))
	e.A.mat<-eigen(A.mat)
	#evals<-e.A.mat$values[order(abs(e.A.mat$values),decreasing=TRUE)][1:10]
	evals.all<-e.A.mat$values[order(abs(e.A.mat$values),decreasing=TRUE)]
	n.evals<-sum(cumsum(evals.all)/sum(evals.all)<0.95)
	evals<-e.A.mat$values[order(abs(e.A.mat$values),decreasing=TRUE)][1:n.evals]
	numer.approx<-get.ganova.approx(evals,K=length(evals),N=B.mc)/n+(SS/n)
	p.val<-mean(ganova.F<numer.approx/denom)

### RETURN RESULTS
	if (return.details==TRUE) {
		return(list(SS=SS,SSE=SSE,numer=numer,denom=denom,ganova.F=ganova.F,evals=evals,evals.all=evals.all,numer.approx=numer.approx,p.val=p.val))
	} else {
		if (return.sample==TRUE) {
			return(list(ganova.F=ganova.F,p.val=p.val,asympt.est.sample=numer.approx[1:return.sample.size]/denom))
		} else {
			return(list(ganova.F=ganova.F,p.val=p.val))
		}
	}

}

mat.to.list<-function(mat) lapply(seq_len(ncol(mat)), function(i) mat[,i])
vec.d<-function(x1,x2) sum(x1-x2)^2


#################################################
#################################################
#DO ANALYSIS
#################################################
#################################################

### ANALYSIS FOR AGE EFFECT - 2 GROUPS

get.upper.tri<-function(mat) mat[upper.tri(mat)]
list.vecs<-lapply(dti.data,get.upper.tri)
dti.mat<-matrix(unlist(list.vecs),ncol=length(list.vecs))

proposed.dmat.time<-system.time(d.mat<-get.d.mat.list(mat.to.list(dti.mat),d=vec.d))[[3]]
proposed.test.time.age<-system.time(proposed.test.age<-do.proposed.test(n=n,dd=(dti.age>med.age),d.mat=d.mat,B.mc=2.5E5,return.details=TRUE))[[3]]+proposed.dmat.time
age.factor<-factor(dti.age>med.age)
competitor.test.time.age<-system.time(adonis.p.age<-adonis(as.dist(d.mat)~age.factor,permutations=2.5E5)$aov[6][[1]][1])[[3]]+proposed.dmat.time

print('Age comparisons with 2 groups:')
print(paste('The proposed test p-value is ', proposed.test.age$p.val ,'...'))
print(paste('The adonis p-value is ', adonis.p.age ,'...'))

print(paste('The computation time was ', proposed.test.time.age,'...'))
print(paste('The adonis computation time was ', competitor.test.time.age,'...'))

### ANALYSIS FOR DISEASE EFFECT

proposed.test.time.dx<-system.time(proposed.test.dx<-do.proposed.test(n=n,dd=(dti.group=='ASD'),d.mat=d.mat,B.mc=2.5E5,return.details=TRUE))[[3]]+proposed.dmat.time
competitor.test.time.dx<-system.time(adonis.p.dx<-adonis(as.dist(d.mat)~dti.group,permutations=2.5E5)$aov[6][[1]][1])[[3]]+proposed.dmat.time

print('Disease comparisons:')
print(paste('The proposed test p-value is ', proposed.test.dx$p.val ,'...'))
print(paste('The adonis p-value is ', adonis.p.dx ,'...'))

print(paste('The computation time was ', proposed.test.time.dx,'...'))
print(paste('The adonis computation time was ', competitor.test.time.dx,'...'))


### ANALYSIS FOR AGE EFFECT - 5 GROUPS

dd.age<-(dti.age>quantile(dti.age,0.2,na.rm=TRUE))+(dti.age>quantile(dti.age,0.4,na.rm=TRUE))+(dti.age>quantile(dti.age,0.6,na.rm=TRUE))+(dti.age>quantile(dti.age,0.8,na.rm=TRUE))
proposed.test.time.age.5<-system.time(proposed.test.age.5<-do.proposed.test(n=n,dd=dd.age,d.mat=d.mat,B.mc=2.5E5,return.details=TRUE))[[3]]+proposed.dmat.time
age.factor<-factor(dd.age)
competitor.test.time.age.5<-system.time(adonis.p.age.5<-adonis(as.dist(d.mat)~age.factor,permutations=2.5E5)$aov[6][[1]][1])[[3]]+proposed.dmat.time

print('Age comparisons with 5 groups:')
print(paste('The proposed test p-value is ', proposed.test.age.5$p.val ,'...'))
print(paste('The adonis p-value is ', adonis.p.age.5 ,'...'))

print(paste('The computation time was ', proposed.test.time.age.5,'...'))
print(paste('The adonis computation time was ', competitor.test.time.age.5,'...'))


### ANALYSIS FOR AGE EFFECT - 10 GROUPS
dd.age<-(dti.age>quantile(dti.age,0.1,na.rm=TRUE))+(dti.age>quantile(dti.age,0.2,na.rm=TRUE))+(dti.age>quantile(dti.age,0.3,na.rm=TRUE))+(dti.age>quantile(dti.age,0.4,na.rm=TRUE))+(dti.age>quantile(dti.age,0.5,na.rm=TRUE))+(dti.age>quantile(dti.age,0.6,na.rm=TRUE))+(dti.age>quantile(dti.age,0.7,na.rm=TRUE))+(dti.age>quantile(dti.age,0.8,na.rm=TRUE))+(dti.age>quantile(dti.age,0.9,na.rm=TRUE))
proposed.test.time.age.10<-system.time(proposed.test.age.10<-do.proposed.test(n=n,dd=dd.age,d.mat=d.mat,B.mc=2.5E5,return.details=TRUE))[[3]]+proposed.dmat.time
age.factor<-factor(dd.age)
competitor.test.time.age.10<-system.time(adonis.p.age.10<-adonis(as.dist(d.mat)~age.factor,permutations=2.5E5)$aov[6][[1]][1])[[3]]+proposed.dmat.time

print('Age comparisons with 10 groups:')
print(paste('The proposed test p-value is ', proposed.test.age.10$p.val ,'...'))
print(paste('The adonis p-value is ', adonis.p.age.10 ,'...'))

print(paste('The computation time was ', proposed.test.time.age.10,'...'))
print(paste('The adonis computation time was ', competitor.test.time.age.10,'...'))

save(proposed.test.age,adonis.p.age,proposed.test.time.age,competitor.test.time.age,proposed.test.age.5,adonis.p.age.5,proposed.test.time.age.5,competitor.test.time.age.5,proposed.test.age.10,adonis.p.age.10,proposed.test.time.age.10,competitor.test.time.age.10,proposed.test.dx,adonis.p.dx,proposed.test.time.dx,competitor.test.time.dx,file=paste0('CURE_DTI_results_2018-11-18.RData'))
