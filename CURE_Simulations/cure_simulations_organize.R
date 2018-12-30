setwd('/project/taki3/Distance_Statistics')

mu.vals<-c(0,7.5,15,30)
ks<-c(20,100,200)
k.vals<-1:3
n.vals<-c(30,50,100,264)
p<-1/2

print(paste('For p=',p,'...'))

rejection.rate.proposed<-array(NA,dim=c(4,3,4))
rejection.rate.competitor<-array(NA,dim=c(4,3,4))
time.proposed<-array(NA,dim=c(4,3,4))
time.competitor<-array(NA,dim=c(4,3,4))

dimnames(rejection.rate.proposed)[[1]]<-paste0('mu=',c(0,7.5,15,30))
dimnames(rejection.rate.proposed)[[2]]<-paste0('k=',c(20,100,200))
dimnames(rejection.rate.proposed)[[3]]<-paste0('n=',c(30,50,100,264))

dimnames(rejection.rate.competitor)[[1]]<-paste0('mu=',c(0,7.5,15,30))
dimnames(rejection.rate.competitor)[[2]]<-paste0('k=',c(20,100,200))
dimnames(rejection.rate.competitor)[[3]]<-paste0('n=',c(30,50,100,264))

dimnames(time.proposed)[[1]]<-paste0('mu=',c(0,7.5,15,30))
dimnames(time.proposed)[[2]]<-paste0('k=',c(20,100,200))
dimnames(time.proposed)[[3]]<-paste0('n=',c(30,50,100,264))

dimnames(time.competitor)[[1]]<-paste0('mu=',c(0,7.5,15,30))
dimnames(time.competitor)[[2]]<-paste0('k=',c(20,100,200))
dimnames(time.competitor)[[3]]<-paste0('n=',c(30,50,100,264))

i1<-1
for (mu in c(0,7.5,15,30)) {
	i2<-1
	for (k in ks) {
		i3<-1
		for (n in c(30,50,100,264)) {
			existing.files<-list.files(pattern=glob2rx(paste0('ganova_simulations_dti_2018*',n,'_mu',mu,'_k',k,'.RData')))
			load(existing.files[length(existing.files)])
			try(rejection.rate.proposed[i1,i2,i3]<-mean(unlist(lapply(sims,function(x) x$proposed.test.p))<0.05,na.rm=TRUE))
			try(rejection.rate.competitor[i1,i2,i3]<-mean(unlist(lapply(sims,function(x) x$adonis.p))<0.05,na.rm=TRUE))
			try(time.competitor[i1,i2,i3]<-mean(unlist(lapply(sims,function(x) x$competitor.test.time)),na.rm=TRUE))
			try(time.proposed[i1,i2,i3]<-mean(unlist(lapply(sims,function(x) x$proposed.test.time)),na.rm=TRUE))
			i3<-i3+1
		}
		i2<-i2+1
	}
	i1<-i1+1
}
save(rejection.rate.proposed,rejection.rate.competitor,time.proposed,time.competitor,file=paste0('ganova_2018-12-29_sims_p',round(p,2),'_CURE.RData'))
