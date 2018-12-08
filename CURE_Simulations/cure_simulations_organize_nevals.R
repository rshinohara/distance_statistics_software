setwd('/project/taki3/Distance_Statistics')

mu.vals<-c(0,7.5,15,30)
ks<-c(20,100,200)
k.vals<-1:3
n.vals<-c(30,50,100,264)
p<-1/2

print(paste('For p=',p,'...'))

nevals<-array(NA,dim=c(4,3,4))

dimnames(nevals)[[1]]<-paste0('mu=',c(0,7.5,15,30))
dimnames(nevals)[[2]]<-paste0('k=',c(20,100,200))
dimnames(nevals)[[3]]<-paste0('n=',c(30,50,100,264))


i1<-1
for (mu in c(0,7.5,15,30)) {
	i2<-1
	for (k in ks) {
		i3<-1
		for (n in c(30,50,100,264)) {
			existing.files<-list.files(pattern=glob2rx(paste0('ganova_simulations_dti_*n',n,'_mu',mu,'_k',k,'.RData')))
			load(existing.files[length(existing.files)])
			try(nevals[i1,i2,i3]<-mean(unlist(lapply(sims,function(x) length(x$evals))),na.rm=TRUE))
			i3<-i3+1
		}
		i2<-i2+1
	}
	i1<-i1+1
}
save(nevals,file=paste0('ganova_2018-12-08_sims_nevals_p',round(p,2),'_CURE.RData'))
