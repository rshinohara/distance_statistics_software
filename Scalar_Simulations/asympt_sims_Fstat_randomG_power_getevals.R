#################################################
#	Estimate number of eigenvalues used.
#	December 8, 2018
#################################################

# setwd('~/Dropbox/Documents/Distance_Statistics/') # laptop
# setwd('/Volumes/RAID Array/Dropbox/Documents/Distance_Statistics/') # desktop
setwd('/project/taki3/Distance_Statistics/') # cluster

mu1<-0
mu2<-1
sig<-2
B<-1E4

### DO SIMULATIONS

for (p in c(1/2,1/3)) {

	print(paste('For p=',p,'...'))

	nevals<-array(NA,dim=c(2,3,4))
	
	dimnames(nevals)[[1]]<-c('euclid','abs')
	dimnames(nevals)[[2]]<-paste0('sig=',c(1,2,3))
	dimnames(nevals)[[3]]<-paste0('n=',c(30,50,100,500))

	i1<-1
	for (sig in c(1,2,3)) {
		i2<-1
		for (n in c(30,50,100,500)) {
			
			dist.name<-"euclid"
			filename<-paste('ganova_power_2018-11-28_n',n,'_sig',sig,dist.name,'_p',round(p,2),'.RData',sep='')
			existing.files<-list.files(pattern=glob2rx(paste0('ganova_power_2018-11-28_n',n,'_sig',sig,dist.name,'_p',round(p,2),'.RData')))

			load(existing.files[length(existing.files)])
			

			try(nevals[1,i1,i2]<-mean(unlist(lapply(sims,function(x) length(x$evals))),na.rm=TRUE))
			
			print(paste('For ',dist.name,', sig=',sig,'n, =',n,' the nevals is ', nevals[1,i1,i2]))

	#######

			dist.name<-"abs"
			filename<-paste('ganova_power_2018-11-28_n',n,'_sig',sig,dist.name,'_p',round(p,2),'.RData',sep='')
			existing.files<-list.files(pattern=glob2rx(paste0('ganova_power_2018-11-28_n',n,'_sig',sig,dist.name,'_p',round(p,2),'.RData')))
			load(existing.files[length(existing.files)])

			try(nevals[2,i1,i2]<-mean(unlist(lapply(sims,function(x) length(x$evals))),na.rm=TRUE))
			
			print(paste('For ',dist.name,', sig=',sig,'n, =',n,' the nevals is ', nevals[2,i1,i2]))
			i2<-i2+1
		}
		i1<-i1+1
	}

	print(nevals)
	save(nevals,file=paste0('ganova_2018-12-08_power_nevals_sims','_p',round(p,2),'.RData'))
}
