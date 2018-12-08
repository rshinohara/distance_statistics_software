#################################################
#	Estimate number of eigenvalues used.
#	December 8, 2018
#################################################

# setwd('~/Dropbox/Documents/Distance_Statistics/') # laptop
# setwd('/Volumes/RAID Array/Dropbox/Documents/Distance_Statistics/') # desktop
setwd('/project/taki3/Distance_Statistics/') # cluster

tau1<-0
sigma<-3E-1

for (p in c(1/2,1/3)) {

	print(paste('For p=',p,'...'))

	nevals<-array(NA,dim=c(3,4))

	dimnames(nevals)[[1]]<-paste0('tau=2',c(0.1,0.3,0.5))
	dimnames(nevals)[[2]]<-paste0('n=',c(30,50,100,250))

	i1<-1
	for (tau2 in c(0.1,0.3,0.5)) {
	### ASSESS SIMULATIONS
		i2<-1
		for (n in c(30,50,100,250)) {

			### DO SIMULATIONS
			existing.files<-list.files(pattern=glob2rx(paste0('ganova_power_function_2018-11-25_n',n,'_sigma',sigma,'_tau1',tau1,'_tau2',tau2,'_p',round(p,2),'.RData')))
			load(existing.files[length(existing.files)])
			
			try(nevals[i1,i2]<-mean(unlist(lapply(sims,function(x) length(x$evals))),na.rm=TRUE))
			
			print(paste('For n=',n,' the nevals is ', nevals[i1,i2] ,'...'))
			
			print(paste('Finished with n=',n,'...'))
			i2<-i2+1
		}
		i1<-i1+1
	}
	save(nevals,file=paste0('ganova_power_nevals_2018-12-08_sims_function_p',round(p,2),'_sigma',sigma,'.RData'))
}

