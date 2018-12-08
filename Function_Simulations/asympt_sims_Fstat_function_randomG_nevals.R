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

	dimnames(nevals)[[1]]<-paste0('sigma=',c(1E-1,3E-1,5E-1))
	dimnames(nevals)[[2]]<-paste0('n=',c(30,50,100,250))

	i1<-1
	for (sigma in c(1E-1,3E-1,5E-1)) {
		i2<-1
		for (n in c(30,50,100,250)) {

			set.seed(976741)
			### DO SIMULATIONS
			existing.files<-list.files(pattern=glob2rx(paste0('ganova_type1err_function_2018-11-25_n',n,'_sigma',sigma,'_p',round(p,2),'.RData')))
			load(existing.files[length(existing.files)])
			
			try(nevals[i1,i2]<-mean(unlist(lapply(sims,function(x) length(x$evals))),na.rm=TRUE))
			
			print(paste('Finished with n=',n,', sigma=,',sigma,'...'))
			i2<-i2+1
		}
		i1<-i1+1
	}
	print(nevals)
	save(nevals,file=paste0('ganova_type1err_nevals_2018-12-08_sims_p',round(p,2),'_function.RData'))
}