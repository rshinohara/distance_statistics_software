#################################################
#	Estimate number of eigenvalues used.
#	December 28, 2018
#################################################

library(RColorBrewer)

# setwd('~/Dropbox/Documents/Distance_Statistics/') # laptop
# setwd('/Volumes/RAID Array/Dropbox/Documents/Distance_Statistics/') # desktop
setwd('/project/taki3/Distance_Statistics/') # cluster

#################################################
# SIMULATE AND COMPUTE FOR SCALAR CASES
#################################################

#Simulate data under the null
B<-1E4
mu<-0
sig<-1

### DO SIMULATIONS
cols<-brewer.pal(3,'Pastel1')
other.cols<-brewer.pal(9,'YlGnBu')
n.cores<-as.numeric(Sys.getenv("LSB_DJOB_NUMPROC"))

for (p in c(1/2,1/3)) {

	print(paste('For p=',p,'...'))

	nevals<-array(NA,dim=c(2,4))

	dimnames(nevals)[[1]]<-c('euclid','abs')
	dimnames(nevals)[[2]]<-paste0('n=',c(30,50,100,500))

	i1<-1 
	for (dist.name in c('euclid','abs')) {
		i2<-1
		for (n in c(30,50,100,500)) {

			load(paste('ganova_type1err_2018-12-19_dist_',dist.name,'_n',n,'_p',round(p,2),'.RData',sep=''))
			try(nevals[i1,i2]<-mean(unlist(lapply(sims,function(x) length(x$evals))),na.rm=TRUE))
			i2<-i2+1
		}
		i1<-i1+1
	}
	
	print(nevals)
	save(nevals,file=paste0('ganova_2018-12-19_type1err_nevals_sims','_p',round(p,2),'.RData'))
}
