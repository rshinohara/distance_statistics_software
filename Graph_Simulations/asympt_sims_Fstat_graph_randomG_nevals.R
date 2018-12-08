#################################################
#	Estimate number of eigenvalues used.
#	December 8, 2018
#################################################

B<-1E3

for (p in c(1/2,1/3)) {

	print(paste('For p=',p,'...'))

	nevals<-array(NA,dim=c(3,4))
	
	dimnames(nevals)[[1]]<-paste0('tau=',c(0.1,0.15,0.2))
	dimnames(nevals)[[2]]<-paste0('n=',c(30,50,100,250))

	i1<-1
	for (tau in c(0.1,0.15,0.2)) {
		i2<-1
		for (n in c(30,50,100,250)) {

			### DO SIMULATIONS
			existing.files<-list.files(pattern=glob2rx(paste0('ganova_type1err_graph_2018-11-28_n',n,'_tau',tau,'_p',round(p,2),'.RData')))
			load(existing.files[length(existing.files)])
			
			try(nevals[i1,i2]<-mean(unlist(lapply(sims,function(x) length(x$evals))),na.rm=TRUE))

			print(paste('Finished with n=',n,', tau=,',tau,'...'))
			i2<-i2+1
		}
		i1<-i1+1
	}
	print(nevals)
	save(nevals,file=paste0('ganova_type1err_nevals_2018-12-08_graphs_sims_p',round(p,2),'.RData'))
}