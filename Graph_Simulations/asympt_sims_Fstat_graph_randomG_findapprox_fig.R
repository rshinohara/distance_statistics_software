par(mfrow=c(3,2))
load('ganova_graph_2017-10-18_approx_times.RData')
x.lims<-c(0.99*min(c(time.proposed)),1.01*max(c(time.proposed)))
y.lims<-c(min(c(p.err.competitor,p.err.proposed)),max(c(p.err.competitor,p.err.proposed)))
eps.x<-(x.lims[2]-x.lims[1])/30;eps.y<-(y.lims[2]-y.lims[1])/30
x.lims[2]<-x.lims[2]+eps.x;y.lims[2]<-y.lims[2]+eps.y
plot(time.proposed,p.err.proposed,type='l',xlim=x.lims,ylim=y.lims,lwd=3,col=cols[2],xlab='Average Test Time (s)',ylab='p-value MSE',main='Proposed Test')
points(time.proposed,p.err.proposed,pch=2)
text(x=time.proposed+eps.x,y=p.err.proposed+eps.y,labels=mc.grid)
x.lims<-c(0.5*min(c(time.competitor)),1.1*max(c(time.competitor)))
plot(time.competitor,p.err.competitor,type='l',xlim=x.lims,ylim=y.lims,lwd=3,col=cols[2],xlab='Average Test Time (s)',ylab='p-value MSE',main='Permutation Test')
points(time.competitor,p.err.competitor,pch=2)
text(x=time.competitor+eps.x,y=p.err.competitor+eps.y,labels=perm.grid)
load('ganova_graph_approx_times_largeeffect_2017-10-18.RData')
x.lims<-c(0.99*min(c(time.proposed)),1.01*max(c(time.proposed)))
y.lims<-c(min(c(p.err.competitor,p.err.proposed)),max(c(p.err.competitor,p.err.proposed)))
eps.x<-(x.lims[2]-x.lims[1])/30;eps.y<-(y.lims[2]-y.lims[1])/30
x.lims[2]<-x.lims[2]+eps.x;y.lims[2]<-y.lims[2]+eps.y
plot(time.proposed,p.err.proposed,type='l',xlim=x.lims,ylim=y.lims,lwd=3,col=cols[2],xlab='Average Test Time (s)',ylab='p-value MSE',main='Proposed Test')
points(time.proposed,p.err.proposed,pch=2)
text(x=time.proposed+eps.x,y=p.err.proposed+eps.y,labels=mc.grid)
x.lims<-c(0.5*min(c(time.competitor)),1.1*max(c(time.competitor)))
plot(time.competitor,p.err.competitor,type='l',xlim=x.lims,ylim=y.lims,lwd=3,col=cols[2],xlab='Average Test Time (s)',ylab='p-value MSE',main='Permutation Test')
points(time.competitor,p.err.competitor,pch=2)
text(x=time.competitor+eps.x,y=p.err.competitor+eps.y,labels=perm.grid)
load('ganova_graph_approx_times_smalleffect_2017-10-18.RData')
x.lims<-c(0.99*min(c(time.proposed)),1.01*max(c(time.proposed)))
y.lims<-c(min(c(p.err.competitor,p.err.proposed)),max(c(p.err.competitor,p.err.proposed)))
eps.x<-(x.lims[2]-x.lims[1])/30;eps.y<-(y.lims[2]-y.lims[1])/30
x.lims[2]<-x.lims[2]+eps.x;y.lims[2]<-y.lims[2]+eps.y
plot(time.proposed,p.err.proposed,type='l',xlim=x.lims,ylim=y.lims,lwd=3,col=cols[2],xlab='Average Test Time (s)',ylab='p-value MSE',main='Proposed Test')
points(time.proposed,p.err.proposed,pch=2)
text(x=time.proposed+eps.x,y=p.err.proposed+eps.y,labels=mc.grid)
x.lims<-c(0.5*min(c(time.competitor)),1.1*max(c(time.competitor)))
plot(time.competitor,p.err.competitor,type='l',xlim=x.lims,ylim=y.lims,lwd=3,col=cols[2],xlab='Average Test Time (s)',ylab='p-value MSE',main='Permutation Test')
points(time.competitor,p.err.competitor,pch=2)
text(x=time.competitor+eps.x,y=p.err.competitor+eps.y,labels=perm.grid)

dev.off()



par(mfrow=c(3,2))
perm.grid<-c(2.5E4,5E4,7.5E4,1E5,2.5E5)
mc.grid<-c(1E5,2.5E5,5E5,7.5E5,1E6)

load('ganova_graph_2017-10-18_approx_times.RData')
df<-rbind(cbind(perm.grid,sqrt(p.err.competitor),time.competitor,'Null','Permutation'),cbind(mc.grid,sqrt(p.err.proposed),time.proposed,'Null','Proposed'))
load('ganova_graph_approx_times_smalleffect_2017-10-18.RData')
df<-rbind(df,cbind(perm.grid,sqrt(p.err.competitor),time.competitor,'Small Effect','Permutation'),cbind(mc.grid,sqrt(p.err.proposed),time.proposed,'Small Effect','Proposed'))
load('ganova_graph_approx_times_largeeffect_2017-10-18.RData')
df<-rbind(df,cbind(perm.grid,sqrt(p.err.competitor),time.competitor,'Large Effect','Permutation'),cbind(mc.grid,sqrt(p.err.proposed),time.proposed,'Large Effect','Proposed'))
df<-data.frame(df); colnames(df)<-c('parameter','rmse','time','eff.size','method')
df$rmse<-as.numeric(levels(df$rmse)[df$rmse])
df$parameter<-as.numeric(levels(df$parameter)[df$parameter])
df$time<-as.numeric(levels(df$time)[df$time])

df$eff.size<-factor(df$eff.size,levels=c('Null','Small Effect','Large Effect'))

library(ggplot2)
means.barplot <- qplot(x=parameter, y=rmse, fill=method, data=df, geom="bar", stat="identity",position="dodge")

# p = ggplot(data = df, aes(x = parameter, y = rmse, fill=method))
# p = p + geom_bar(stat='identity',position="dodge")
# p = p + facet_grid(~eff.size, scale='free_x')
# p = p + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5,colour='gray50'))
# p


p = ggplot(data = df, aes(x = time, y = rmse, colour=method))
p = p + geom_line() + geom_point() + geom_text(aes(label=parameter,hjust=-0.2, vjust=0), size=2)
p = p + xlab('Computation Time (s)') + ylab('p-value RMSE')  + scale_colour_discrete(name="Method")
p = p + expand_limits(x = c(0, max(df$time)*1.1)) + facet_grid(~eff.size, scale='free_x')
p = p + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5,colour='gray50'))
p

x.lims<-c(0.99*min(c(time.proposed)),1.01*max(c(time.proposed)))
y.lims<-c(min(c(p.err.competitor,p.err.proposed)),max(c(p.err.competitor,p.err.proposed)))
eps.x<-(x.lims[2]-x.lims[1])/30;eps.y<-(y.lims[2]-y.lims[1])/30
x.lims[2]<-x.lims[2]+eps.x;y.lims[2]<-y.lims[2]+eps.y
plot(time.proposed,p.err.proposed,type='l',xlim=x.lims,ylim=y.lims,lwd=3,col=cols[2],xlab='Average Test Time (s)',ylab='p-value MSE',main='Proposed Test')
points(time.proposed,p.err.proposed,pch=2)
text(x=time.proposed+eps.x,y=p.err.proposed+eps.y,labels=mc.grid)
x.lims<-c(0.5*min(c(time.competitor)),1.1*max(c(time.competitor)))
plot(time.competitor,p.err.competitor,type='l',xlim=x.lims,ylim=y.lims,lwd=3,col=cols[2],xlab='Average Test Time (s)',ylab='p-value MSE',main='Permutation Test')
points(time.competitor,p.err.competitor,pch=2)
text(x=time.competitor+eps.x,y=p.err.competitor+eps.y,labels=perm.grid)

x.lims<-c(0.99*min(c(time.proposed)),1.01*max(c(time.proposed)))
y.lims<-c(min(c(p.err.competitor,p.err.proposed)),max(c(p.err.competitor,p.err.proposed)))
eps.x<-(x.lims[2]-x.lims[1])/30;eps.y<-(y.lims[2]-y.lims[1])/30
x.lims[2]<-x.lims[2]+eps.x;y.lims[2]<-y.lims[2]+eps.y
plot(time.proposed,p.err.proposed,type='l',xlim=x.lims,ylim=y.lims,lwd=3,col=cols[2],xlab='Average Test Time (s)',ylab='p-value MSE',main='Proposed Test')
points(time.proposed,p.err.proposed,pch=2)
text(x=time.proposed+eps.x,y=p.err.proposed+eps.y,labels=mc.grid)
x.lims<-c(0.5*min(c(time.competitor)),1.1*max(c(time.competitor)))
plot(time.competitor,p.err.competitor,type='l',xlim=x.lims,ylim=y.lims,lwd=3,col=cols[2],xlab='Average Test Time (s)',ylab='p-value MSE',main='Permutation Test')
points(time.competitor,p.err.competitor,pch=2)
text(x=time.competitor+eps.x,y=p.err.competitor+eps.y,labels=perm.grid)
load('ganova_graph_approx_times_smalleffect_2017-10-18.RData')
x.lims<-c(0.99*min(c(time.proposed)),1.01*max(c(time.proposed)))
y.lims<-c(min(c(p.err.competitor,p.err.proposed)),max(c(p.err.competitor,p.err.proposed)))
eps.x<-(x.lims[2]-x.lims[1])/30;eps.y<-(y.lims[2]-y.lims[1])/30
x.lims[2]<-x.lims[2]+eps.x;y.lims[2]<-y.lims[2]+eps.y
plot(time.proposed,p.err.proposed,type='l',xlim=x.lims,ylim=y.lims,lwd=3,col=cols[2],xlab='Average Test Time (s)',ylab='p-value MSE',main='Proposed Test')
points(time.proposed,p.err.proposed,pch=2)
text(x=time.proposed+eps.x,y=p.err.proposed+eps.y,labels=mc.grid)
x.lims<-c(0.5*min(c(time.competitor)),1.1*max(c(time.competitor)))
plot(time.competitor,p.err.competitor,type='l',xlim=x.lims,ylim=y.lims,lwd=3,col=cols[2],xlab='Average Test Time (s)',ylab='p-value MSE',main='Permutation Test')
points(time.competitor,p.err.competitor,pch=2)
text(x=time.competitor+eps.x,y=p.err.competitor+eps.y,labels=perm.grid)
