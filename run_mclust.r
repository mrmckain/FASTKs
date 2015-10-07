
runmclust<-function(ksfile,comp, species) {
	kaks<-read.table(file=ksfile, header=T);
	ks<-kaks[(kaks[,colnames(kaks)=="Ks"]<=3.0),colnames(kaks)=="Ks"];
	ks<-ks[ks>=0.02];
	ks<-log(ks);
	logksclust<-Mclust(ks, G=1:comp);
	ks<-exp(ks);
	ksclust<-Mclust(ks,G=1:comp);
	if(logksclust$bic < ksclust$bic){
		ksclust<-logksclust;}
	ksclust$G;
	bin<-0.05*seq(0,60);
	
	h.kst<-hist(ks,breaks=bin,plot=F)
	#h.kst;
	nc<-h.kst$counts;
	maxnc<-max(nc)+0.5*max(nc);
	vx<-h.kst$mids;
	ntot<-sum(nc);
	par(mai=c(0.5,0.5,0,0)); #set margin for plot bottom,left top, right
	par(pin=c(2.5,2.5)); #plot dimension in inches
	g<-fitden(ksclust,ks);
	h<-ntot*2.5/sum(g);
	vx<-seq(1,100)*0.03;
	ymax<-max(nc)+5;
	myxlab=expression(italic('Ks'));
	#titlename=expression(italic(species));
	barplot(nc,space=0.25,offset=0,width=0.04,xlim=c(0,3),ylim=c(0,maxnc), main=species, las=1, xlab=myxlab, ylab="Number of Pairs");
	#barplot(nc,density=0,space=0,offset=0,width=0.1,xlim=c(0,2),ylim=c(0,ymax),ylab="");
	axis(1); # add x-axis	#return(ks, ksfile.clust);
	color<-c('blue','green','purple','black', 'orange','yellow','darkgreen','darksalmon','deepskyblue','deeppink');
	for ( i in 1:ksclust$G) {
	lines(vx,g[,i]*h,lwd=2,col=color[i]);
	}
}
	
fitden<-function(ksclust,ks){
		comp=ksclust$G;
		mu<-matrix(nrow=1, ncol=comp);
		var<-matrix(nrow=1,ncol=comp);
		pi<-matrix(nrow=1, ncol=comp);
		for (i in 1:comp){
				mu[1,i]<-mean(ks[ksclust$classification==i]);
				var[1,i]<-var(ks[ksclust$classification==i]);
				var[1,i]<-var[1,i]/mu[1,i]^2;
				mu[1,i]=log(mu[1,i]);
				pi[1,i]=length(ks[ksclust$classification==i])/length(ks);
				}
		vx<-seq(1,100)*0.03;
		fx<-matrix(0,100,comp);
		for (i in 1:100) {
	    	for (j in 1:comp) {
	      		fx[i,j]<-pi[j]*dlnorm(vx[i],meanlog=mu[j],sdlog=(sqrt(var[j])));
	    	}
		}
		fx;
}
