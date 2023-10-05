#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
#Functions to support quantifying CU biological status. Code is originally from Korman and Engligh 2013 and has been slightly adapted for the CU snapshots#
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

#--------------#
#Functions#
#--------------#
"getCL"=function(x,a,b){
	#Compute a matrix which has CLs for each spawner level in stock-recruit plot

	nr=length(x);nc=length(a)	#define number of rows (# spawner level) and columns (# posterior simulations) for predicted recruits 
	
	r=matrix(nrow=nr,ncol=nc)
	for (j in 1:nc) {r[,j]=x*exp(a[j]-b[j]*x)}

	cl=matrix(nrow=nr,ncol=2)	#Get 95% CLs for the predicted recruits
	for (i in 1:nr) {cl[i,]=quantile(r[i,],prob=c(0.025,0.975))}
	
	return(cl)
}

"NLL_Sgen1"=function(estSgen,aa,b,Smsy){	#minize this function to get stock size that can recover to Smsy in one generation
	pSmsy=estSgen*exp(aa-b*estSgen)
	NLL=(log(pSmsy)-log(Smsy))^2
	return(NLL)
}

"NLL_Sgen2"=function(estSgen,aa,b,Smsy){	#minize this function to get stock size that can recover to Smsy in two generations
	g1=estSgen*exp(aa-b*estSgen)+1
	pSmsy=g1*exp(aa-b*g1)+1
	NLL=(log(pSmsy)-log(Smsy))^2
	return(NLL)
}

"getProd"=function(aa,b,AvgEsc,AvgER){
	#aa <- mypost[,6] 
	#b <- mypost[,12] 
	#Compute production parameters and some simple benchmarks
	b1=aa/b
	Prod=round(exp(aa),digits=1);Rmax=round(b1*exp(aa-1)/aa,digits=0);Smax=round(b1/aa,digits=0);
	Smsy=round(b1*(0.5-0.07*aa),digits=0); Uopt=round(0.5*aa-0.07*aa^2,digits=2)
	point8Smsy=(round(0.8*b1*(0.5-0.07*aa),digits=0))
	PropSo=0.1*b1	#10% of unfished equilibrium escapement
	
	#compute Sgen benchmarks by iteration and resulting status
	Sgen1=vector(length=nsims);Sgen2=Sgen1;status=Sgen1
	pOver=0
	Sgen1Errors=0;Sgen2Errors=0
	for (i in 1:nsims) {	#nsims
		
		Sgen1[i]=NA;Sgen2[i]=NA;status[i]=NA
		
		#Estimate Sgen1, which is stock size that allows recovery to Smsy in one generation
		init=Smsy[i]*.5	
		fit=optim(par=init,fn=NLL_Sgen1,aa=aa[i],b=b[i],Smsy=Smsy[i],method="L-BFGS-B",hessian=F,lower=1,upper=Smsy[i],control=list(fnscale=1))	#500
		if(fit$convergence==0){
			Sgen1[i]=fit$par
		} else {
			#print(c(i,fit$message))	#fit$message
			Sgen1Errors=Sgen1Errors+1
		}
		
 	 if(is.na(Sgen1[i])==T){Sgen1[i]<-0}

		#Compute status for each sim based on average of last X years of available escapement data where X == generation length
		if(is.na(AvgEsc)==F ){ # removed from if statement
			if(AvgEsc<Sgen1[i]) {
				status[i]=1
			} else if (AvgEsc<point8Smsy[i]){
				status[i]=2
			} else {
				status[i]=3
			}
		}
		
		if(is.na(AvgER)==F) if(AvgER>Uopt[i]) pOver=pOver+1
	}
	pOver=pOver/nsims	#proportion of posterior of Uopt that was less than the average exploitation rate
	
	print(c("Sgen1Error",Sgen1Errors))
	
	stat=matrix(nrow=10,ncol=3)
	i=1;stat[i,1]= median(PropSo,na.rm=T);stat[i,2:3]=quantile(PropSo,prob=c(0.025,0.975),na.rm=T)
	#i=2;stat[i,1]= mean(Sgen2,na.rm=T);stat[i,2:3]=quantile(Sgen2,prob=c(0.025,0.975),na.rm=T)
	i=2;stat[i,1]= median(Sgen1,na.rm=T);stat[i,2:3]=quantile(Sgen1,prob=c(0.025,0.975),na.rm=T)
	i=3;stat[i,1]= median(point8Smsy);stat[i,2:3]=quantile(point8Smsy,prob=c(0.025,0.975))
	i=4;stat[i,1]= median(Smax);stat[i,2:3]=quantile(Smax,prob=c(0.025,0.975))
	i=5;stat[i,1]= median(Prod);stat[i,2:3]=quantile(Prod,prob=c(0.025,0.975))
	i=6;stat[i,1]= median(Uopt);stat[i,2:3]=quantile(Uopt,prob=c(0.025,0.975))
	i=7;stat[i,1]= length(which(status==1))/nsims;stat[i,2]=length(which(status==2))/nsims;stat[i,3]=length(which(status==3))/nsims
	i=8;stat[i,1]= pOver
	i=9;stat[i,1]= 1-exp(-mean(aa));stat[i,2]=1-exp(-quantile(aa,prob=0.025));stat[i,3]=1-exp(-quantile(aa,prob=0.975))
	i=10;stat[i,1]= AvgEsc/mean(b1);stat[i,2]=AvgEsc/quantile(b1,prob=0.025);stat[i,3]=AvgEsc/quantile(b1,prob=0.975)
	
	return(stat)
}
