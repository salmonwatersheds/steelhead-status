rm(list=ls(all=TRUE))
graphics.off()
if(exists(".SavedPlots",where=1)==T){rm(.SavedPlots,pos=1)}
windows(record=T);par(las=1)

"getCL"=function(x,a,b){
	#Compute a matrix which has CLs for each spawner level in stock-recruit plot

	nr=length(x);nc=length(a)	#define number of rows (# spawner level) and columns (# posterior simulations) for predicted recruits 
	
	r=matrix(nrow=nr,ncol=nc)
	for (j in 1:nc) {r[,j]=x*exp(a[j]-b[j]*x)}

	cl=matrix(nrow=nr,ncol=2)	#Get 95% CLs for the predicted recruits
	for (i in 1:nr) {cl[i,]=quantile(r[i,],prob=c(0.025,0.975))}
	
	return(cl)
}

"NLL_Sgen1"=function(estSgen,a,b,Smsy){	#minize this function to get stock size that can recover to Smsy in one generation
	pSmsy=estSgen*exp(a-b*estSgen)
	NLL=(log(pSmsy)-log(Smsy))^2
	return(NLL)
}

"NLL_Sgen2"=function(estSgen,a,b,Smsy){	#minize this function to get stock size that can recover to Smsy in two generations
	g1=estSgen*exp(a-b*estSgen)+1
	pSmsy=g1*exp(a-b*g1)+1
	NLL=(log(pSmsy)-log(Smsy))^2
	return(NLL)
}



"getProd"=function(a,b,AvgEsc,AvgER){
	
	#Compute production parameters and some simple benchmarks
	b1=a/b
	Prod=round(exp(a),digits=1);Rmax=round(b1*exp(a-1)/a,digits=0);Smax=round(b1/a,digits=0);
	Smsy=round(b1*(0.5-0.07*a),digits=0); Uopt=round(0.5*a-0.07*a^2,digits=2)
	PropSo=0.1*b1	#10% of unfished equilibrium escapement
	
	#compute Sgen benchmarks by iteration and resulting status
	Sgen1=vector(length=nsims);Sgen2=Sgen1;status=Sgen1
	pOver=0
	Sgen1Errors=0;Sgen2Errors=0
	for (i in 1:nsims) {	#nsims
		
		Sgen1[i]=NA;Sgen2[i]=NA;status[i]=NA
		
		#Estimate Sgen1, which is stock size that allows recovery to Smsy in one generation
		init=Smsy[i]*.5	
		fit=optim(par=init,fn=NLL_Sgen1,a=a[i],b=b[i],Smsy=Smsy[i],method="L-BFGS-B",hessian=F,lower=1,upper=Smsy[i],control=list(fnscale=1))	#500
		if(fit$convergence==0){
			Sgen1[i]=fit$par
		} else {
			#print(c(i,fit$message))	#fit$message
			Sgen1Errors=Sgen1Errors+1
		}
		
		#Estimate Sgen2, which is stock size that allows recovery to Smsy in two generations
		#init=Sgen1[i]*.1
		#init=Smsy[i]*.1
		#fit=optim(par=init,fn=NLL_Sgen2,a=a[i],b=b[i],Smsy=Smsy[i],method="L-BFGS-B",hessian=F,lower=1,upper=Smsy[i],control=list(fnscale=1))	#300
		#if(fit$convergence==0){
		#	Sgen2[i]=fit$par
		#} else {
			#print(c(i,fit$message))
		#	Sgen2Errors=Sgen2Errors+1
		#}


		#check to make sure all the results make sense
		#if(Sgen2[i]>=Sgen1[i] | Sgen1[i]>=Smsy[i])print(c(i,a[i],b[i],Sgen2[i],Sgen1[i],Smsy[i]))
		
		#Compute status for each sim based on average of last 5 years of available escapement data
		if(is.na(AvgEsc)==F & is.na(Sgen1[i])==F){
			if(AvgEsc<Sgen1[i]) {
				status[i]=1
			} else if (AvgEsc<Smsy[i]){
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
	i=1;stat[i,1]=mean(PropSo,na.rm=T);stat[i,2:3]=quantile(PropSo,prob=c(0.025,0.975),na.rm=T)
	#i=2;stat[i,1]=mean(Sgen2,na.rm=T);stat[i,2:3]=quantile(Sgen2,prob=c(0.025,0.975),na.rm=T)
	i=2;stat[i,1]=mean(Sgen1,na.rm=T);stat[i,2:3]=quantile(Sgen1,prob=c(0.025,0.975),na.rm=T)
	i=3;stat[i,1]=mean(Smsy);stat[i,2:3]=quantile(Smsy,prob=c(0.025,0.975))
	i=4;stat[i,1]=mean(Smax);stat[i,2:3]=quantile(Smax,prob=c(0.025,0.975))
	i=5;stat[i,1]=mean(Prod);stat[i,2:3]=quantile(Prod,prob=c(0.025,0.975))
	i=6;stat[i,1]=mean(Uopt);stat[i,2:3]=quantile(Uopt,prob=c(0.025,0.975))
	i=7;stat[i,1]=length(which(status==1))/nsims;stat[i,2]=length(which(status==2))/nsims;stat[i,3]=length(which(status==3))/nsims
	i=8;stat[i,1]=pOver
	i=9;stat[i,1]=1-exp(-mean(a));stat[i,2]=1-exp(-quantile(a,prob=0.025));;stat[i,3]=1-exp(-quantile(a,prob=0.975))
	i=10;stat[i,1]=AvgEsc/mean(b1);stat[i,2]=AvgEsc/quantile(b1,prob=0.025);stat[i,3]=AvgEsc/quantile(b1,prob=0.975)
	
	return(stat)
}


###########################################################
########## START GRAPHICS PROGRAM HERE ###################

DoStats=T
fnBench="Bench.out"	
fnStatus="Status.out"	
	
#Read in control file which specifies data input filename
i=1;fndata=scan(file="model.ctl",nlines=1,skip=i,what="character")
i=i+2;HBMflag=scan(file="model.ctl",nlines=1,skip=i)
i=i+2;x=scan(file="model.ctl",nlines=1,skip=i)
FBYr=x[1];MinSRpts=x[2]

########Read in Stock-Recruit Data ###########
MaxStocks=scan(file=fndata,nlines=1,skip=1)

#First pass through to pull up all potential stocks, but they may not have enought SR points
d0<-read.table(file=fndata,header=T,skip=3+MaxStocks)
d=subset(d0,is.na(Rec)==F & is.na(Esc)==F & BY>=FBYr)
StNames=unique(d$CU)
Nstocks=length(StNames)
Nyrs=vector(length=Nstocks);NEscYrs=Nyrs
GdSt=StNames;GdSt[1:Nstocks]=rep(NA,Nstocks)
for (i in 1:Nstocks) {			#Get number of years of valid data for each stock
	d1=subset(d,CU==StNames[i])
	Nyrs[i]=dim(d1)[1]
	if (Nyrs[i]>=MinSRpts) GdSt[i]=StNames[i]	#List of stocks with at least min # of SR points (MinSRpts)
}

#Second pass through data to exclude stocks that don't have enought SR points
d0<-read.table(file=fndata,header=T,skip=3+MaxStocks)
d=subset(d0,is.na(Rec)==F & is.na(Esc)==F & BY>=FBYr & is.na(match(CU,GdSt))==F)
StNames=unique(d$CU)
Nstocks=length(StNames)
Nyrs=vector(length=Nstocks)
for (i in 1:Nstocks) {			#Get number of years of valid data for each stock
	d1=subset(d,CU==StNames[i])
	Nyrs[i]=dim(d1)[1]
	
	d1a=subset(d0,CU==StNames[i])	#Need to create an escapement array for all years with escapement data (some of these years don't have recruit data)
	NEscYrs[i]=dim(d1a)[1]
}

MaxYrs=max(Nyrs)			#Assign data to appropriate elements of S and R arrays
S=matrix(nrow=MaxYrs,ncol=Nstocks);R=S	
SR_BY=matrix(nrow=MaxYrs,ncol=Nstocks)
for (i in 1:Nstocks) {
	d1=subset(d,CU==StNames[i])
	S[1:Nyrs[i],i]=d1$Esc
	R[1:Nyrs[i],i]=d1$Rec
	SR_BY[1:Nyrs[i],i]=d1$BY
}

#Also creat arrays holding all escapement data and years of escapement (note some of these years don't have recruitment data) 
MaxEscYrs=max(NEscYrs)			
Esc=matrix(nrow=MaxEscYrs,ncol=Nstocks);BY=Esc;ER=Esc;
for (i in 1:Nstocks) {
	d1a=subset(d0,CU==StNames[i])
	BY[1:NEscYrs[i],i]=d1a$BY
	Esc[1:NEscYrs[i],i]=d1a$Esc
	ER[1:NEscYrs[i],i]=d1a$ER
}

#Set priors on b
d0=read.table(file=fndata,header=T,skip=2,nrows=MaxStocks)
d1=subset(d0,is.na(match(CU,StNames))==F)
prSmax=d1$prSmax;prCV=d1$prCV
prmub=log(1/prSmax)	#convert mean prior on Smax to log b for winbugs model
prtaub=1/prCV^2				#convert from cv to tau
prsdb=1/sqrt(prtaub)
##########################################

fnpost="post.out"
d2<-read.table(file=fnpost,header=T)
nsims=dim(d2)[1]


#### Dump Production Parameter Stats##############################################
AvgEsc=vector(length=Nstocks);AvgER=AvgEsc
Umsy=matrix(nrow=Nstocks,ncol=3);HER=Umsy
Sgen=Umsy;Smsy=Umsy;EscSo=Umsy

write(file=fnBench,x=c("Stat", "CU", "Mean","LCL","UCL"),ncolumns=5,append=F)
write(file=fnStatus,x=c("CU","AvgEsc", "pRed", "pAmber","pGreen","AvgER","Uopt","pOver","S/So"),ncolumns=9,append=F)

for (i in 1:Nstocks) {
	print(as.character(StNames[i]))
	
	#Get appropriate a and b values for current stock
	ii=which(names(d2)==paste("a.",i,sep=""));jj=which(names(d2)==paste("b.",i,sep=""))
	
	#Compute average escapement and Exploitation Rate between 2004 and 2008 for computation of status in getProd function and other plots
	#j=NEscYrs[i]-5+1;k=NEscYrs[i];	AvgEsc[i]=mean(Esc[j:k,i]);AvgER[i]=mean(ER[j:k,i])
	#x=which(BY[,i]>2004 & BY[,i]<2009)
	x=which(BY[,i]>2005)
	AvgEsc[i]=mean(Esc[x,i],na.rm=T);AvgER[i]=mean(ER[x,i],na.rm=T)	#print(c(i,AvgEsc[i],AvgER[i]))
	
	
	HER[i,1]=mean(ER[,i],na.rm=T)	#average ER over period of record
	HER[i,2:3]=quantile(ER[,i],prob=c(0.025,0.975),na.rm=T)
	
	if (DoStats==T){
		bench=round(getProd(d2[,ii],d2[,jj],AvgEsc[i],AvgER[i]),digits=2)
		
		write.table(file=fnBench,cbind(as.character(StNames[i]),bench[c(1:6,9),]),col.names=F,row.names=c("0.1*So","Sgen1","Smsy","Smax","Prod","Uopt","Umax"),append=T)
		write.table(file=fnStatus,cbind(as.character(StNames[i]),AvgEsc[i],t(bench[7,]),AvgER[i],bench[6,1],bench[8,1],bench[10,1]),col.names=F,row.names=F,append=T)
		
		Umsy[i,1:3]=bench[6,1:3]
		Sgen[i,1:3]=bench[2,1:3]
		Smsy[i,1:3]=bench[3,1:3]
		EscSo[i,1:3]=bench[10,1:3]
				
	}
}



############## Plot Full Escapement Time Series #############################################
#ngrows=1;ngcols=1
#ngrows=6;ngcols=3		#SX
#ngrows=4;ngcols=2		#CN
ngrows=2;ngcols=2		#CO,CM
#ngrows=3;ngcols=2		#PK

par(mfcol=c(ngrows,ngcols),xaxs="i",mai=c(.35,.35,.1,.3),omi=c(0.25,0.35,0.25,0))

for (i in 1:Nstocks) {
	if (DoStats==T) {ymax=max(Smsy[i,1]/1000,Esc[1:NEscYrs[i],i]/1000,na.rm=T)} else {ymax=max(Esc[1:NEscYrs[i],i]/1000,na.rm=T)}	

	plot(BY[1:NEscYrs[i],i],Esc[1:NEscYrs[i],i]/1000,ylim=c(0,ymax),xlim=c(1950,2010),type='o',pch=19,bty='n',xlab="Spawners ('000s)",ylab="Recryuts ('000s)",main=StNames[i])#xlim=c(BY[1,i]-1,BY[NEscYrs[i],i]+1)
	
	if (DoStats==T) abline(h=Sgen[i,1]/1000,lty=2,col="red",lwd=2)
	if (DoStats==T) abline(h=Smsy[i,1]/1000,lty=3,col="green",lwd=2)
}
mtext("Brood Year",side = 1,line = 0, outer = T,cex=1,font=2)
par(srt=90);mtext("Escapement ('000s)",side=2, las=3, padj=-1.2,outer=T,cex=1,font=2);	par(srt=0)



############## Plot SR curve based on mean a and b, and 95% CLs along with data###############
par(mfcol=c(ngrows,ngcols),xaxs="i",mai=c(.35,.35,.1,.3),omi=c(0.25,0.35,0.25,0))
#par(mfcol=c(1,1),mai=c(.9,.75,.3,.2),omi=c(0.1,0.1,0.1,0.1))

a=vector(length=Nstocks);b=a
res=matrix(nrow=MaxYrs,ncol=Nstocks)

for (i in 1:Nstocks) {#Nstocks

	Sx=seq(0,max(S[1:Nyrs[i],i]),length.out=200)
	
	j=which(names(d2)==paste("a.",i,sep=""));mu_a=mean(d2[,j])
	k=which(names(d2)==paste("b.",i,sep=""));mu_b=mean(d2[,k])
	
	pR=Sx*exp(mu_a-mu_b*Sx)
	
	SmaxFromPR=prSmax[i]/1000
	#ymax=1.1*max(max(S[1:Nyrs[i],i])/1000,SmaxFromPR)
	ymax=1.1*max(S[1:Nyrs[i],i])/1000
	
	plot(S[1:Nyrs[i],i]/1000,R[1:Nyrs[i],i]/1000,pch=19,xlab="Spawners ('000s)",ylab="Recruits ('000s)",bty='l',xlim=c(0,ymax),ylim=c(0,1.1*max(R[1:Nyrs[i],i])/1000),main=StNames[i])
	
	lines(Sx/1000,pR/1000,lty=1,lwd=2)
	abline(a=0,b=1,lty=3)		#1:1 line
	cl=getCL(Sx,d2[,j],d2[,k])	#plot 95% CLs of SR curve
	lines(Sx/1000,cl[,1]/1000,lty=2,lwd=2);lines(Sx/1000,cl[,2]/1000,lty=2,lwd=2)	
	
	predRec=S[1:Nyrs[i],i]*exp(mu_a-mu_b*S[1:Nyrs[i],i])

	res[1:Nyrs[i],i]=log(R[1:Nyrs[i],i]/1000)-log(predRec/1000)
	
	#Compute linear regression based estimate of parameters and plot
	reg=lm(log(R[1:Nyrs[i],i]/S[1:Nyrs[i],i])~S[1:Nyrs[i],i])
	a[i]=as.double(reg$coefficients[1]);b[i]=as.double(abs(reg$coefficients[2]))#;tau=as.double(1/sd(reg$residuals)^2)
	pR=Sx*exp(a[i]-b[i]*Sx)
	lines(Sx/1000,pR/1000,lty=1,col="gray",lwd=2)
	
	#Plot vertical line to show Smax based on PR model
	#if(i<3 | i>7) abline(v=SmaxFromPR,lty=3,col="red")	#exclude PR prior line for Babine sockeye
	abline(v=SmaxFromPR,lty=3,col="red")
	
	#if(i==1) legend("topleft",legend=c("Bayesian Fit with Prior","95% CLs for Bayesian Fit","Regression Fit - no Prior"),lty=c(1,2,1),lwd=c(2,2,2),col=c("black","black","gray"),bty='n')
}
mtext("Spawners ('000s)",side = 1,line = 0, outer = T,cex=1,font=2)
par(srt=90);mtext("Recruits ('000s)",side=2, las=3, padj=-1.2,outer=T,cex=1,font=2);	par(srt=0)



######### Plot Residuals of SR Curve as a function of time ##########################
par(mfcol=c(ngrows,ngcols),xaxs="i",mai=c(.35,.35,.1,.3),omi=c(0.25,0.35,0.25,0))
for (i in 1:Nstocks) {
	plot(SR_BY[1:Nyrs[i],i],res[1:Nyrs[i],i],type='p',bty='l',pch=19,xlim=c(1950,2010),main=StNames[i])
	fit=lm(res[1:Nyrs[i],i]~SR_BY[1:Nyrs[i],i])

	abline(coef=fit$coefficients,lty=3)
	print(i);print(summary(fit))
}
mtext("Brood Year",side = 1,line = 0, outer = T,cex=1,font=2)
par(srt=90);mtext("Log Residuals (obs-pred)",side=2, las=3, padj=-1.2,outer=T,cex=1,font=2);	par(srt=0)


############# Plot a vs b Scatterplot###############
par(mfcol=c(ngrows,ngcols),xaxs="i",mai=c(.35,.35,.1,.3),omi=c(0.25,0.35,0.25,0))
for (i in 1:Nstocks) {
	ii=which(names(d2)==paste("a.",i,sep=""))
	jj=which(names(d2)==paste("b.",i,sep=""))
	
	plot(d2[,ii],d2[,jj],type='p',bty='l',xlab="",ylab="",xaxt='n',yaxt='n',main=StNames[i])
}
mtext(expression(alpha),side = 1,line = 0, outer = T,cex=1,font=2)
par(srt=90);mtext(expression(beta),side=2, las=3, padj=-1.2,outer=T,cex=1,font=2);	par(srt=0)



########## Plot prior and posterior for b, and also converted to Smax ###########
par(mfcol=c(ngrows,ngcols),xaxs="i",mai=c(.35,.35,.1,.3),omi=c(0.25,0.35,0.25,0))
nbreaks=25
for (i in 1:Nstocks) {	#1:Nstocks
	ii=which(names(d2)==paste("b.",i,sep=""))

	
	x=seq(min(d2[,ii]),max(d2[,ii]),length=nbreaks)		#b
	y=dlnorm(x,prmub[i],prsdb[i]);y1=nsims*y/sum(y)
	ymax=max(y1,max(hist(d2[,ii],breaks=nbreaks,plot=F)$counts))
	hist(d2[,ii],breaks=nbreaks,ylim=c(0,ymax),main=StNames[i],xlab="",ylab="",yaxt='n')
	
	lines(x,y1)			#all other species
	#if(i<3 | i>6) lines(x,y1)	#exclude Babine for SX

	
	if(i==2) legend("topright",legend=c("posterior","prior"),lty=c(NA,1),pch=c(22,NA),bty='n')

	#x=seq(min(1/d2$b),max(1/d2$b),length=nbreaks)					#Smax
	#y=dlnorm(x,log(1E6*prSmax[1]),prb[2]);y1=nsims*y/sum(y)
	#ymax=max(y1,max(hist(1/d2$b,breaks=nbreaks,plot=F)$counts))
	#hist(1/d2$b,breaks=nbreaks,ylim=c(0,ymax),main="",xlab="Smax ('000s)",ylab="",yaxt='n')
	#lines(x,y1)
}
mtext(expression(beta),side = 1,line = 0, outer = T,cex=1,font=2)
par(srt=90);mtext("Probability",side=2, las=3, padj=-1.2,outer=T,cex=1,font=2);	par(srt=0)



######Plot mean hyper distribtion for a and stock-specific estimates with CLs###############
#par(mfcol=c(2,1),mai=c(.5,.5,.1,.3),omi=c(0.25,0.35,0.35,0))
par(mfcol=c(1,1),mai=c(.9,.75,.3,.2),omi=c(0.1,0.1,0.1,0.1))
#x=seq(0,3,length.out=100)	#SX
x=seq(0,5,length.out=100)	#CN

#mean of hyper
mu_a=mean(d2$mu_a);sd_a=mean(d2$sd_a)
y=dlnorm(x,mu_a,sd_a)
#lognormal distribution that fits the independent estimates of a
m=mean(log(a));s=sd(log(a))
y1=dlnorm(x,m,s)
ymax=max(y,y1)

plot(x,y,ylim=c(0,ymax),type='l',bty='l',lty=1,lwd=2,xlab=expression(alpha),ylab="Probability")
lines(x,y1,lty=2,lwd=2)


ymin=min(y);ymax=max(y)
y=ymin
for (i in 1:Nstocks) {
	ii=which(names(d2)==paste("a.",i,sep=""))
	mu=mean(d2[,ii]);lcl=quantile(d2[,ii],prob=0.025);ucl=quantile(d2[,ii],prob=0.975)

	points(mu,y,pch=19)
	lines(x=c(lcl,ucl),y=c(y,y),col="gray",lty=1,lwd=0.5)
	points(a[i],y,pch=21)
	

	#text(0.25,y,StNames[i],cex=0.8)	#sx
	text(4,y,StNames[i],cex=0.8)	#sx
	y=y+(ymax-ymin)/(Nstocks-1)
}



####### Plot a random set of hyper distributions for a ###################
ntrials=100
ip=round(runif(n=ntrials,min=1,max=nsims),digits=0)
i=1
y=dlnorm(x,d2$mu_a[ip[i]],d2$sd_a[ip[i]])
plot(x,y,type='l',bty='l',col="gray",xlab=expression(alpha),ylab="Probability",ylim=c(0,1.5))
for (i in 2:ntrials) {
	y=dlnorm(x,d2$mu_a[ip[i]],d2$sd_a[ip[i]])

	lines(x,y,lty=1,col="gray")
}
y=dlnorm(x,mu_a,sd_a)
lines(x,y,col="black",lwd=2)

mtext(expression(alpha),side = 1,line = 0, outer = T,cex=1,font=2)
par(srt=90);mtext("Probability",side=2, las=3, padj=-1.2,outer=T,cex=1,font=2);	par(srt=0)


#### Plot priors and posteriors for hyper-parameters. Note the priors are hardwired in WB code#########
par(mfcol=c(2,1),mai=c(.95,.5,.25,.3),omi=c(0.25,0.35,0.35,0))
hist(d2$mu_a,prob=T)
sdn=1/sqrt(1.0E-6)
x1=seq(min(d2$mu_a),max(d2$mu_a),length.out=100);yscale=1/max(dnorm(x1,0.5,sdn));y1=yscale*dnorm(x1,0.5,sdn)
lines(x1,y1)
hist(d2$sd_a,prob=T)
x2=seq(min(d2$sd_a),max(d2$sd_a),length.out=100);y2=dgamma(x,0.5,0.5)
lines(x2,y2)


####### Plot random draws of a from hyper distribution, and a converted to Uopt #############
par(mfcol=c(2,1),mai=c(.95,.5,.25,.3),omi=c(0.25,0.35,0.35,0))
ntrials=5000
Uopt=vector(length=ntrials);ra=Uopt
ipick=round(runif(n=ntrials,min=1,max=nsims),digits=0)
for (i in 1:ntrials) {
	mu_a=d2$mu_a[ipick[i]]
	sd_a=d2$sd_a[ipick[i]]
	ra[i]=rlnorm(1,mu_a,sd_a)
	Uopt[i]=0.5*ra[i]-0.07*ra[i]^2
	if(Uopt[i]<0) Uopt[i]=0
}
hist(ra,xlim=c(0,5),breaks=250,xlab=expression(alpha),ylab="",prob=T,main="")
hist(Uopt,xlim=c(0,1),breaks=50,xlab="Harvest Rate for MSY",ylab="",main="",prob=T)
par(srt=90);mtext("Probability",side=2, las=3, padj=-1.2,outer=T,cex=1,font=2);	par(srt=0)


#### Plot average historical exploitation rate and optimal rate from SR curve with 95% clss
par(mfcol=c(1,1),mai=c(.9,.75,.3,.2),omi=c(0.1,0.1,0.1,0.1))
plot(Umsy[,1],HER[,1],type='p',pch=19,xlim=c(0,1),ylim=c(0,1),bty='l',ylab="Historical Total Exploitation Rate",xlab="Exploitation Rate for MSY")
rndy=runif(n=Nstocks,min=0.9,max=1.1)
for (i in 1:Nstocks) {
	lines(y=c(HER[i,1],HER[i,1]),x=c(Umsy[i,2],Umsy[i,3]),col="gray",lty=1,lwd=0.5)
	lines(x=c(Umsy[i,1],Umsy[i,1]),y=c(HER[i,2],HER[i,3]),lty=1,lwd=0.5,col="gray")

	#text(x=Umsy[i,1]*0.95,y=HER[i,1]*1.08+(log(i)*.1),StNames[i],cex=0.8)
	text(x=Umsy[i,1],y=HER[i,1]*rndy[i],StNames[i],cex=0.8)
}
abline(a=0,b=1,lty=3)		#1:1 line


############## Plot Exploitation Time Series and Uopt #############################################
par(mfcol=c(ngrows,ngcols),xaxs="i",mai=c(.35,.35,.1,.3),omi=c(0.25,0.35,0.25,0))
for (i in 1:Nstocks) {
	ymin=min(Umsy[i,2],ER[,i],na.rm=T);ymax=max(Umsy[i,3],ER[,i],na.rm=T)
	
	plot(BY[1:NEscYrs[i],i],ER[1:NEscYrs[i],i],ylim=c(ymin,ymax),xlim=c(1950,2010),type='o',pch=19,bty='n',xlab="",ylab="",main=StNames[i])#xlim=c(BY[1,i]-1,BY[NEscYrs[i],i]+1)
	abline(h=Umsy[i,1],lty=2)
	abline(h=Umsy[i,2],lty=3)
	abline(h=Umsy[i,3],lty=3)
}
mtext("Year",side = 1,line = 0, outer = T,cex=1,font=2)
par(srt=90);mtext("Exploitation Rate",side=2, las=3, padj=-1.2,outer=T,cex=1,font=2);	par(srt=0)



############ Single plot showing abundance and exploitation Status #############################
par(mfcol=c(1,1),mai=c(.9,.75,.3,.2),omi=c(0.1,0.1,0.1,0.1))

Xtype=1
xmin=-0.05;xmax=10;ymin=0;ymax=2.5	#SX, CO, CN, PKe, PKo
#xmin=-0.05;xmax=2;ymin=0;ymax=2		#CN

for (i in 1:Nstocks) {
	if(Xtype==1){
		x=AvgEsc[i]/Sgen[i,1];xlcl=AvgEsc[i]/Sgen[i,3];xucl=AvgEsc[i]/Sgen[i,2]
	}else{
		x=EscSo[i,1];xlcl=EscSo[i,2];xucl=EscSo[i,3]
	}

	
	y=AvgER[i]/Umsy[i,1];ylcl=AvgER[i]/(Umsy[i,3]+.0001);yucl=AvgER[i]/(Umsy[i,2]+.0001);
	
	if (i==1) {
		plot(x,y,pch=19,bty='l',xlab=switch(Xtype,"AvgEsc/Sgen1","AvgEsc/So"),ylab="AvgER/Uopt",xlim=c(xmin,xmax),ylim=c(ymin,ymax))
	} else {
		points(x,y,pch=19)
	}

	lines(x=c(xlcl,xucl),y=c(y,y),col="gray");lines(x=c(x,x),y=c(ylcl,yucl),col="gray")
	text(x=x*0.95,y=y*1.08,StNames[i],cex=0.8)
	#text(x=x,y=y*rndy[i],StNames[i],cex=0.8)
}
abline(v=1,lty=3,lwd=2);abline(h=1,lty=3,lwd=2)
#text(x=1,y=0,"red - underfished",cex=0.7,font=2);text(x=4.5,y=0,"amber or green - underfished",cex=0.7,font=2)
#text(x=1,y=1.5,"red - overfished",cex=0.7,font=2);text(x=4.5,y=ymax,"amber or green - overfished",cex=0.7,font=2)



