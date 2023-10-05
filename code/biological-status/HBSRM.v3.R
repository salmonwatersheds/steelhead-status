#Heiarchical bayesian SR analysis for North Coast salmon conservation units
#Code adpated from Korman and English 2013

#Load required libraries and source code
setwd("~//Dropbox (Salmon Watersheds)/X Drive/1_PROJECTS/Skeena Updates/Skeena 2017 status assessment/HBM and status/")

source("linear_SR_func.R")

library(R2jags) 
library(modeest)       

#set species and constrainst on analysis (first brood year and min # of SR datapoints)
Species <- "PK"
fndata <- paste(Species,"_SRdata.txt",sep="")
FBYr <- -99 # set first brood year, "-99" for no contraint
MinSRpts <- 3 # set minimum # of SR data points required to be included in the analysis

#---------------------------------------------------------------#
# 1. Read in Stock-Recruit Data 
#---------------------------------------------------------------#

MaxStocks=scan(file=fndata,nlines=1,skip=1)

#First pass through to get list of stocks and determine how many years of SR points for each
d0<-read.table(file=fndata,header=T,skip=3+MaxStocks, fill=TRUE)
d=subset(d0,is.na(Rec)==F & is.na(Esc)==F & BY>=FBYr)
StNames=unique(d$CU)
Nstocks=length(StNames)
Nyrs=vector(length=Nstocks)
GdSt=StNames;GdSt[1:Nstocks]=rep(NA,Nstocks)
for (i in 1:Nstocks) {			#Get number of years of valid data for each stock
	d1=subset(d,CU==StNames[i])
	Nyrs[i]=dim(d1)[1]
	if (Nyrs[i]>=MinSRpts) GdSt[i]=StNames[i]
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
}

MaxYrs=max(Nyrs)			#Assign data to appropriate elements of S and R arrays
S=matrix(nrow=MaxYrs,ncol=Nstocks);R=S	
for (i in 1:Nstocks) {
	d1=subset(d,CU==StNames[i])
	S[1:Nyrs[i],i]=d1$Esc
	R[1:Nyrs[i],i]=d1$Rec
}

#Set priors on b
d0=read.table(file=fndata,header=T,skip=2,nrows=MaxStocks)
d1=subset(d0,is.na(match(CU,StNames))==F)
prSmax=d1$prSmax;
prCV=d1$prCV

prmub=log(1/prSmax)	#convert mean prior on Smax to log b for winbugs model
prtaub=1/prCV^2				#convert from cv to tau

####Estimate a and b by linreg and plot
quartz()
LNRS=log(R/S)
inipars=LinReg(Nyrs,LNRS,S,R,StNames)	

#------------------------------------------------------------------------------#
#  Bayes model
#------------------------------------------------------------------------------#
modelFilename = "Bayes_SR_model.txt"
  cat("
model{

#Hyper priors
mu_a~dnorm(0.5,1.0E-6)
tau_a~dgamma(0.5,0.5)
sd_a<-pow(tau_a,-0.5)

for(i in 1:Nstocks) {	
	a[i]~dlnorm(mu_a,tau_a) #Hyper distribution on alpha

	b[i]~dlnorm(prmub[i],prtaub[i])I(1.0E-5,)	#prior on stock-independent b

	sd[i]~dunif(0.05,10)
	tau[i]<-pow(sd[i],-0.5)	

}

for(i in 1:Nstocks) {	
	for(j in 1:Nyrs[i]) {
		LNRS[j,i]~dnorm(Pred[j,i],tau[i])
		Pred[j,i]<-a[i]-b[i]*S[j,i]
	}
}


}
", fill=TRUE, file=modelFilename)

#------------------------------------------------------------------------------#
#  Jags inputs
#------------------------------------------------------------------------------#
jags.data = list("Nstocks","Nyrs","LNRS","S","prmub","prtaub","cov")

jags.parms = c("a","b","b2","sd","mu_a","sd_a","mu_b2","sd_b2")

#------------------------------------------------------------------------------#
#   Run Model
#------------------------------------------------------------------------------#
print("Running Parallel")
ptm = proc.time()
jagsfit.p <- jags.parallel(data=jags.data,  parameters.to.save=jags.parms,n.thin=10,
              n.iter=100000, model.file=modelFilename,n.burnin = 5000,n.chains=6)
endtime = proc.time()-ptm
endtime[3]/60
post = as.mcmc(jagsfit.p)
mypost = as.matrix(post, chain=F)

##### INFERENCE #####
gelman.diag(post, multivariate = F)
model.probs <- round(cbind(est = colMeans(mypost),sd = apply(mypost,2,sd),ci = t(apply(mypost,2,quantile,c(.025,.975)))),digits=8)
model.probs
write.table(mypost,file=paste(Species,".post.out",sep=""),col.names=T,row.names=F)

#------------------------------------------------------------------------------#
#   Plot SR relationship
#------------------------------------------------------------------------------#

sx.post <- read.table(file="SX.post.out",header=T)
pk.post <- read.table(file="PK.post.out",header=T)
cm.post <- read.table(file="CM.post.out",header=T)
co.post <- read.table(file="CO.post.out",header=T)
cn.post <- read.table(file="CN.post.out",header=T)

MaxStocks=scan(file=fndata,nlines=1,skip=1)

#First pass through to get list of stocks and determine how many years of SR points for each
sx.rs<-read.table(file=fndata,header=T,skip=3+MaxStocks, fill=TRUE)



par(mar=c(3.5,3.5,2,2)+.1,mfrow=c(1,1))
#----------------
#sockeye
#---------------
#long
Long.cu <- subset(sx.rs, CU=="Long")

i <- "Long"

alpha <- sx.post[,16]; beta <- sx.post[,65]

max(Long.cu$Esc, na.rm=T)

spw <- seq(0,max(Long.cu$Esc, na.rm=T),100)
spw
recruits <- array(NA,dim=c(1,length(spw),10000))

for(j in 1:10000){
	iter <- sample(length(alpha),1)
	recruits[,,j] <- spw*exp(alpha[iter]-beta[iter]* spw )

}

rec.m <- apply(recruits,c(1,2),median)/1000
rec.u <- apply(recruits,c(1,2),quantile,probs=0.975)/1000
rec.l <- apply(recruits,c(1,2),quantile,probs=0.025)/1000

jpeg("LongSR.jpg", width=8, height=4,units="in",res=200)

par(mfrow=c(1,1),mar=c(1.5,2.5,1,1.5),oma=c(2,2,2,1),new=F)
plot(Long.cu$Esc/1000,Long.cu$Rec/1000,bty='l',xlab="",ylab="",xaxt="n",yaxt="n")
axis(1)
axis(2,las=2)
mtext("Sockeye - Long",3,cex=0.85,line=1)

points(spw/1000, rec.m,type="l" )

polygon(c(spw/1000,rev(spw/1000)),c(rec.u,rev(rec.l)),col=grey(0.9),border=NA)
points(spw/1000, rec.m,type="l",lwd=2 )
points(Long.cu$Esc/1000,Long.cu$Rec/1000)


abline(v=41.4,col="red",lwd=2)
abline(v=29.8,col="red",lty=2)
abline(v=76.7,col="red",lty=2)

abline(v=82.8,col="dark green",lwd=2)
abline(v=59.7,col="dark green",lty=2)
abline(v=153.4,col="dark green",lty=2)
mtext("Spawners (000s)",1,outer=T,line=1)
mtext("Recruits (000s)",2,outer=T,line=1)

dev.off()

#owikeno
Owikeno.cu <- subset(sx.rs, CU=="Owikeno")

alpha <- sx.post[,17]; beta <- sx.post[,66]

max(Owikeno.cu$Esc, na.rm=T)

spw <- seq(0,max(Owikeno.cu$Esc, na.rm=T),10000)
spw
recruits <- array(NA,dim=c(1,length(spw),10000))

for(j in 1:10000){
  iter <- sample(length(alpha),1)
  recruits[,,j] <- spw*exp(alpha[iter]-beta[iter]* spw )
  
}

rec.m <- apply(recruits,c(1,2),median)/1000
rec.u <- apply(recruits,c(1,2),quantile,probs=0.975, na.rm=T)/1000
rec.l <- apply(recruits,c(1,2),quantile,probs=0.025, na.rm=T)/1000

pdf("OwikenoSR.pdf", width=6, height=4)

par(mfrow=c(1,1),mar=c(1.5,2.5,1,1.5),oma=c(2,2,2,1),new=F)
plot(Owikeno.cu$Esc/1000,Owikeno.cu$Rec/1000,bty='l',xlab="",ylab="",xaxt="n",yaxt="n")
axis(1)
axis(2,las=2)
mtext("Sockeye - Owikeno",3,cex=0.85,line=1)

points(spw/1000, rec.m,type="l" )

polygon(c(spw/1000,rev(spw/1000)),c(rec.u,rev(rec.l)),col=grey(0.9),border=NA)
points(spw/1000, rec.m,type="l",lwd=2 )
points(Owikeno.cu$Esc/1000,Owikeno.cu$Rec/1000)


abline(v=275614.25/1000,col="red",lwd=2)
abline(v=0,col="red",lty=2)
abline(v=37500985.95/1000,col="red",lty=2)

		


abline(v=589008.5/1000,col="dark green",lwd=2)
abline(v=340622.55/1000,col="dark green",lty=2)
abline(v=229476000000/1000,col="dark green",lty=2)
mtext("Spawners (000s)",1,outer=T,line=1)
mtext("Recruits (000s)",2,outer=T,line=1)

dev.off()

#backland
Backland.cu <- subset(sx.rs, CU=="Backland")


alpha <- sx.post[,18]; beta <- sx.post[,67]

max(Backland.cu$Esc, na.rm=T)

spw <- seq(0,max(Long.cu$Esc, na.rm=T),100)
spw
recruits <- array(NA,dim=c(1,length(spw),10000))

for(j in 1:10000){
  iter <- sample(length(alpha),1)
  recruits[,,j] <- spw*exp(alpha[iter]-beta[iter]* spw )
  
}

rec.m <- apply(recruits,c(1,2),median)/1000
rec.u <- apply(recruits,c(1,2),quantile,probs=0.975)/1000
rec.l <- apply(recruits,c(1,2),quantile,probs=0.025)/1000

pdf("BacklandSR.pdf", width=6, height=4)

par(mfrow=c(1,1),mar=c(1.5,2.5,1,1.5),oma=c(2,2,2,1),new=F)
plot(Backland.cu$Esc/1000,Backland.cu$Rec/1000,bty='l',ylim=c(0,1.5),xlab="",ylab="",xaxt="n",yaxt="n")
axis(1)
axis(2,las=2)
mtext("Sockeye - Backland",3,cex=0.85,line=1)

points(spw/1000, rec.m,type="l" )

polygon(c(spw/1000,rev(spw/1000)),c(rec.u,rev(rec.l)),col=grey(0.9),border=NA)
points(spw/1000, rec.m,type="l",lwd=2 )
points(Backland.cu$Esc/1000,Backland.cu$Rec/1000)


abline(v=248.62/1000,col="red",lwd=2)
abline(v=36/1000,col="red",lty=2)
abline(v=1884/1000,col="red",lty=2)

abline(v=513/1000,col="dark green",lwd=2)
abline(v=75/1000,col="dark green",lty=2)
abline(v=4026/1000,col="dark green",lty=2)
mtext("Spawners (000s)",1,outer=T,line=1)
mtext("Recruits (000s)",2,outer=T,line=1)

dev.off()


#backland
Canoona.cu <- subset(sx.rs, CU=="Canoona")


alpha <- sx.post[,19]; beta <- sx.post[,68]

max(Canoona.cu$Esc, na.rm=T)

spw <- seq(0,max(Canoona.cu$Esc, na.rm=T),100)
spw
recruits <- array(NA,dim=c(1,length(spw),10000))

for(j in 1:10000){
  iter <- sample(length(alpha),1)
  recruits[,,j] <- spw*exp(alpha[iter]-beta[iter]* spw )
  
}

rec.m <- apply(recruits,c(1,2),median)/1000
rec.u <- apply(recruits,c(1,2),quantile,probs=0.975)/1000
rec.l <- apply(recruits,c(1,2),quantile,probs=0.025)/1000

pdf("CanoonaSR.pdf", width=6, height=4)

par(mfrow=c(1,1),mar=c(1.5,2.5,1,1.5),oma=c(2,2,2,1),new=F)
plot(Canoona.cu$Esc/1000,Canoona.cu$Rec/1000,bty='l',xlab="",ylab="",xaxt="n",yaxt="n")
axis(1)
axis(2,las=2)
mtext("Sockeye - Canoona",3,cex=0.85,line=1)

points(spw/1000, rec.m,type="l" )

polygon(c(spw/1000,rev(spw/1000)),c(rec.u,rev(rec.l)),col=grey(0.9),border=NA)
points(spw/1000, rec.m,type="l",lwd=2 )
points(Canoona.cu$Esc/1000,Canoona.cu$Rec/1000)


abline(v=543/1000,col="red",lwd=2)
abline(v=324/1000,col="red",lty=2)
abline(v=1021/1000,col="red",lty=2)

abline(v=2439/1000,col="dark green",lwd=2)
abline(v=1998/1000,col="dark green",lty=2)
abline(v=3308/1000,col="dark green",lty=2)
mtext("Spawners (000s)",1,outer=T,line=1)
mtext("Recruits (000s)",2,outer=T,line=1)

dev.off()


#evelyn
Evelyn.cu <- subset(sx.rs, CU=="Evelyn")


alpha <- sx.post[,20]; beta <- sx.post[,69]

max(Evelyn.cu$Esc, na.rm=T)

spw <- seq(0,max(Evelyn.cu$Esc, na.rm=T),100)
spw
recruits <- array(NA,dim=c(1,length(spw),10000))

for(j in 1:10000){
  iter <- sample(length(alpha),1)
  recruits[,,j] <- spw*exp(alpha[iter]-beta[iter]* spw )
  
}

rec.m <- apply(recruits,c(1,2),median)/1000
rec.u <- apply(recruits,c(1,2),quantile,probs=0.975)/1000
rec.l <- apply(recruits,c(1,2),quantile,probs=0.025)/1000

pdf("EvelynSR.pdf", width=6, height=4)

par(mfrow=c(1,1),mar=c(1.5,2.5,1,1.5),oma=c(2,2,2,1),new=F)
plot(Evelyn.cu$Esc/1000,Evelyn.cu$Rec/1000,bty='l',xlab="",ylab="",xaxt="n",yaxt="n")
axis(1)
axis(2,las=2)
mtext("Sockeye - Evelyn",3,cex=0.85,line=1)

points(spw/1000, rec.m,type="l" )

polygon(c(spw/1000,rev(spw/1000)),c(rec.u,rev(rec.l)),col=grey(0.9),border=NA)
points(spw/1000, rec.m,type="l",lwd=2 )
points(Evelyn.cu$Esc/1000,Evelyn.cu$Rec/1000)


abline(v=690/1000,col="red",lwd=2)
abline(v=340/1000,col="red",lty=2)
abline(v=7414/1000,col="red",lty=2)

abline(v=1793/1000,col="dark green",lwd=2)
abline(v=1216/1000,col="dark green",lty=2)
abline(v=96329/1000,col="dark green",lty=2)
mtext("Spawners (000s)",1,outer=T,line=1)
mtext("Recruits (000s)",2,outer=T,line=1)

dev.off()

#kainet
Kainet.cu <- subset(sx.rs, CU=="Kainet_Creek")

alpha <- sx.post[,21]; beta <- sx.post[,70]

max(Kainet.cu$Esc, na.rm=T)

spw <- seq(0,max(Kainet.cu$Esc, na.rm=T),100)
spw
recruits <- array(NA,dim=c(1,length(spw),10000))

for(j in 1:10000){
  iter <- sample(length(alpha),1)
  recruits[,,j] <- spw*exp(alpha[iter]-beta[iter]* spw )
  
}

rec.m <- apply(recruits,c(1,2),median)/1000
rec.u <- apply(recruits,c(1,2),quantile,probs=0.975)/1000
rec.l <- apply(recruits,c(1,2),quantile,probs=0.025)/1000

pdf("KainetSR.pdf", width=6, height=4)

par(mfrow=c(1,1),mar=c(1.5,2.5,1,1.5),oma=c(2,2,2,1),new=F)
plot(Kainet.cu$Esc/1000,Kainet.cu$Rec/1000,bty='l',xlab="",ylab="",xaxt="n",yaxt="n")
axis(1)
axis(2,las=2)
mtext("Sockeye - Kainet",3,cex=0.85,line=1)

points(spw/1000, rec.m,type="l" )

polygon(c(spw/1000,rev(spw/1000)),c(rec.u,rev(rec.l)),col=grey(0.9),border=NA)
points(spw/1000, rec.m,type="l",lwd=2 )
points(Kainet.cu$Esc/1000,Kainet.cu$Rec/1000)

abline(v=228/1000,col="red",lwd=2)
abline(v=113/1000,col="red",lty=2)
abline(v=606/1000,col="red",lty=2)

abline(v=1439/1000,col="dark green",lwd=2)
abline(v=1105/1000,col="dark green",lty=2)
abline(v=2437/1000,col="dark green",lty=2)
mtext("Spawners (000s)",1,outer=T,line=1)
mtext("Recruits (000s)",2,outer=T,line=1)

dev.off()

#kitlope
Kitlope.cu <- subset(sx.rs, CU=="Kitlope")

alpha <- sx.post[,22]; beta <- sx.post[,71]

max(Kitlope.cu$Esc, na.rm=T)

spw <- seq(0,max(Kitlope.cu$Esc, na.rm=T),100)
spw
recruits <- array(NA,dim=c(1,length(spw),10000))

for(j in 1:10000){
  iter <- sample(length(alpha),1)
  recruits[,,j] <- spw*exp(alpha[iter]-beta[iter]* spw )
  
}

rec.m <- apply(recruits,c(1,2),median)/1000
rec.u <- apply(recruits,c(1,2),quantile,probs=0.975)/1000
rec.l <- apply(recruits,c(1,2),quantile,probs=0.025)/1000

pdf("KitlopeSR.pdf", width=6, height=4)

par(mfrow=c(1,1),mar=c(1.5,2.5,1,1.5),oma=c(2,2,2,1),new=F)
plot(Kitlope.cu$Esc/1000,Kitlope.cu$Rec/1000,bty='l',xlab="",ylab="",xaxt="n",yaxt="n")
axis(1)
axis(2,las=2)
mtext("Sockeye - Kitlope",3,cex=0.85,line=1)

points(spw/1000, rec.m,type="l" )

polygon(c(spw/1000,rev(spw/1000)),c(rec.u,rev(rec.l)),col=grey(0.9),border=NA)
points(spw/1000, rec.m,type="l",lwd=2 )
points(Kitlope.cu$Esc/1000,Kitlope.cu$Rec/1000)

abline(v=23257/1000,col="red",lwd=2)
abline(v=16600/1000,col="red",lty=2)
abline(v=42939/1000,col="red",lty=2)

abline(v=46542/1000,col="dark green",lwd=2)
abline(v=33264/1000,col="dark green",lty=2)
abline(v=86485/1000,col="dark green",lty=2)
mtext("Spawners (000s)",1,outer=T,line=1)
mtext("Recruits (000s)",2,outer=T,line=1)

dev.off()

#bloomfield
Bloomfield.cu <- subset(sx.rs, CU=="Bloomfield")

alpha <- sx.post[,24]; beta <- sx.post[,73]

max(Bloomfield.cu$Esc, na.rm=T)

spw <- seq(0,max(Bloomfield.cu$Esc, na.rm=T),100)
spw
recruits <- array(NA,dim=c(1,length(spw),10000))

for(j in 1:10000){
  iter <- sample(length(alpha),1)
  recruits[,,j] <- spw*exp(alpha[iter]-beta[iter]* spw )
  
}

rec.m <- apply(recruits,c(1,2),median)/1000
rec.u <- apply(recruits,c(1,2),quantile,probs=0.975)/1000
rec.l <- apply(recruits,c(1,2),quantile,probs=0.025)/1000

pdf("BloomfieldSR.pdf", width=6, height=4)

par(mfrow=c(1,1),mar=c(1.5,2.5,1,1.5),oma=c(2,2,2,1),new=F)
plot(Bloomfield.cu$Esc/1000,Bloomfield.cu$Rec/1000,bty='l',xlab="",ylab="",xaxt="n",yaxt="n")
axis(1)
axis(2,las=2)
mtext("Sockeye - Bloomfield",3,cex=0.85,line=1)

points(spw/1000, rec.m,type="l" )

polygon(c(spw/1000,rev(spw/1000)),c(rec.u,rev(rec.l)),col=grey(0.9),border=NA)
points(spw/1000, rec.m,type="l",lwd=2 )
points(Bloomfield.cu$Esc/1000,Bloomfield.cu$Rec/1000)
541.98	327.19	1161.88	1248	869	2294.05

abline(v=541/1000,col="red",lwd=2)
abline(v=327/1000,col="red",lty=2)
abline(v=1162/1000,col="red",lty=2)

abline(v=1248/1000,col="dark green",lwd=2)
abline(v=869/1000,col="dark green",lty=2)
abline(v=2294/1000,col="dark green",lty=2)
mtext("Spawners (000s)",1,outer=T,line=1)
mtext("Recruits (000s)",2,outer=T,line=1)

dev.off()


#curtis inlet
Curtis.cu <- subset(sx.rs, CU=="Curtis_Inlet")

alpha <- sx.post[,25]; beta <- sx.post[,74]

max(Curtis.cu$Esc, na.rm=T)

spw <- seq(0,max(Curtis.cu$Esc, na.rm=T),100)
spw
recruits <- array(NA,dim=c(1,length(spw),10000))

for(j in 1:10000){
  iter <- sample(length(alpha),1)
  recruits[,,j] <- spw*exp(alpha[iter]-beta[iter]* spw )
  
}

rec.m <- apply(recruits,c(1,2),median)/1000
rec.u <- apply(recruits,c(1,2),quantile,probs=0.975)/1000
rec.l <- apply(recruits,c(1,2),quantile,probs=0.025)/1000

pdf("CurtisInletSR.pdf", width=6, height=4)

par(mfrow=c(1,1),mar=c(1.5,2.5,1,1.5),oma=c(2,2,2,1),new=F)
plot(Curtis.cu$Esc/1000,Curtis.cu$Rec/1000,bty='l',xlab="",ylab="",xaxt="n",yaxt="n")
axis(1)
axis(2,las=2)
mtext("Sockeye - Curtis Inlet",3,cex=0.85,line=1)

points(spw/1000, rec.m,type="l" )

polygon(c(spw/1000,rev(spw/1000)),c(rec.u,rev(rec.l)),col=grey(0.9),border=NA)
points(spw/1000, rec.m,type="l",lwd=2 )
points(Curtis.cu$Esc/1000,Curtis.cu$Rec/1000)
4259.5	2379.43	13326.68	8798	5787	31563.6


abline(v=4259/1000,col="red",lwd=2)
abline(v=2379/1000,col="red",lty=2)
abline(v=13327/1000,col="red",lty=2)

abline(v=8798/1000,col="dark green",lwd=2)
abline(v=5787/1000,col="dark green",lty=2)
abline(v=31563/1000,col="dark green",lty=2)
mtext("Spawners (000s)",1,outer=T,line=1)
mtext("Recruits (000s)",2,outer=T,line=1)

dev.off()

#devon
Devon.cu <- subset(sx.rs, CU=="Devon")

alpha <- sx.post[,26]; beta <- sx.post[,75]

max(Devon.cu$Esc, na.rm=T)

spw <- seq(0,max(Devon.cu$Esc, na.rm=T),100)
spw
recruits <- array(NA,dim=c(1,length(spw),10000))

for(j in 1:10000){
  iter <- sample(length(alpha),1)
  recruits[,,j] <- spw*exp(alpha[iter]-beta[iter]* spw )
  
}

rec.m <- apply(recruits,c(1,2),median)/1000
rec.u <- apply(recruits,c(1,2),quantile,probs=0.975)/1000
rec.l <- apply(recruits,c(1,2),quantile,probs=0.025)/1000

pdf("DevonSR.pdf", width=6, height=4)

par(mfrow=c(1,1),mar=c(1.5,2.5,1,1.5),oma=c(2,2,2,1),new=F)
plot(Devon.cu$Esc/1000,Devon.cu$Rec/1000,bty='l',xlab="",ylab="",xaxt="n",yaxt="n")
axis(1)
axis(2,las=2)
mtext("Sockeye - Devon",3,cex=0.85,line=1)

points(spw/1000, rec.m,type="l" )

polygon(c(spw/1000,rev(spw/1000)),c(rec.u,rev(rec.l)),col=grey(0.9),border=NA)
points(spw/1000, rec.m,type="l",lwd=2 )
points(Devon.cu$Esc/1000,Devon.cu$Rec/1000)
1790.68	1245.93	2634	4283	3589	5341.03



abline(v=1790/1000,col="red",lwd=2)
abline(v=1245/1000,col="red",lty=2)
abline(v=2634/1000,col="red",lty=2)

abline(v=4283/1000,col="dark green",lwd=2)
abline(v=3589/1000,col="dark green",lty=2)
abline(v=5341/1000,col="dark green",lty=2)
mtext("Spawners (000s)",1,outer=T,line=1)
mtext("Recruits (000s)",2,outer=T,line=1)

dev.off()

#freeda
Freeda.cu <- subset(sx.rs, CU=="Freeda")

alpha <- sx.post[,27]; beta <- sx.post[,76]

max(Freeda.cu$Esc, na.rm=T)

spw <- seq(0,max(Freeda.cu$Esc, na.rm=T),100)
spw
recruits <- array(NA,dim=c(1,length(spw),10000))

for(j in 1:10000){
  iter <- sample(length(alpha),1)
  recruits[,,j] <- spw*exp(alpha[iter]-beta[iter]* spw )
  
}

rec.m <- apply(recruits,c(1,2),median)/1000
rec.u <- apply(recruits,c(1,2),quantile,probs=0.975)/1000
rec.l <- apply(recruits,c(1,2),quantile,probs=0.025)/1000

pdf("FreedaSR.pdf", width=6, height=4)

par(mfrow=c(1,1),mar=c(1.5,2.5,1,1.5),oma=c(2,2,2,1),new=F)
plot(Freeda.cu$Esc/1000,Freeda.cu$Rec/1000,bty='l',xlab="",ylab="",xaxt="n",yaxt="n")
axis(1)
axis(2,las=2)
mtext("Sockeye - Freeda",3,cex=0.85,line=1)

points(spw/1000, rec.m,type="l" )

polygon(c(spw/1000,rev(spw/1000)),c(rec.u,rev(rec.l)),col=grey(0.9),border=NA)
points(spw/1000, rec.m,type="l",lwd=2 )
points(Freeda.cu$Esc/1000,Freeda.cu$Rec/1000)
241.46	127.9	523.1	603	401	1104.03




abline(v=241/1000,col="red",lwd=2)
abline(v=127/1000,col="red",lty=2)
abline(v=523/1000,col="red",lty=2)

abline(v=603/1000,col="dark green",lwd=2)
abline(v=401/1000,col="dark green",lty=2)
abline(v=1104/1000,col="dark green",lty=2)
mtext("Spawners (000s)",1,outer=T,line=1)
mtext("Recruits (000s)",2,outer=T,line=1)

dev.off()

#Hartley bay
Hartley.cu <- subset(sx.rs, CU=="Hartley_Bay")

alpha <- sx.post[,28]; beta <- sx.post[,77]

max(Hartley.cu$Esc, na.rm=T)

spw <- seq(0,max(Hartley.cu$Esc, na.rm=T),100)
spw
recruits <- array(NA,dim=c(1,length(spw),10000))

for(j in 1:10000){
  iter <- sample(length(alpha),1)
  recruits[,,j] <- spw*exp(alpha[iter]-beta[iter]* spw )
  
}

rec.m <- apply(recruits,c(1,2),median)/1000
rec.u <- apply(recruits,c(1,2),quantile,probs=0.975)/1000
rec.l <- apply(recruits,c(1,2),quantile,probs=0.025)/1000

pdf("HartleySR.pdf", width=6, height=4)

par(mfrow=c(1,1),mar=c(1.5,2.5,1,1.5),oma=c(2,2,2,1),new=F)
plot(Hartley.cu$Esc/1000,Hartley.cu$Rec/1000,bty='l',xlab="",ylab="",xaxt="n",yaxt="n")
axis(1)
axis(2,las=2)
mtext("Sockeye - Hartley Bay",3,cex=0.85,line=1)

points(spw/1000, rec.m,type="l" )

polygon(c(spw/1000,rev(spw/1000)),c(rec.u,rev(rec.l)),col=grey(0.9),border=NA)
points(spw/1000, rec.m,type="l",lwd=2 )
points(Hartley.cu$Esc/1000,Hartley.cu$Rec/1000)
692.42	0	5075388.84	1920.5	1060	4433004922

abline(v=692/1000,col="red",lwd=2)
abline(v=0/1000,col="red",lty=2)
abline(v=5075388/1000,col="red",lty=2)

abline(v=1920/1000,col="dark green",lwd=2)
abline(v=1060/1000,col="dark green",lty=2)
abline(v=4433004922/1000,col="dark green",lty=2)
mtext("Spawners (000s)",1,outer=T,line=1)
mtext("Recruits (000s)",2,outer=T,line=1)

dev.off()


#Kadjusdis River
Kad.cu <- subset(sx.rs, CU=="Kadjusdis_River")

alpha <- sx.post[,29]; beta <- sx.post[,78]

max(Kad.cu$Esc, na.rm=T)

spw <- seq(0,max(Kad.cu$Esc, na.rm=T),100)
spw
recruits <- array(NA,dim=c(1,length(spw),10000))

for(j in 1:10000){
  iter <- sample(length(alpha),1)
  recruits[,,j] <- spw*exp(alpha[iter]-beta[iter]* spw )
  
}

rec.m <- apply(recruits,c(1,2),median)/1000
rec.u <- apply(recruits,c(1,2),quantile,probs=0.975)/1000
rec.l <- apply(recruits,c(1,2),quantile,probs=0.025)/1000

pdf("KadSR.pdf", width=6, height=4)

par(mfrow=c(1,1),mar=c(1.5,2.5,1,1.5),oma=c(2,2,2,1),new=F)
plot(Kad.cu$Esc/1000,Kad.cu$Rec/1000,bty='l',xlab="",ylab="",xaxt="n",yaxt="n")
axis(1)
axis(2,las=2)
mtext("Sockeye - Kadjusdis River",3,cex=0.85,line=1)

points(spw/1000, rec.m,type="l" )

polygon(c(spw/1000,rev(spw/1000)),c(rec.u,rev(rec.l)),col=grey(0.9),border=NA)
points(spw/1000, rec.m,type="l",lwd=2 )
points(Kad.cu$Esc/1000,Kad.cu$Rec/1000)
1329.22	566.59	4549.49	4872	3018	12378.08


abline(v=1329/1000,col="red",lwd=2)
abline(v=566/1000,col="red",lty=2)
abline(v=4549/1000,col="red",lty=2)

abline(v=4872/1000,col="dark green",lwd=2)
abline(v=3018/1000,col="dark green",lty=2)
abline(v=12378/1000,col="dark green",lty=2)
mtext("Spawners (000s)",1,outer=T,line=1)
mtext("Recruits (000s)",2,outer=T,line=1)

dev.off()


#keecha
Keecha.cu <- subset(sx.rs, CU=="Keecha")

alpha <- sx.post[,30]; beta <- sx.post[,79]

max(Keecha.cu$Esc, na.rm=T)

spw <- seq(0,max(Keecha.cu$Esc, na.rm=T),100)
spw
recruits <- array(NA,dim=c(1,length(spw),10000))

for(j in 1:10000){
  iter <- sample(length(alpha),1)
  recruits[,,j] <- spw*exp(alpha[iter]-beta[iter]* spw )
  
}

rec.m <- apply(recruits,c(1,2),median)/1000
rec.u <- apply(recruits,c(1,2),quantile,probs=0.975)/1000
rec.l <- apply(recruits,c(1,2),quantile,probs=0.025)/1000

pdf("KeechaSR.pdf", width=6, height=4)

par(mfrow=c(1,1),mar=c(1.5,2.5,1,1.5),oma=c(2,2,2,1),new=F)
plot(Keecha.cu$Esc/1000,Keecha.cu$Rec/1000,bty='l',xlab="",ylab="",xaxt="n",yaxt="n")
axis(1)
axis(2,las=2)
mtext("Sockeye - Keecha",3,cex=0.85,line=1)

points(spw/1000, rec.m,type="l" )

polygon(c(spw/1000,rev(spw/1000)),c(rec.u,rev(rec.l)),col=grey(0.9),border=NA)
points(spw/1000, rec.m,type="l",lwd=2 )
points(Keecha.cu$Esc/1000,Keecha.cu$Rec/1000)
1989.96	1085.6	6048.54	4200	2772	12680.23



abline(v=1989/1000,col="red",lwd=2)
abline(v=1086/1000,col="red",lty=2)
abline(v=6048/1000,col="red",lty=2)

abline(v=4200/1000,col="dark green",lwd=2)
abline(v=2772/1000,col="dark green",lty=2)
abline(v=12680/1000,col="dark green",lty=2)
mtext("Spawners (000s)",1,outer=T,line=1)
mtext("Recruits (000s)",2,outer=T,line=1)

dev.off()



#koeye
Koeye.cu <- subset(sx.rs, CU=="Koeye")

alpha <- sx.post[,31]; beta <- sx.post[,80]

max(Koeye.cu$Esc, na.rm=T)

spw <- seq(0,max(Koeye.cu$Esc, na.rm=T),100)
spw
recruits <- array(NA,dim=c(1,length(spw),10000))

for(j in 1:10000){
  iter <- sample(length(alpha),1)
  recruits[,,j] <- spw*exp(alpha[iter]-beta[iter]* spw )
  
}

rec.m <- apply(recruits,c(1,2),median)/1000
rec.u <- apply(recruits,c(1,2),quantile,probs=0.975)/1000
rec.l <- apply(recruits,c(1,2),quantile,probs=0.025)/1000

pdf("KoeyeSR.pdf", width=6, height=4)

par(mfrow=c(1,1),mar=c(1.5,2.5,1,1.5),oma=c(2,2,2,1),new=F)
plot(Koeye.cu$Esc/1000,Koeye.cu$Rec/1000,bty='l',xlab="",ylab="",xaxt="n",yaxt="n")
axis(1)
axis(2,las=2)
mtext("Sockeye - Koeye",3,cex=0.85,line=1)

points(spw/1000, rec.m,type="l" )

polygon(c(spw/1000,rev(spw/1000)),c(rec.u,rev(rec.l)),col=grey(0.9),border=NA)
points(spw/1000, rec.m,type="l",lwd=2 )
points(Koeye.cu$Esc/1000,Koeye.cu$Rec/1000)
2548.92	1032.09	10349.15	6924.5	3852	20824




abline(v=2548/1000,col="red",lwd=2)
abline(v=1032/1000,col="red",lty=2)
abline(v=10349/1000,col="red",lty=2)

abline(v=6924/1000,col="dark green",lwd=2)
abline(v=3852/1000,col="dark green",lty=2)
abline(v=20824/1000,col="dark green",lty=2)
mtext("Spawners (000s)",1,outer=T,line=1)
mtext("Recruits (000s)",2,outer=T,line=1)

dev.off()

#kooryet
Kooryet.cu <- subset(sx.rs, CU=="Kooryet")

alpha <- sx.post[,32]; beta <- sx.post[,81]

max(Kooryet.cu$Esc, na.rm=T)

spw <- seq(0,max(Kooryet.cu$Esc, na.rm=T),100)
spw
recruits <- array(NA,dim=c(1,length(spw),10000))

for(j in 1:10000){
  iter <- sample(length(alpha),1)
  recruits[,,j] <- spw*exp(alpha[iter]-beta[iter]* spw )
  
}

rec.m <- apply(recruits,c(1,2),median)/1000
rec.u <- apply(recruits,c(1,2),quantile,probs=0.975)/1000
rec.l <- apply(recruits,c(1,2),quantile,probs=0.025)/1000

pdf("KooryetSR.pdf", width=6, height=4)

par(mfrow=c(1,1),mar=c(1.5,2.5,1,1.5),oma=c(2,2,2,1),new=F)
plot(Kooryet.cu$Esc/1000,Kooryet.cu$Rec/1000,bty='l',xlab="",ylab="",xaxt="n",yaxt="n")
axis(1)
axis(2,las=2)
mtext("Sockeye - Kooryet",3,cex=0.85,line=1)

points(spw/1000, rec.m,type="l" )

polygon(c(spw/1000,rev(spw/1000)),c(rec.u,rev(rec.l)),col=grey(0.9),border=NA)
points(spw/1000, rec.m,type="l",lwd=2 )
points(Kooryet.cu$Esc/1000,Kooryet.cu$Rec/1000)
2022.17	1097.42	8854.77	4238	2771	29093.85

abline(v=2022/1000,col="red",lwd=2)
abline(v=1097/1000,col="red",lty=2)
abline(v=8854/1000,col="red",lty=2)

abline(v=4238/1000,col="dark green",lwd=2)
abline(v=2771/1000,col="dark green",lty=2)
abline(v=29093/1000,col="dark green",lty=2)
mtext("Spawners (000s)",1,outer=T,line=1)
mtext("Recruits (000s)",2,outer=T,line=1)

dev.off()




#Kwakwa Creek

Kwakwa.cu <- subset(sx.rs, CU=="Kwakwa_Creek")

alpha <- sx.post[,33]; beta <- sx.post[,82]

max(Kwakwa.cu$Esc, na.rm=T)

spw <- seq(0,max(Kwakwa.cu$Esc, na.rm=T),100)
spw
recruits <- array(NA,dim=c(1,length(spw),10000))

for(j in 1:10000){
  iter <- sample(length(alpha),1)
  recruits[,,j] <- spw*exp(alpha[iter]-beta[iter]* spw )
  
}

rec.m <- apply(recruits,c(1,2),median)/1000
rec.u <- apply(recruits,c(1,2),quantile,probs=0.975)/1000
rec.l <- apply(recruits,c(1,2),quantile,probs=0.025)/1000

pdf("KwakwaSR.pdf", width=6, height=4)

par(mfrow=c(1,1),mar=c(1.5,2.5,1,1.5),oma=c(2,2,2,1),new=F)
plot(Kwakwa.cu$Esc/1000,Kwakwa.cu$Rec/1000,bty='l',xlab="",ylab="",xaxt="n",yaxt="n")
axis(1)
axis(2,las=2)
mtext("Sockeye - Kwakwa Creek",3,cex=0.85,line=1)

points(spw/1000, rec.m,type="l" )

polygon(c(spw/1000,rev(spw/1000)),c(rec.u,rev(rec.l)),col=grey(0.9),border=NA)
points(spw/1000, rec.m,type="l",lwd=2 )
points(Kwakwa.cu$Esc/1000,Kwakwa.cu$Rec/1000)
1360.18	819.35	2912.77	3511	2535.98	6166


abline(v=1360/1000,col="red",lwd=2)
abline(v=819/1000,col="red",lty=2)
abline(v=2912/1000,col="red",lty=2)

abline(v=3511/1000,col="dark green",lwd=2)
abline(v=2525/1000,col="dark green",lty=2)
abline(v=6166/1000,col="dark green",lty=2)
mtext("Spawners (000s)",1,outer=T,line=1)
mtext("Recruits (000s)",2,outer=T,line=1)

dev.off()



#Lowe/Simpson/Weir


Lowe.cu <- subset(sx.rs, CU=="Lowe/Simpson/Weir")

alpha <- sx.post[,35]; beta <- sx.post[,84]

max(Lowe.cu$Esc, na.rm=T)

spw <- seq(0,max(Lowe.cu$Esc, na.rm=T),100)
spw
recruits <- array(NA,dim=c(1,length(spw),10000))

for(j in 1:10000){
  iter <- sample(length(alpha),1)
  recruits[,,j] <- spw*exp(alpha[iter]-beta[iter]* spw )
  
}

rec.m <- apply(recruits,c(1,2),median)/1000
rec.u <- apply(recruits,c(1,2),quantile,probs=0.975)/1000
rec.l <- apply(recruits,c(1,2),quantile,probs=0.025)/1000

pdf("LoweSR.pdf", width=6, height=4)

par(mfrow=c(1,1),mar=c(1.5,2.5,1,1.5),oma=c(2,2,2,1),new=F)
plot(Lowe.cu$Esc/1000,Lowe.cu$Rec/1000,bty='l',xlab="",ylab="",xaxt="n",yaxt="n")
axis(1)
axis(2,las=2)
mtext("Sockeye - Lowe/Simpson/Weir",3,cex=0.85,line=1)

points(spw/1000, rec.m,type="l" )

polygon(c(spw/1000,rev(spw/1000)),c(rec.u,rev(rec.l)),col=grey(0.9),border=NA)
points(spw/1000, rec.m,type="l",lwd=2 )
points(Lowe.cu$Esc/1000,Lowe.cu$Rec/1000)
3674.38	2560.63	5674.5	7221	5040	10739



abline(v=3674/1000,col="red",lwd=2)
abline(v=2560/1000,col="red",lty=2)
abline(v=5674/1000,col="red",lty=2)

abline(v=7221/1000,col="dark green",lwd=2)
abline(v=5040/1000,col="dark green",lty=2)
abline(v=10739/1000,col="dark green",lty=2)
mtext("Spawners (000s)",1,outer=T,line=1)
mtext("Recruits (000s)",2,outer=T,line=1)

dev.off()

#Mary Cove Creek

Mary.cu <- subset(sx.rs, CU=="Mary_Cove_Creek")

alpha <- sx.post[,36]; beta <- sx.post[,85]

max(Mary.cu$Esc, na.rm=T)

spw <- seq(0,max(Mary.cu$Esc, na.rm=T),100)
spw
recruits <- array(NA,dim=c(1,length(spw),10000))

for(j in 1:10000){
  iter <- sample(length(alpha),1)
  recruits[,,j] <- spw*exp(alpha[iter]-beta[iter]* spw )
  
}

rec.m <- apply(recruits,c(1,2),median)/1000
rec.u <- apply(recruits,c(1,2),quantile,probs=0.975)/1000
rec.l <- apply(recruits,c(1,2),quantile,probs=0.025)/1000

pdf("MaryCoveSR.pdf", width=6, height=4)

par(mfrow=c(1,1),mar=c(1.5,2.5,1,1.5),oma=c(2,2,2,1),new=F)
plot(Mary.cu$Esc/1000,Mary.cu$Rec/1000,bty='l',xlab="",ylab="",xaxt="n",yaxt="n")
axis(1)
axis(2,las=2)
mtext("Sockeye - Mary Cove Creek",3,cex=0.85,line=1)

points(spw/1000, rec.m,type="l" )

polygon(c(spw/1000,rev(spw/1000)),c(rec.u,rev(rec.l)),col=grey(0.9),border=NA)
points(spw/1000, rec.m,type="l",lwd=2 )
points(Mary.cu$Esc/1000,Mary.cu$Rec/1000)
3793.61	0	4843115.85	9023.5	4643	4959859117


abline(v=3793/1000,col="red",lwd=2)
abline(v=0/1000,col="red",lty=2)
abline(v=4843115/1000,col="red",lty=2)

abline(v=9023/1000,col="dark green",lwd=2)
abline(v=4643/1000,col="dark green",lty=2)
abline(v=4959859117/1000,col="dark green",lty=2)
mtext("Spawners (000s)",1,outer=T,line=1)
mtext("Recruits (000s)",2,outer=T,line=1)

dev.off()


#Namu

Namu.cu <- subset(sx.rs, CU=="Namu")

alpha <- sx.post[,38]; beta <- sx.post[,87]

max(Namu.cu$Esc, na.rm=T)

spw <- seq(0,max(Namu.cu$Esc, na.rm=T),100)
spw
recruits <- array(NA,dim=c(1,length(spw),10000))

for(j in 1:10000){
  iter <- sample(length(alpha),1)
  recruits[,,j] <- spw*exp(alpha[iter]-beta[iter]* spw )
  
}

rec.m <- apply(recruits,c(1,2),median)/1000
rec.u <- apply(recruits,c(1,2),quantile,probs=0.975)/1000
rec.l <- apply(recruits,c(1,2),quantile,probs=0.025)/1000

pdf("NamuSR.pdf", width=6, height=4)

par(mfrow=c(1,1),mar=c(1.5,2.5,1,1.5),oma=c(2,2,2,1),new=F)
plot(Namu.cu$Esc/1000,Namu.cu$Rec/1000,bty='l',xlab="",ylab="",xaxt="n",yaxt="n")
axis(1)
axis(2,las=2)
mtext("Sockeye - Namu",3,cex=0.85,line=1)

points(spw/1000, rec.m,type="l" )

polygon(c(spw/1000,rev(spw/1000)),c(rec.u,rev(rec.l)),col=grey(0.9),border=NA)
points(spw/1000, rec.m,type="l",lwd=2 )
points(Namu.cu$Esc/1000,Namu.cu$Rec/1000)
1479.96	0	4001068.7	3799	2300	2435537683


abline(v=1479/1000,col="red",lwd=2)
abline(v=0/1000,col="red",lty=2)
abline(v=4001068/1000,col="red",lty=2)

abline(v=3799/1000,col="dark green",lwd=2)
abline(v=2300/1000,col="dark green",lty=2)
abline(v=2435537683/1000,col="dark green",lty=2)
mtext("Spawners (000s)",1,outer=T,line=1)
mtext("Recruits (000s)",2,outer=T,line=1)

dev.off()


#Port John

Port.cu <- subset(sx.rs, CU=="Port_John")

alpha <- sx.post[,39]; beta <- sx.post[,88]

max(Port.cu$Esc, na.rm=T)

spw <- seq(0,max(Port.cu$Esc, na.rm=T),100)
spw
recruits <- array(NA,dim=c(1,length(spw),10000))

for(j in 1:10000){
  iter <- sample(length(alpha),1)
  recruits[,,j] <- spw*exp(alpha[iter]-beta[iter]* spw )
  
}

rec.m <- apply(recruits,c(1,2),median)/1000
rec.u <- apply(recruits,c(1,2),quantile,probs=0.975)/1000
rec.l <- apply(recruits,c(1,2),quantile,probs=0.025)/1000

pdf("PortJohnSR.pdf", width=6, height=4)

par(mfrow=c(1,1),mar=c(1.5,2.5,1,1.5),oma=c(2,2,2,1),new=F)
plot(Port.cu$Esc/1000,Port.cu$Rec/1000,bty='l',xlab="",ylab="",xaxt="n",yaxt="n")
axis(1)
axis(2,las=2)
mtext("Sockeye - Port John",3,cex=0.85,line=1)

points(spw/1000, rec.m,type="l" )

polygon(c(spw/1000,rev(spw/1000)),c(rec.u,rev(rec.l)),col=grey(0.9),border=NA)
points(spw/1000, rec.m,type="l",lwd=2 )
points(Port.cu$Esc/1000,Port.cu$Rec/1000)
413.24	185.2	1249.19	992.5	564	2605


abline(v=413/1000,col="red",lwd=2)
abline(v=185/1000,col="red",lty=2)
abline(v=1249/1000,col="red",lty=2)

abline(v=992/1000,col="dark green",lwd=2)
abline(v=564/1000,col="dark green",lty=2)
abline(v=2605/1000,col="dark green",lty=2)
mtext("Spawners (000s)",1,outer=T,line=1)
mtext("Recruits (000s)",2,outer=T,line=1)

dev.off()

#Roderick

Roderick.cu <- subset(sx.rs, CU=="Roderick")

alpha <- sx.post[,40]; beta <- sx.post[,89]

max(Roderick.cu$Esc, na.rm=T)

spw <- seq(0,max(Roderick.cu$Esc, na.rm=T),100)
spw
recruits <- array(NA,dim=c(1,length(spw),10000))

for(j in 1:10000){
  iter <- sample(length(alpha),1)
  recruits[,,j] <- spw*exp(alpha[iter]-beta[iter]* spw )
  
}

rec.m <- apply(recruits,c(1,2),median)/1000
rec.u <- apply(recruits,c(1,2),quantile,probs=0.975)/1000
rec.l <- apply(recruits,c(1,2),quantile,probs=0.025)/1000

pdf("RoderickSR.pdf", width=6, height=4)

par(mfrow=c(1,1),mar=c(1.5,2.5,1,1.5),oma=c(2,2,2,1),new=F)
plot(Roderick.cu$Esc/1000,Roderick.cu$Rec/1000,bty='l',xlab="",ylab="",xaxt="n",yaxt="n")
axis(1)
axis(2,las=2)
mtext("Sockeye - Roderick",3,cex=0.85,line=1)

points(spw/1000, rec.m,type="l" )

polygon(c(spw/1000,rev(spw/1000)),c(rec.u,rev(rec.l)),col=grey(0.9),border=NA)
points(spw/1000, rec.m,type="l",lwd=2 )
points(Roderick.cu$Esc/1000,Roderick.cu$Rec/1000)
604.97	343.88	20577.39	1094	679	439263.15


abline(v=605/1000,col="red",lwd=2)
abline(v=343/1000,col="red",lty=2)
abline(v=20577/1000,col="red",lty=2)

abline(v=1094/1000,col="dark green",lwd=2)
abline(v=679/1000,col="dark green",lty=2)
abline(v=439263/1000,col="dark green",lty=2)
mtext("Spawners (000s)",1,outer=T,line=1)
mtext("Recruits (000s)",2,outer=T,line=1)

dev.off()


#Tankeeah River

Tankeeah.cu <- subset(sx.rs, CU=="Tankeeah_River")

alpha <- sx.post[,41]; beta <- sx.post[,90]

max(Tankeeah.cu$Esc, na.rm=T)

spw <- seq(0,max(Tankeeah.cu$Esc, na.rm=T),100)
spw
recruits <- array(NA,dim=c(1,length(spw),10000))

for(j in 1:10000){
  iter <- sample(length(alpha),1)
  recruits[,,j] <- spw*exp(alpha[iter]-beta[iter]* spw )
  
}

rec.m <- apply(recruits,c(1,2),median)/1000
rec.u <- apply(recruits,c(1,2),quantile,probs=0.975)/1000
rec.l <- apply(recruits,c(1,2),quantile,probs=0.025)/1000

pdf("TankeeahSR.pdf", width=6, height=4)

par(mfrow=c(1,1),mar=c(1.5,2.5,1,1.5),oma=c(2,2,2,1),new=F)
plot(Tankeeah.cu$Esc/1000,Tankeeah.cu$Rec/1000,bty='l',xlab="",ylab="",xaxt="n",yaxt="n")
axis(1)
axis(2,las=2)
mtext("Sockeye - Tankeeah River",3,cex=0.85,line=1)

points(spw/1000, rec.m,type="l" )

polygon(c(spw/1000,rev(spw/1000)),c(rec.u,rev(rec.l)),col=grey(0.9),border=NA)
points(spw/1000, rec.m,type="l",lwd=2 )
points(Tankeeah.cu$Esc/1000,Tankeeah.cu$Rec/1000)
1456.71	0	2230647.06	4845	3078	824152418.8



abline(v=1456/1000,col="red",lwd=2)
abline(v=0/1000,col="red",lty=2)
abline(v=2230647/1000,col="red",lty=2)

abline(v=4845/1000,col="dark green",lwd=2)
abline(v=3078/1000,col="dark green",lty=2)
abline(v=824152418/1000,col="dark green",lty=2)
mtext("Spawners (000s)",1,outer=T,line=1)
mtext("Recruits (000s)",2,outer=T,line=1)

dev.off()

#Tsimtack/Moore/Roger

Tsim.cu <- subset(sx.rs, CU=="Tsimtack/Moore/Roger")

alpha <- sx.post[,42]; beta <- sx.post[,91]

max(Tsim.cu$Esc, na.rm=T)

spw <- seq(0,max(Tsim.cu$Esc, na.rm=T),100)
spw
recruits <- array(NA,dim=c(1,length(spw),10000))

for(j in 1:10000){
  iter <- sample(length(alpha),1)
  recruits[,,j] <- spw*exp(alpha[iter]-beta[iter]* spw )
  
}

rec.m <- apply(recruits,c(1,2),median)/1000
rec.u <- apply(recruits,c(1,2),quantile,probs=0.975)/1000
rec.l <- apply(recruits,c(1,2),quantile,probs=0.025)/1000

pdf("TsimtackMooreRogerSR.pdf", width=6, height=4)

par(mfrow=c(1,1),mar=c(1.5,2.5,1,1.5),oma=c(2,2,2,1),new=F)
plot(Tsim.cu$Esc/1000,Tsim.cu$Rec/1000,bty='l',xlab="",ylab="",xaxt="n",yaxt="n")
axis(1)
axis(2,las=2)
mtext("Sockeye - Tsimtack/Moore/Roger",3,cex=0.85,line=1)

points(spw/1000, rec.m,type="l" )

polygon(c(spw/1000,rev(spw/1000)),c(rec.u,rev(rec.l)),col=grey(0.9),border=NA)
points(spw/1000, rec.m,type="l",lwd=2 )
points(Tsim.cu$Esc/1000,Tsim.cu$Rec/1000)
976.88	493.47	2593.28	2165	1427	4949.03



abline(v=976/1000,col="red",lwd=2)
abline(v=493/1000,col="red",lty=2)
abline(v=2593/1000,col="red",lty=2)

abline(v=2165/1000,col="dark green",lwd=2)
abline(v=1427/1000,col="dark green",lty=2)
abline(v=4949/1000,col="dark green",lty=2)
mtext("Spawners (000s)",1,outer=T,line=1)
mtext("Recruits (000s)",2,outer=T,line=1)

dev.off()

#Yeo

Yeo.cu <- subset(sx.rs, CU=="Yeo")

alpha <- sx.post[,43]; beta <- sx.post[,92]

max(Yeo.cu$Esc, na.rm=T)

spw <- seq(0,max(Yeo.cu$Esc, na.rm=T),100)
spw
recruits <- array(NA,dim=c(1,length(spw),10000))

for(j in 1:10000){
  iter <- sample(length(alpha),1)
  recruits[,,j] <- spw*exp(alpha[iter]-beta[iter]* spw )
  
}

rec.m <- apply(recruits,c(1,2),median)/1000
rec.u <- apply(recruits,c(1,2),quantile,probs=0.975)/1000
rec.l <- apply(recruits,c(1,2),quantile,probs=0.025)/1000

pdf("YeoSR.pdf", width=6, height=4)

par(mfrow=c(1,1),mar=c(1.5,2.5,1,1.5),oma=c(2,2,2,1),new=F)
plot(Yeo.cu$Esc/1000,Yeo.cu$Rec/1000,bty='l',xlab="",ylab="",xaxt="n",yaxt="n")
axis(1)
axis(2,las=2)
mtext("Sockeye - Yeo",3,cex=0.85,line=1)

points(spw/1000, rec.m,type="l" )

polygon(c(spw/1000,rev(spw/1000)),c(rec.u,rev(rec.l)),col=grey(0.9),border=NA)
points(spw/1000, rec.m,type="l",lwd=2 )
points(Yeo.cu$Esc/1000,Yeo.cu$Rec/1000)
314.69	164.29	688.31	1026	760	1683


abline(v=314/1000,col="red",lwd=2)
abline(v=164/1000,col="red",lty=2)
abline(v=688/1000,col="red",lty=2)

abline(v=1026/1000,col="dark green",lwd=2)
abline(v=760/1000,col="dark green",lty=2)
abline(v=1683/1000,col="dark green",lty=2)
mtext("Spawners (000s)",1,outer=T,line=1)
mtext("Recruits (000s)",2,outer=T,line=1)

dev.off()

###############
#################
##############
#CHUM

#Smith_Inlet

Smith.cu <- subset(sx.rs, CU=="Smith_Inlet")

alpha <- cm.post[,10]; beta <- cm.post[,23]

max(Smith.cu$Esc, na.rm=T)

spw <- seq(0,max(Smith.cu$Esc, na.rm=T),100)
spw
recruits <- array(NA,dim=c(1,length(spw),10000))

for(j in 1:10000){
  iter <- sample(length(alpha),1)
  recruits[,,j] <- spw*exp(alpha[iter]-beta[iter]* spw )
  
}

rec.m <- apply(recruits,c(1,2),median)/1000
rec.u <- apply(recruits,c(1,2),quantile,probs=0.975)/1000
rec.l <- apply(recruits,c(1,2),quantile,probs=0.025)/1000

pdf("CMSmithSR.pdf", width=6, height=4)

par(mfrow=c(1,1),mar=c(1.5,2.5,1,1.5),oma=c(2,2,2,1),new=F)
plot(Smith.cu$Esc/1000,Smith.cu$Rec/1000,bty='l',xlab="",ylab="",xaxt="n",yaxt="n")
axis(1)
axis(2,las=2)
mtext("Chum - Smith Inlet",3,cex=0.85,line=1)

points(spw/1000, rec.m,type="l" )

polygon(c(spw/1000,rev(spw/1000)),c(rec.u,rev(rec.l)),col=grey(0.9),border=NA)
points(spw/1000, rec.m,type="l",lwd=2 )
points(Smith.cu$Esc/1000,Smith.cu$Rec/1000)
16057	6654.64	21031.01	32241	26655.97	42062.03


abline(v=16057/1000,col="red",lwd=2)
abline(v=6654/1000,col="red",lty=2)
abline(v=21031/1000,col="red",lty=2)

abline(v=32241/1000,col="dark green",lwd=2)
abline(v=26655/1000,col="dark green",lty=2)
abline(v=42062/1000,col="dark green",lty=2)
mtext("Spawners (000s)",1,outer=T,line=1)
mtext("Recruits (000s)",2,outer=T,line=1)

dev.off()


#Smith_Inlet

Rivers.cu <- subset(sx.rs, CU=="Rivers_Inlet")

alpha <- cm.post[,11]; beta <- cm.post[,24]

max(Rivers.cu$Esc, na.rm=T)

spw <- seq(0,max(Rivers.cu$Esc, na.rm=T),100)
spw
recruits <- array(NA,dim=c(1,length(spw),10000))

for(j in 1:10000){
  iter <- sample(length(alpha),1)
  recruits[,,j] <- spw*exp(alpha[iter]-beta[iter]* spw )
  
}


rec.m <- apply(recruits,c(1,2),median)/1000
rec.u <- apply(recruits,c(1,2),quantile,probs=0.975)/1000
rec.l <- apply(recruits,c(1,2),quantile,probs=0.025)/1000

pdf("CMRiversSR.pdf", width=6, height=4)

par(mfrow=c(1,1),mar=c(1.5,2.5,1,1.5),oma=c(2,2,2,1),new=F)
plot(Rivers.cu$Esc/1000,Rivers.cu$Rec/1000,bty='l',xlab="",ylab="",xaxt="n",yaxt="n")
axis(1)
axis(2,las=2)
mtext("Chum - Rivers Inlet",3,cex=0.85,line=1)

points(spw/1000, rec.m,type="l" )

polygon(c(spw/1000,rev(spw/1000)),c(rec.u,rev(rec.l)),col=grey(0.9),border=NA)
points(spw/1000, rec.m,type="l",lwd=2 )
points(Rivers.cu$Esc/1000,Rivers.cu$Rec/1000)
13717.5	6167.69	19205	28029	21777.92	38415.03



abline(v=13717/1000,col="red",lwd=2)
abline(v=6167/1000,col="red",lty=2)
abline(v=19205/1000,col="red",lty=2)

abline(v=28029/1000,col="dark green",lwd=2)
abline(v=21777/1000,col="dark green",lty=2)
abline(v=38415/1000,col="dark green",lty=2)
mtext("Spawners (000s)",1,outer=T,line=1)
mtext("Recruits (000s)",2,outer=T,line=1)

dev.off()



#Spiller-Fitz-Hugh-Burke


Spiller.cu <- subset(sx.rs, CU=="Spiller-Fitz-Hugh-Burke")

alpha <- cm.post[,12]; beta <- cm.post[,25]

max(Spiller.cu$Esc, na.rm=T)

spw <- seq(0,max(Spiller.cu$Esc, na.rm=T),100)
spw
recruits <- array(NA,dim=c(1,length(spw),10000))

for(j in 1:10000){
  iter <- sample(length(alpha),1)
  recruits[,,j] <- spw*exp(alpha[iter]-beta[iter]* spw )
  
}


rec.m <- apply(recruits,c(1,2),median)/1000
rec.u <- apply(recruits,c(1,2),quantile,probs=0.975)/1000
rec.l <- apply(recruits,c(1,2),quantile,probs=0.025)/1000

pdf("CMSpillerSR.pdf", width=6, height=4)

par(mfrow=c(1,1),mar=c(1.5,2.5,1,1.5),oma=c(2,2,2,1),new=F)
plot(Spiller.cu$Esc/1000,Spiller.cu$Rec/1000,bty='l',xlab="",ylab="",xaxt="n",yaxt="n")
axis(1)
axis(2,las=2)
mtext("Chum - Spiller-Fitz-Hugh-Burke",3,cex=0.85,line=1)

points(spw/1000, rec.m,type="l" )

polygon(c(spw/1000,rev(spw/1000)),c(rec.u,rev(rec.l)),col=grey(0.9),border=NA)
points(spw/1000, rec.m,type="l",lwd=2 )
points(Spiller.cu$Esc/1000,Spiller.cu$Rec/1000)
88361.25	72710.46	122025.36	176722.5	145420.92	244050.73


abline(v=88361/1000,col="red",lwd=2)
abline(v=72710/1000,col="red",lty=2)
abline(v=122025/1000,col="red",lty=2)

abline(v=176722/1000,col="dark green",lwd=2)
abline(v=145420/1000,col="dark green",lty=2)
abline(v=244050/1000,col="dark green",lty=2)
mtext("Spawners (000s)",1,outer=T,line=1)
mtext("Recruits (000s)",2,outer=T,line=1)

dev.off()



#Bella Colla-Dean Rivers


BCDean.cu <- subset(sx.rs, CU=="Bella_Colla-Dean_Rivers")

alpha <- cm.post[,13]; beta <- cm.post[,26]

max(BCDean.cu$Esc, na.rm=T)

spw <- seq(0,max(BCDean.cu$Esc, na.rm=T),100)
spw
recruits <- array(NA,dim=c(1,length(spw),10000))

for(j in 1:10000){
  iter <- sample(length(alpha),1)
  recruits[,,j] <- spw*exp(alpha[iter]-beta[iter]* spw )
  
}


rec.m <- apply(recruits,c(1,2),median)/1000
rec.u <- apply(recruits,c(1,2),quantile,probs=0.975)/1000
rec.l <- apply(recruits,c(1,2),quantile,probs=0.025)/1000

pdf("CMBellaCoolaDeanSR.pdf", width=6, height=4)

par(mfrow=c(1,1),mar=c(1.5,2.5,1,1.5),oma=c(2,2,2,1),new=F)
plot(BCDean.cu$Esc/1000,BCDean.cu$Rec/1000,bty='l',xlab="",ylab="",xaxt="n",yaxt="n")
axis(1)
axis(2,las=2)
mtext("Chum - Bella Coola - Dean",3,cex=0.85,line=1)

points(spw/1000, rec.m,type="l" )

polygon(c(spw/1000,rev(spw/1000)),c(rec.u,rev(rec.l)),col=grey(0.9),border=NA)
points(spw/1000, rec.m,type="l",lwd=2 )
points(BCDean.cu$Esc/1000,BCDean.cu$Rec/1000)
44521.5	35034.47	68421.01	89050.5	70092.98	136969.15



abline(v=44521/1000,col="red",lwd=2)
abline(v=35034/1000,col="red",lty=2)
abline(v=68421/1000,col="red",lty=2)

abline(v=89050/1000,col="dark green",lwd=2)
abline(v=70092/1000,col="dark green",lty=2)
abline(v=136969/1000,col="dark green",lty=2)
mtext("Spawners (000s)",1,outer=T,line=1)
mtext("Recruits (000s)",2,outer=T,line=1)

dev.off()


#Bella Coola Late


BCLate.cu <- subset(sx.rs, CU=="Bella_Coola_River-Late")

alpha <- cm.post[,2]; beta <- cm.post[,15]

max(BCLate.cu$Esc, na.rm=T)

spw <- seq(0,max(BCLate.cu$Esc, na.rm=T),100)
spw
recruits <- array(NA,dim=c(1,length(spw),10000))

for(j in 1:10000){
  iter <- sample(length(alpha),1)
  recruits[,,j] <- spw*exp(alpha[iter]-beta[iter]* spw )
  
}


rec.m <- apply(recruits,c(1,2),median)/1000
rec.u <- apply(recruits,c(1,2),quantile,probs=0.975)/1000
rec.l <- apply(recruits,c(1,2),quantile,probs=0.025)/1000

jpeg("CMBellaCoolaLateSR1.jpg", width=6, height=4, units="in", res=200)

par(mfrow=c(1,1),mar=c(1.5,2.5,1,1.5),oma=c(2,2,2,1),new=F)
plot(BCLate.cu$Esc/1000,BCLate.cu$Rec/1000,bty='l',xlab="",ylab="",xaxt="n",yaxt="n")
axis(1)
axis(2,las=2)
mtext("Chum - Bella Coola - Late",3,cex=0.85,line=1)

#points(spw/1000, rec.m,type="l" )

#polygon(c(spw/1000,rev(spw/1000)),c(rec.u,rev(rec.l)),col=grey(0.9),border=NA)
#points(spw/1000, rec.m,type="l",lwd=2 )
#points(BCLate.cu$Esc/1000,BCLate.cu$Rec/1000)
59761	49163.97	79898.04	119522	98327.95	159796.08
#abline(a=0,b=1,lty=2,lwd=2)

#abline(v=59761/1000,col="red",lwd=2)
#abline(v=49163/1000,col="red",lty=2)
#abline(v=79898/1000,col="red",lty=2)

#abline(v=119522/1000,col="dark green",lwd=2)
#abline(v=98327/1000,col="dark green",lty=2)
#abline(v=159796/1000,col="dark green",lty=2)
mtext("Spawners (000s)",1,outer=T,line=1)
mtext("Recruits (000s)",2,outer=T,line=1)

dev.off()



#Hecate low


HecLow.cu <- subset(sx.rs, CU=="Hecate_Lowlands")

alpha <- cm.post[,3]; beta <- cm.post[,16]

max(HecLow.cu$Esc, na.rm=T)

spw <- seq(0,max(HecLow.cu$Esc, na.rm=T),10000)
spw
recruits <- array(NA,dim=c(1,length(spw),10000))

for(j in 1:10000){
  iter <- sample(length(alpha),1)
  recruits[,,j] <- spw*exp(alpha[iter]-beta[iter]* spw )
  
}


rec.m <- apply(recruits,c(1,2),median)/1000
rec.u <- apply(recruits,c(1,2),quantile,probs=0.975)/1000
rec.l <- apply(recruits,c(1,2),quantile,probs=0.025)/1000

pdf("CMHecateLowSR.pdf", width=6, height=4)

par(mfrow=c(1,1),mar=c(1.5,2.5,1,1.5),oma=c(2,2,2,1),new=F)
plot(HecLow.cu$Esc/1000,HecLow.cu$Rec/1000,bty='l',xlab="",ylab="",xaxt="n",yaxt="n")
axis(1)
axis(2,las=2)
mtext("Chum - Hecate Lowlands",3,cex=0.85,line=1)

points(spw/1000, rec.m,type="l" )

polygon(c(spw/1000,rev(spw/1000)),c(rec.u,rev(rec.l)),col=grey(0.9),border=NA)
points(spw/1000, rec.m,type="l",lwd=2 )
points(HecLow.cu$Esc/1000,HecLow.cu$Rec/1000)
41375.5	34786	54000.08	82751	69572	108000.15


abline(v=41375/1000,col="red",lwd=2)
abline(v=34786/1000,col="red",lty=2)
abline(v=54000/1000,col="red",lty=2)

abline(v=82751/1000,col="dark green",lwd=2)
abline(v=69572/1000,col="dark green",lty=2)
abline(v=108000/1000,col="dark green",lty=2)
mtext("Spawners (000s)",1,outer=T,line=1)
mtext("Recruits (000s)",2,outer=T,line=1)

dev.off()

#Mussel-Kynoch


Mussel.cu <- subset(sx.rs, CU=="Mussel-Kynock")

alpha <- cm.post[,4]; beta <- cm.post[,17]

max(Mussel.cu$Esc, na.rm=T)

spw <- seq(0,max(HecLow.cu$Esc, na.rm=T),10000)
spw
recruits <- array(NA,dim=c(1,length(spw),10000))

for(j in 1:10000){
  iter <- sample(length(alpha),1)
  recruits[,,j] <- spw*exp(alpha[iter]-beta[iter]* spw )
  
}


rec.m <- apply(recruits,c(1,2),median)/1000
rec.u <- apply(recruits,c(1,2),quantile,probs=0.975)/1000
rec.l <- apply(recruits,c(1,2),quantile,probs=0.025)/1000

pdf("CMMusselKynockSR.pdf", width=6, height=4)

par(mfrow=c(1,1),mar=c(1.5,2.5,1,1.5),oma=c(2,2,2,1),new=F)
plot(Mussel.cu$Esc/1000,Mussel.cu$Rec/1000,bty='l',xlab="",ylab="",xaxt="n",yaxt="n")
axis(1)
axis(2,las=2)
mtext("Chum - Mussel-Kynock",3,cex=0.85,line=1)

points(spw/1000, rec.m,type="l" )

polygon(c(spw/1000,rev(spw/1000)),c(rec.u,rev(rec.l)),col=grey(0.9),border=NA)
points(spw/1000, rec.m,type="l",lwd=2 )
points(Mussel.cu$Esc/1000,Mussel.cu$Rec/1000)
32366.75	28037.99	39412.51	64733.5	56075.98	78825.03



abline(v=32366/1000,col="red",lwd=2)
abline(v=28037/1000,col="red",lty=2)
abline(v=39412/1000,col="red",lty=2)

abline(v=64733/1000,col="dark green",lwd=2)
abline(v=56075/1000,col="dark green",lty=2)
abline(v=78825/1000,col="dark green",lty=2)
mtext("Spawners (000s)",1,outer=T,line=1)
mtext("Recruits (000s)",2,outer=T,line=1)

dev.off()


#Douglas-Gardner

Douglas.cu <- subset(sx.rs, CU=="Douglas-Gardner")

alpha <- cm.post[,5]; beta <- cm.post[,18]

max(Douglas.cu$Esc, na.rm=T)

spw <- seq(0,max(Douglas.cu$Esc, na.rm=T),10000)
spw
recruits <- array(NA,dim=c(1,length(spw),10000))

for(j in 1:10000){
  iter <- sample(length(alpha),1)
  recruits[,,j] <- spw*exp(alpha[iter]-beta[iter]* spw )
  
}


rec.m <- apply(recruits,c(1,2),median)/1000
rec.u <- apply(recruits,c(1,2),quantile,probs=0.975)/1000
rec.l <- apply(recruits,c(1,2),quantile,probs=0.025)/1000

pdf("CMDouglasGardnerSR.pdf", width=6, height=4)

par(mfrow=c(1,1),mar=c(1.5,2.5,1,1.5),oma=c(2,2,2,1),new=F)
plot(Douglas.cu$Esc/1000,Douglas.cu$Rec/1000,bty='l',xlab="",ylab="",xaxt="n",yaxt="n")
axis(1)
axis(2,las=2)
mtext("Chum - Douglas-Gardner",3,cex=0.85,line=1)

points(spw/1000, rec.m,type="l" )

polygon(c(spw/1000,rev(spw/1000)),c(rec.u,rev(rec.l)),col=grey(0.9),border=NA)
points(spw/1000, rec.m,type="l",lwd=2 )
points(Douglas.cu$Esc/1000,Douglas.cu$Rec/1000)
91657	72728.44	129551.06	183315	145459.93	259107.53


abline(v=91657/1000,col="red",lwd=2)
abline(v=72728/1000,col="red",lty=2)
abline(v=129551/1000,col="red",lty=2)

abline(v=183315/1000,col="dark green",lwd=2)
abline(v=145459/1000,col="dark green",lty=2)
abline(v=259107/1000,col="dark green",lty=2)
mtext("Spawners (000s)",1,outer=T,line=1)
mtext("Recruits (000s)",2,outer=T,line=1)

dev.off()


##########
###########
###########
#Chinook



#Rivers Inlet

Rivers.cu <- subset(sx.rs, CU=="Rivers_Inlet")

alpha <- cn.post[,3]; beta <- cn.post[,19]

max(Rivers.cu$Esc, na.rm=T)

spw <- seq(0,max(Rivers.cu$Esc, na.rm=T),100)
spw
recruits <- array(NA,dim=c(1,length(spw),10000))

for(j in 1:10000){
  iter <- sample(length(alpha),1)
  recruits[,,j] <- spw*exp(alpha[iter]-beta[iter]* spw )
  
}


rec.m <- apply(recruits,c(1,2),median)/1000
rec.u <- apply(recruits,c(1,2),quantile,probs=0.975)/1000
rec.l <- apply(recruits,c(1,2),quantile,probs=0.025)/1000

pdf("CNRiversInletSR.pdf", width=6, height=4)

par(mfrow=c(1,1),mar=c(1.5,2.5,1,1.5),oma=c(2,2,2,1),new=F)
plot(Rivers.cu$Esc/1000,Rivers.cu$Rec/1000,bty='l',xlab="",ylab="",xaxt="n",yaxt="n")
axis(1)
axis(2,las=2)
mtext("Chinook - Rivers Inlet",3,cex=0.85,line=1)

points(spw/1000, rec.m,type="l" )

polygon(c(spw/1000,rev(spw/1000)),c(rec.u,rev(rec.l)),col=grey(0.9),border=NA)
points(spw/1000, rec.m,type="l",lwd=2 )
points(Rivers.cu$Esc/1000,Rivers.cu$Rec/1000)
464.43	271.75	844.97	1118	848	1656



abline(v=464/1000,col="red",lwd=2)
abline(v=271/1000,col="red",lty=2)
abline(v=844/1000,col="red",lty=2)

abline(v=1118/1000,col="dark green",lwd=2)
abline(v=848/1000,col="dark green",lty=2)
abline(v=1656/1000,col="dark green",lty=2)
mtext("Spawners (000s)",1,outer=T,line=1)
mtext("Recruits (000s)",2,outer=T,line=1)

dev.off()

#Wannock

Wannock.cu <- subset(sx.rs, CU=="Wannock")

alpha <- cn.post[,4]; beta <- cn.post[,20]

max(Wannock.cu$Esc, na.rm=T)

spw <- seq(0,max(Wannock.cu$Esc, na.rm=T),100)
spw
recruits <- array(NA,dim=c(1,length(spw),10000))

for(j in 1:10000){
  iter <- sample(length(alpha),1)
  recruits[,,j] <- spw*exp(alpha[iter]-beta[iter]* spw )
  
}


rec.m <- apply(recruits,c(1,2),median)/1000
rec.u <- apply(recruits,c(1,2),quantile,probs=0.975)/1000
rec.l <- apply(recruits,c(1,2),quantile,probs=0.025)/1000

pdf("CNWannockSR.pdf", width=6, height=4)

par(mfrow=c(1,1),mar=c(1.5,2.5,1,1.5),oma=c(2,2,2,1),new=F)
plot(Wannock.cu$Esc/1000,Wannock.cu$Rec/1000,bty='l',xlab="",ylab="",xaxt="n",yaxt="n")
axis(1)
axis(2,las=2)
mtext("Chinook - Wannock",3,cex=0.85,line=1)

points(spw/1000, rec.m,type="l" )

polygon(c(spw/1000,rev(spw/1000)),c(rec.u,rev(rec.l)),col=grey(0.9),border=NA)
points(spw/1000, rec.m,type="l",lwd=2 )
points(Wannock.cu$Esc/1000,Wannock.cu$Rec/1000)
817.71	457.56	1549.75	3830	3225.98	4910


abline(v=817/1000,col="red",lwd=2)
abline(v=457/1000,col="red",lty=2)
abline(v=1549/1000,col="red",lty=2)

abline(v=3830/1000,col="dark green",lwd=2)
abline(v=3225/1000,col="dark green",lty=2)
abline(v=4910/1000,col="dark green",lty=2)
mtext("Spawners (000s)",1,outer=T,line=1)
mtext("Recruits (000s)",2,outer=T,line=1)

dev.off()

#BC Bentick

BC.cu <- subset(sx.rs, CU=="Bella_Coola-Bentinck")

alpha <- cn.post[,5]; beta <- cn.post[,21]

max(BC.cu$Esc, na.rm=T)

spw <- seq(0,max(BC.cu$Esc, na.rm=T),100)
spw
recruits <- array(NA,dim=c(1,length(spw),10000))

for(j in 1:10000){
  iter <- sample(length(alpha),1)
  recruits[,,j] <- spw*exp(alpha[iter]-beta[iter]* spw )
  
}


rec.m <- apply(recruits,c(1,2),median)/1000
rec.u <- apply(recruits,c(1,2),quantile,probs=0.975)/1000
rec.l <- apply(recruits,c(1,2),quantile,probs=0.025)/1000

pdf("CNBCBentickSR.pdf", width=6, height=4)

par(mfrow=c(1,1),mar=c(1.5,2.5,1,1.5),oma=c(2,2,2,1),new=F)
plot(BC.cu$Esc/1000,BC.cu$Rec/1000,bty='l',xlab="",ylab="",xaxt="n",yaxt="n")
axis(1)
axis(2,las=2)
mtext("Chinook - Bella Coola-Bentinck",3,cex=0.85,line=1)

points(spw/1000, rec.m,type="l" )

polygon(c(spw/1000,rev(spw/1000)),c(rec.u,rev(rec.l)),col=grey(0.9),border=NA)
points(spw/1000, rec.m,type="l",lwd=2 )
points(BC.cu$Esc/1000,BC.cu$Rec/1000)
5033.14	1971.76	15015.03	15188	10945.98	30030.05



abline(v=5033/1000,col="red",lwd=2)
abline(v=1971/1000,col="red",lty=2)
abline(v=15015/1000,col="red",lty=2)

abline(v=15188/1000,col="dark green",lwd=2)
abline(v=10945/1000,col="dark green",lty=2)
abline(v=30030/1000,col="dark green",lty=2)
mtext("Spawners (000s)",1,outer=T,line=1)
mtext("Recruits (000s)",2,outer=T,line=1)

dev.off()




#BC Bentick

Dean.cu <- subset(sx.rs, CU=="Dean_River")

alpha <- cn.post[,6]; beta <- cn.post[,22]

max(Dean.cu$Esc, na.rm=T)

spw <- seq(0,max(Dean.cu$Esc, na.rm=T),100)
spw
recruits <- array(NA,dim=c(1,length(spw),10000))

for(j in 1:10000){
  iter <- sample(length(alpha),1)
  recruits[,,j] <- spw*exp(alpha[iter]-beta[iter]* spw )
  
}


rec.m <- apply(recruits,c(1,2),median)/1000
rec.u <- apply(recruits,c(1,2),quantile,probs=0.975)/1000
rec.l <- apply(recruits,c(1,2),quantile,probs=0.025)/1000

pdf("CNDeanSR.pdf", width=6, height=4)

par(mfrow=c(1,1),mar=c(1.5,2.5,1,1.5),oma=c(2,2,2,1),new=F)
plot(Dean.cu$Esc/1000,Dean.cu$Rec/1000,bty='l',xlab="",ylab="",xaxt="n",yaxt="n")
axis(1)
axis(2,las=2)
mtext("Chinook - Dean River",3,cex=0.85,line=1)

points(spw/1000, rec.m,type="l" )

polygon(c(spw/1000,rev(spw/1000)),c(rec.u,rev(rec.l)),col=grey(0.9),border=NA)
points(spw/1000, rec.m,type="l",lwd=2 )
points(Dean.cu$Esc/1000,Dean.cu$Rec/1000)
350.42	183.97	661.57	1083	913	1400


abline(v=350/1000,col="red",lwd=2)
abline(v=184/1000,col="red",lty=2)
abline(v=661/1000,col="red",lty=2)

abline(v=1083/1000,col="dark green",lwd=2)
abline(v=913/1000,col="dark green",lty=2)
abline(v=1400/1000,col="dark green",lty=2)
mtext("Spawners (000s)",1,outer=T,line=1)
mtext("Recruits (000s)",2,outer=T,line=1)

dev.off()



#NCC late

NCClate.cu <- subset(sx.rs, CU=="NCC-late_timing")

alpha <- cn.post[,7]; beta <- cn.post[,23]

max(NCClate.cu$Esc, na.rm=T)

spw <- seq(0,max(NCClate.cu$Esc, na.rm=T),100)
spw
recruits <- array(NA,dim=c(1,length(spw),10000))

for(j in 1:10000){
  iter <- sample(length(alpha),1)
  recruits[,,j] <- spw*exp(alpha[iter]-beta[iter]* spw )
  
}


rec.m <- apply(recruits,c(1,2),median)/1000
rec.u <- apply(recruits,c(1,2),quantile,probs=0.975)/1000
rec.l <- apply(recruits,c(1,2),quantile,probs=0.025)/1000

pdf("CNNCClateSR.pdf", width=6, height=4)

par(mfrow=c(1,1),mar=c(1.5,2.5,1,1.5),oma=c(2,2,2,1),new=F)
plot(NCClate.cu$Esc/1000,NCClate.cu$Rec/1000,bty='l',xlab="",ylab="",xaxt="n",yaxt="n")
axis(1)
axis(2,las=2)
mtext("Chinook - NCC late timing",3,cex=0.85,line=1)

points(spw/1000, rec.m,type="l" )

polygon(c(spw/1000,rev(spw/1000)),c(rec.u,rev(rec.l)),col=grey(0.9),border=NA)
points(spw/1000, rec.m,type="l",lwd=2 )
points(NCClate.cu$Esc/1000,NCClate.cu$Rec/1000)
328.95	74.66	1964.21	1154	510



abline(v=329/1000,col="red",lwd=2)
abline(v=75/1000,col="red",lty=2)
abline(v=1964/1000,col="red",lty=2)

abline(v=1154/1000,col="dark green",lwd=2)
abline(v=510/1000,col="dark green",lty=2)
abline(v=5612/1000,col="dark green",lty=2)
mtext("Spawners (000s)",1,outer=T,line=1)
mtext("Recruits (000s)",2,outer=T,line=1)

dev.off()


#NCC late

NCCearly.cu <- subset(sx.rs, CU=="NCC-early_timing")

alpha <- cn.post[,8]; beta <- cn.post[,24]

max(NCCearly.cu$Esc, na.rm=T)

spw <- seq(0,max(NCCearly.cu$Esc, na.rm=T),100)
spw
recruits <- array(NA,dim=c(1,length(spw),10000))

for(j in 1:10000){
  iter <- sample(length(alpha),1)
  recruits[,,j] <- spw*exp(alpha[iter]-beta[iter]* spw )
  
}


rec.m <- apply(recruits,c(1,2),median)/1000
rec.u <- apply(recruits,c(1,2),quantile,probs=0.975)/1000
rec.l <- apply(recruits,c(1,2),quantile,probs=0.025)/1000

pdf("CNNCCearlySR.pdf", width=6, height=4)

par(mfrow=c(1,1),mar=c(1.5,2.5,1,1.5),oma=c(2,2,2,1),new=F)
plot(NCCearly.cu$Esc/1000,NCCearly.cu$Rec/1000,bty='l',xlab="",ylab="",xaxt="n",yaxt="n")
axis(1)
axis(2,las=2)
mtext("Chinook - NCC early timing",3,cex=0.85,line=1)

points(spw/1000, rec.m,type="l" )

polygon(c(spw/1000,rev(spw/1000)),c(rec.u,rev(rec.l)),col=grey(0.9),border=NA)
points(spw/1000, rec.m,type="l",lwd=2 )
points(NCCearly.cu$Esc/1000,NCCearly.cu$Rec/1000)
161.11	104.43	247.72	428	355	536


abline(v=161/1000,col="red",lwd=2)
abline(v=104/1000,col="red",lty=2)
abline(v=247/1000,col="red",lty=2)

abline(v=428/1000,col="dark green",lwd=2)
abline(v=355/1000,col="dark green",lty=2)
abline(v=536/1000,col="dark green",lty=2)
mtext("Spawners (000s)",1,outer=T,line=1)
mtext("Recruits (000s)",2,outer=T,line=1)

dev.off()



########### coho
#######
##########
##########
##########


#Smith

Smith.cu <- subset(sx.rs, CU=="Smith_Inlet")

alpha <- co.post[,12]; beta <- co.post[,26]

max(Smith.cu$Esc, na.rm=T)

spw <- seq(0,max(Smith.cu$Esc, na.rm=T),100)
spw
recruits <- array(NA,dim=c(1,length(spw),10000))

for(j in 1:10000){
  iter <- sample(length(alpha),1)
  recruits[,,j] <- spw*exp(alpha[iter]-beta[iter]* spw )
  
}


rec.m <- apply(recruits,c(1,2),median)/1000
rec.u <- apply(recruits,c(1,2),quantile,probs=0.975)/1000
rec.l <- apply(recruits,c(1,2),quantile,probs=0.025)/1000

pdf("COSmithSR.pdf", width=6, height=4)

par(mfrow=c(1,1),mar=c(1.5,2.5,1,1.5),oma=c(2,2,2,1),new=F)
plot(Smith.cu$Esc/1000,Smith.cu$Rec/1000,bty='l',xlab="",ylab="",xaxt="n",yaxt="n")
axis(1)
axis(2,las=2)
mtext("Coho - Smith Inlet",3,cex=0.85,line=1)

points(spw/1000, rec.m,type="l" )

polygon(c(spw/1000,rev(spw/1000)),c(rec.u,rev(rec.l)),col=grey(0.9),border=NA)
points(spw/1000, rec.m,type="l",lwd=2 )
points(Smith.cu$Esc/1000,Smith.cu$Rec/1000)
4185	0	22523731.76	8008	3634.97	27303428386



abline(v=4185/1000,col="red",lwd=2)
abline(v=0/1000,col="red",lty=2)
abline(v=22523731.76/1000,col="red",lty=2)

abline(v=8008/1000,col="dark green",lwd=2)
abline(v=3634/1000,col="dark green",lty=2)
abline(v=27303428386/1000,col="dark green",lty=2)
mtext("Spawners (000s)",1,outer=T,line=1)
mtext("Recruits (000s)",2,outer=T,line=1)

dev.off()


#Rivers

Rivers.cu <- subset(sx.rs, CU=="Rivers_Inlet")

alpha <- co.post[,13]; beta <- co.post[,27]

max(Rivers.cu$Esc, na.rm=T)

spw <- seq(0,max(Rivers.cu$Esc, na.rm=T),100)
spw
recruits <- array(NA,dim=c(1,length(spw),10000))

for(j in 1:10000){
  iter <- sample(length(alpha),1)
  recruits[,,j] <- spw*exp(alpha[iter]-beta[iter]* spw )
  
}


rec.m <- apply(recruits,c(1,2),median)/1000
rec.u <- apply(recruits,c(1,2),quantile,probs=0.975)/1000
rec.l <- apply(recruits,c(1,2),quantile,probs=0.025)/1000

pdf("CORiversSR.pdf", width=6, height=4)

par(mfrow=c(1,1),mar=c(1.5,2.5,1,1.5),oma=c(2,2,2,1),new=F)
plot(Rivers.cu$Esc/1000,Rivers.cu$Rec/1000,bty='l',xlab="",ylab="",xaxt="n",yaxt="n")
axis(1)
axis(2,las=2)
mtext("Coho - Rivers Inlet",3,cex=0.85,line=1)

points(spw/1000, rec.m,type="l" )

polygon(c(spw/1000,rev(spw/1000)),c(rec.u,rev(rec.l)),col=grey(0.9),border=NA)
points(spw/1000, rec.m,type="l",lwd=2 )
points(Rivers.cu$Esc/1000,Rivers.cu$Rec/1000)
23916	0	3076644.39	49060	28736.92	1783533110

abline(v=23916/1000,col="red",lwd=2)
abline(v=0/1000,col="red",lty=2)
abline(v=3076644/1000,col="red",lty=2)

abline(v=49060/1000,col="dark green",lwd=2)
abline(v=28736/1000,col="dark green",lty=2)
abline(v=1783533110/1000,col="dark green",lty=2)
mtext("Spawners (000s)",1,outer=T,line=1)
mtext("Recruits (000s)",2,outer=T,line=1)

dev.off()



#Bella Coola-Dean Rivers

BCDean.cu <- subset(sx.rs, CU=="Bella_Coola-Dean_Rivers")

alpha <- co.post[,14]; beta <- co.post[,28]

max(BCDean.cu$Esc, na.rm=T)

spw <- seq(0,max(BCDean.cu$Esc, na.rm=T),100)
spw
recruits <- array(NA,dim=c(1,length(spw),10000))

for(j in 1:10000){
  iter <- sample(length(alpha),1)
  recruits[,,j] <- spw*exp(alpha[iter]-beta[iter]* spw )
  
}


rec.m <- apply(recruits,c(1,2),median)/1000
rec.u <- apply(recruits,c(1,2),quantile,probs=0.975)/1000
rec.l <- apply(recruits,c(1,2),quantile,probs=0.025)/1000

pdf("COBellaCoolaDeanSR.pdf", width=6, height=4)

par(mfrow=c(1,1),mar=c(1.5,2.5,1,1.5),oma=c(2,2,2,1),new=F)
plot(BCDean.cu$Esc/1000,BCDean.cu$Rec/1000,bty='l',xlab="",ylab="",xaxt="n",yaxt="n")
axis(1)
axis(2,las=2)
mtext("Coho - Bella Coola-Dean Rivers",3,cex=0.85,line=1)

points(spw/1000, rec.m,type="l" )

polygon(c(spw/1000,rev(spw/1000)),c(rec.u,rev(rec.l)),col=grey(0.9),border=NA)
points(spw/1000, rec.m,type="l",lwd=2 )
points(BCDean.cu$Esc/1000,BCDean.cu$Rec/1000)
13337.5	6339.29	18285.01	26747	21789	36570.03

abline(v=13337/1000,col="red",lwd=2)
abline(v=6339/1000,col="red",lty=2)
abline(v=18285/1000,col="red",lty=2)

abline(v=26747/1000,col="dark green",lwd=2)
abline(v=21789/1000,col="dark green",lty=2)
abline(v=36570/1000,col="dark green",lty=2)
mtext("Spawners (000s)",1,outer=T,line=1)
mtext("Recruits (000s)",2,outer=T,line=1)

dev.off()



#Bella Coola-Dean Rivers

Mussel.cu <- subset(sx.rs, CU=="Mussel-Kynoch")

alpha <- co.post[,2]; beta <- co.post[,16]

max(Mussel.cu$Esc, na.rm=T)

spw <- seq(0,max(Mussel.cu$Esc, na.rm=T),100)
spw
recruits <- array(NA,dim=c(1,length(spw),10000))

for(j in 1:10000){
  iter <- sample(length(alpha),1)
  recruits[,,j] <- spw*exp(alpha[iter]-beta[iter]* spw )
  
}


rec.m <- apply(recruits,c(1,2),median)/1000
rec.u <- apply(recruits,c(1,2),quantile,probs=0.975)/1000
rec.l <- apply(recruits,c(1,2),quantile,probs=0.025)/1000

pdf("COMusselKynockSR.pdf", width=6, height=4)

par(mfrow=c(1,1),mar=c(1.5,2.5,1,1.5),oma=c(2,2,2,1),new=F)
plot(Mussel.cu$Esc/1000,Mussel.cu$Rec/1000,bty='l',xlab="",ylab="",xaxt="n",yaxt="n")
axis(1)
axis(2,las=2)
mtext("Coho - Mussel-Kynoch",3,cex=0.85,line=1)

points(spw/1000, rec.m,type="l" )

polygon(c(spw/1000,rev(spw/1000)),c(rec.u,rev(rec.l)),col=grey(0.9),border=NA)
points(spw/1000, rec.m,type="l",lwd=2 )
points(Mussel.cu$Esc/1000,Mussel.cu$Rec/1000)
3112	1884.95	7667.32	6154	4397.97	13882.1


abline(v=3112/1000,col="red",lwd=2)
abline(v=1884/1000,col="red",lty=2)
abline(v=7667/1000,col="red",lty=2)

abline(v=6154/1000,col="dark green",lwd=2)
abline(v=4397/1000,col="dark green",lty=2)
abline(v=13882/1000,col="dark green",lty=2)
mtext("Spawners (000s)",1,outer=T,line=1)
mtext("Recruits (000s)",2,outer=T,line=1)

dev.off()



#Hecate Strait Mainland

Hecate.cu <- subset(sx.rs, CU=="Hecate_Strait_Mainland")

alpha <- co.post[,3]; beta <- co.post[,17]

max(Hecate.cu$Esc, na.rm=T)

spw <- seq(0,max(Hecate.cu$Esc, na.rm=T),100)
spw
recruits <- array(NA,dim=c(1,length(spw),10000))

for(j in 1:10000){
  iter <- sample(length(alpha),1)
  recruits[,,j] <- spw*exp(alpha[iter]-beta[iter]* spw )
  
}


rec.m <- apply(recruits,c(1,2),median)/1000
rec.u <- apply(recruits,c(1,2),quantile,probs=0.975)/1000
rec.l <- apply(recruits,c(1,2),quantile,probs=0.025)/1000

pdf("COHecateSR.pdf", width=6, height=4)

par(mfrow=c(1,1),mar=c(1.5,2.5,1,1.5),oma=c(2,2,2,1),new=F)
plot(Hecate.cu$Esc/1000,Hecate.cu$Rec/1000,bty='l',xlab="",ylab="",xaxt="n",yaxt="n")
axis(1)
axis(2,las=2)
mtext("Coho - Hecate Strait Mainland",3,cex=0.85,line=1)

points(spw/1000, rec.m,type="l" )

polygon(c(spw/1000,rev(spw/1000)),c(rec.u,rev(rec.l)),col=grey(0.9),border=NA)
points(spw/1000, rec.m,type="l",lwd=2 )
points(Hecate.cu$Esc/1000,Hecate.cu$Rec/1000)
67715.25	51290.48	115570.68	135502.5	102842.98	234091.15

abline(v=67715/1000,col="red",lwd=2)
abline(v=51290/1000,col="red",lty=2)
abline(v=115570/1000,col="red",lty=2)

abline(v=135502/1000,col="dark green",lwd=2)
abline(v=102842/1000,col="dark green",lty=2)
abline(v=234091/1000,col="dark green",lty=2)
mtext("Spawners (000s)",1,outer=T,line=1)
mtext("Recruits (000s)",2,outer=T,line=1)

dev.off()


#Brim-Wahoo

Brim.cu <- subset(sx.rs, CU=="Brim-Wahoo")

alpha <- co.post[,4]; beta <- co.post[,18]

max(Brim.cu$Esc, na.rm=T)

spw <- seq(0,max(Brim.cu$Esc, na.rm=T),100)
spw
recruits <- array(NA,dim=c(1,length(spw),10000))

for(j in 1:10000){
  iter <- sample(length(alpha),1)
  recruits[,,j] <- spw*exp(alpha[iter]-beta[iter]* spw )
  
}


rec.m <- apply(recruits,c(1,2),median)/1000
rec.u <- apply(recruits,c(1,2),quantile,probs=0.975)/1000
rec.l <- apply(recruits,c(1,2),quantile,probs=0.025)/1000

pdf("COBrimWahooSR.pdf", width=6, height=4)

par(mfrow=c(1,1),mar=c(1.5,2.5,1,1.5),oma=c(2,2,2,1),new=F)
plot(Brim.cu$Esc/1000,Brim.cu$Rec/1000,bty='l',xlab="",ylab="",xaxt="n",yaxt="n")
axis(1)
axis(2,las=2)
mtext("Coho - Brim-Wahoo",3,cex=0.85,line=1)

points(spw/1000, rec.m,type="l" )

polygon(c(spw/1000,rev(spw/1000)),c(rec.u,rev(rec.l)),col=grey(0.9),border=NA)
points(spw/1000, rec.m,type="l",lwd=2 )
points(Brim.cu$Esc/1000,Brim.cu$Rec/1000)
1373.01	721.64	3435.53	4465	3368	7781.03


abline(v=1373/1000,col="red",lwd=2)
abline(v=721/1000,col="red",lty=2)
abline(v=3435/1000,col="red",lty=2)

abline(v=4465/1000,col="dark green",lwd=2)
abline(v=3368/1000,col="dark green",lty=2)
abline(v=7781/1000,col="dark green",lty=2)
mtext("Spawners (000s)",1,outer=T,line=1)
mtext("Recruits (000s)",2,outer=T,line=1)

dev.off()


#N Coastal Streams

NCC.cu <- subset(sx.rs, CU=="N_Coastal_Streams")

alpha <- co.post[,6]; beta <- co.post[,20]

max(NCC.cu$Esc, na.rm=T)

spw <- seq(0,max(NCC.cu$Esc, na.rm=T),100)
spw
recruits <- array(NA,dim=c(1,length(spw),10000))

for(j in 1:10000){
  iter <- sample(length(alpha),1)
  recruits[,,j] <- spw*exp(alpha[iter]-beta[iter]* spw )
  
}


rec.m <- apply(recruits,c(1,2),median)/1000
rec.u <- apply(recruits,c(1,2),quantile,probs=0.975)/1000
rec.l <- apply(recruits,c(1,2),quantile,probs=0.025)/1000

pdf("COCoastalStreamsSR.pdf", width=6, height=4)

par(mfrow=c(1,1),mar=c(1.5,2.5,1,1.5),oma=c(2,2,2,1),new=F)
plot(NCC.cu$Esc/1000,NCC.cu$Rec/1000,bty='l',xlab="",ylab="",xaxt="n",yaxt="n")
axis(1)
axis(2,las=2)
mtext("Coho - Northern Coastal Streams",3,cex=0.85,line=1)

points(spw/1000, rec.m,type="l" )

polygon(c(spw/1000,rev(spw/1000)),c(rec.u,rev(rec.l)),col=grey(0.9),border=NA)
points(spw/1000, rec.m,type="l",lwd=2 )
points(NCC.cu$Esc/1000,NCC.cu$Rec/1000)
43411	35523.95	61143.52	86838	71089.95	122587.45



abline(v=43411/1000,col="red",lwd=2)
abline(v=35523/1000,col="red",lty=2)
abline(v=61143/1000,col="red",lty=2)

abline(v=86838/1000,col="dark green",lwd=2)
abline(v=71089/1000,col="dark green",lty=2)
abline(v=122587/1000,col="dark green",lty=2)
mtext("Spawners (000s)",1,outer=T,line=1)
mtext("Recruits (000s)",2,outer=T,line=1)

dev.off()


#Douglas Channel-Kitimat Arm

Doug.cu <- subset(sx.rs, CU=="Douglas_Channel-Kitimat_Arm")

alpha <- co.post[,5]; beta <- co.post[,19]

max(Doug.cu$Esc, na.rm=T)

spw <- seq(0,max(Doug.cu$Esc, na.rm=T),100)
spw
recruits <- array(NA,dim=c(1,length(spw),10000))

for(j in 1:10000){
  iter <- sample(length(alpha),1)
  recruits[,,j] <- spw*exp(alpha[iter]-beta[iter]* spw )
  
}


rec.m <- apply(recruits,c(1,2),median)/1000
rec.u <- apply(recruits,c(1,2),quantile,probs=0.975)/1000
rec.l <- apply(recruits,c(1,2),quantile,probs=0.025)/1000

jpeg("CODougKitimatSR5.jpg", width=6, height=4,units="in",res=200)
#pdf("CODouglasKitimatSR.pdf", width=6, height=4)

par(mfrow=c(1,1),mar=c(3.5,3.5,2,2),oma=c(2,2,2,1),new=F)
plot(Doug.cu$Esc/1000,Doug.cu$Rec/1000,bty='l',xlab="",ylab="",xaxt="n",yaxt="n")
axis(1)
axis(2,las=2)
mtext("Coho - Douglas Channel-Kitimat Arm",3,cex=0.85,line=1)

points(spw/1000, rec.m,type="l" )

polygon(c(spw/1000,rev(spw/1000)),c(rec.u,rev(rec.l)),col=grey(0.9),border=NA)
points(spw/1000, rec.m,type="l",lwd=2 )
points(Doug.cu$Esc/1000,Doug.cu$Rec/1000)
22123	4749.39	30294.01	45371	37536.97	60588.03


#abline(a=0,b=1,lwd=2,lty=2)
abline(v=22123/1000,col="red",lwd=2)
abline(v=4749/1000,col="red",lty=2)
abline(v=30294/1000,col="red",lty=2)

abline(v=45371/1000,col="dark green",lwd=2)
abline(v=37536/1000,col="dark green",lty=2)
abline(v=60588/1000,col="dark green",lty=2)
mtext("Spawners (000s)",1,outer=T,line=1)
mtext("Recruits (000s)",2,outer=T,line=1)

dev.off()



#############
###########
###Pink
#######

#Hecate Lowlands

Hec.cu <- subset(sx.rs, CU=="Hecate_Lowlands_even")

alpha <- pk.post[,9]; beta <- pk.post[,20]

max(Hec.cu$Esc, na.rm=T)

spw <- seq(0,max(Hec.cu$Esc, na.rm=T),10000)
spw
recruits <- array(NA,dim=c(1,length(spw),10000))

for(j in 1:10000){
  iter <- sample(length(alpha),1)
  recruits[,,j] <- spw*exp(alpha[iter]-beta[iter]* spw )
  
}


rec.m <- apply(recruits,c(1,2),median)/1000
rec.u <- apply(recruits,c(1,2),quantile,probs=0.975)/1000
rec.l <- apply(recruits,c(1,2),quantile,probs=0.025)/1000

pdf("PKHecateLowEvenSR.pdf", width=6, height=4)

par(mfrow=c(1,1),mar=c(1.5,2.5,1,1.5),oma=c(2,2,2,1),new=F)
plot(Hec.cu$Esc/1000,Hec.cu$Rec/1000,bty='l',xlab="",ylab="",xaxt="n",yaxt="n")
axis(1)
axis(2,las=2)
mtext("Pink - Hecate Lowlands Even",3,cex=0.85,line=1)

points(spw/1000, rec.m,type="l" )

polygon(c(spw/1000,rev(spw/1000)),c(rec.u,rev(rec.l)),col=grey(0.9),border=NA)
points(spw/1000, rec.m,type="l",lwd=2 )
points(Hec.cu$Esc/1000,Hec.cu$Rec/1000)
409621	0	2.71801E+13	1449367	508087	5.43602E+13


abline(v=409621/1000,col="red",lwd=2)
abline(v=0/1000,col="red",lty=2)
abline(v=2718010000000/1000,col="red",lty=2)

abline(v=1449367/1000,col="dark green",lwd=2)
abline(v=508087/1000,col="dark green",lty=2)
abline(v=5436020000000/1000,col="dark green",lty=2)
mtext("Spawners (000s)",1,outer=T,line=1)
mtext("Recruits (000s)",2,outer=T,line=1)

dev.off()


#Hecate Strait-Fjords

Hec.cu <- subset(sx.rs, CU=="Hecate_Strait-Fjords_even")

alpha <- pk.post[,10]; beta <- pk.post[,21]

max(Hec.cu$Esc, na.rm=T)

spw <- seq(0,max(Hec.cu$Esc, na.rm=T),10000)
spw
recruits <- array(NA,dim=c(1,length(spw),10000))

for(j in 1:10000){
  iter <- sample(length(alpha),1)
  recruits[,,j] <- spw*exp(alpha[iter]-beta[iter]* spw )
  
}


rec.m <- apply(recruits,c(1,2),median)/1000
rec.u <- apply(recruits,c(1,2),quantile,probs=0.975)/1000
rec.l <- apply(recruits,c(1,2),quantile,probs=0.025)/1000

pdf("PKHecateFjordsEvenSR.pdf", width=6, height=4)

par(mfrow=c(1,1),mar=c(1.5,2.5,1,1.5),oma=c(2,2,2,1),new=F)
plot(Hec.cu$Esc/1000,Hec.cu$Rec/1000,bty='l',xlab="",ylab="",xaxt="n",yaxt="n")
axis(1)
axis(2,las=2)
mtext("Pink - Hecate Strait-Fjords Even",3,cex=0.85,line=1)

points(spw/1000, rec.m,type="l" )

polygon(c(spw/1000,rev(spw/1000)),c(rec.u,rev(rec.l)),col=grey(0.9),border=NA)
points(spw/1000, rec.m,type="l",lwd=2 )
points(Hec.cu$Esc/1000,Hec.cu$Rec/1000)
1710253	0	7.60759E+13	5595447	2178434.52	1.52152E+14


abline(v=1710253/1000,col="red",lwd=2)
abline(v=0/1000,col="red",lty=2)
abline(v=2718010000000/1000,col="red",lty=2)

abline(v=5595447/1000,col="dark green",lwd=2)
abline(v=2178434/1000,col="dark green",lty=2)
abline(v=5436020000000/1000,col="dark green",lty=2)
mtext("Spawners (000s)",1,outer=T,line=1)
mtext("Recruits (000s)",2,outer=T,line=1)

dev.off()



#Homathko-Klinaklini-Rivers-Smith-Bella Coola Dean

Hec.cu <- subset(sx.rs, CU=="Homathko-Klinaklini-Rivers-Smith-Bella_Coola_Dean_odd")

alpha <- pk.post[,11]; beta <- pk.post[,22]

max(Hec.cu$Esc, na.rm=T)

spw <- seq(0,max(Hec.cu$Esc, na.rm=T),10000)
spw
recruits <- array(NA,dim=c(1,length(spw),10000))

for(j in 1:10000){
  iter <- sample(length(alpha),1)
  recruits[,,j] <- spw*exp(alpha[iter]-beta[iter]* spw )
  
}


rec.m <- apply(recruits,c(1,2),median)/1000
rec.u <- apply(recruits,c(1,2),quantile,probs=0.975)/1000
rec.l <- apply(recruits,c(1,2),quantile,probs=0.025)/1000

pdf("PKHomathkoOddSR.pdf", width=6, height=4)

par(mfrow=c(1,1),mar=c(1.5,2.5,1,1.5),oma=c(2,2,2,1),new=F)
plot(Hec.cu$Esc/1000,Hec.cu$Rec/1000,bty='l',xlab="",ylab="",xaxt="n",yaxt="n")
axis(1)
axis(2,las=2)
mtext("Pink - Homathko-Klinaklini-Rivers-Smith-Bella Coola Dean Odd",3,cex=0.85,line=1)

points(spw/1000, rec.m,type="l" )

polygon(c(spw/1000,rev(spw/1000)),c(rec.u,rev(rec.l)),col=grey(0.9),border=NA)
points(spw/1000, rec.m,type="l",lwd=2 )
points(Hec.cu$Esc/1000,Hec.cu$Rec/1000)
492353.75	334789.42	1130347.93	989185	693272.78	3075538.73


abline(v=492353/1000,col="red",lwd=2)
abline(v=334789/1000,col="red",lty=2)
abline(v=1130347/1000,col="red",lty=2)

abline(v=989185/1000,col="dark green",lwd=2)
abline(v=693272/1000,col="dark green",lty=2)
abline(v=3075538/1000,col="dark green",lty=2)
mtext("Spawners (000s)",1,outer=T,line=1)
mtext("Recruits (000s)",2,outer=T,line=1)

dev.off()



#Hecate Strait-Lowlands

Hec.cu <- subset(sx.rs, CU=="Hecate_Strait-Lowlands_odd")

alpha <- pk.post[,2]; beta <- pk.post[,13]

max(Hec.cu$Esc, na.rm=T)

spw <- seq(0,max(Hec.cu$Esc, na.rm=T),10000)
spw
recruits <- array(NA,dim=c(1,length(spw),10000))

for(j in 1:10000){
  iter <- sample(length(alpha),1)
  recruits[,,j] <- spw*exp(alpha[iter]-beta[iter]* spw )
  
}


rec.m <- apply(recruits,c(1,2),median)/1000
rec.u <- apply(recruits,c(1,2),quantile,probs=0.975)/1000
rec.l <- apply(recruits,c(1,2),quantile,probs=0.025)/1000

pdf("PKHecLowOddSR.pdf", width=6, height=4)

par(mfrow=c(1,1),mar=c(1.5,2.5,1,1.5),oma=c(2,2,2,1),new=F)
plot(Hec.cu$Esc/1000,Hec.cu$Rec/1000,bty='l',xlab="",ylab="",xaxt="n",yaxt="n")
axis(1)
axis(2,las=2)
mtext("Pink - Hecate Strait-Lowlands Odd",3,cex=0.85,line=1)

points(spw/1000, rec.m,type="l" )

polygon(c(spw/1000,rev(spw/1000)),c(rec.u,rev(rec.l)),col=grey(0.9),border=NA)
points(spw/1000, rec.m,type="l",lwd=2 )
points(Hec.cu$Esc/1000,Hec.cu$Rec/1000)
246283	185075.71	525848.25	494166	377094.83	1310036.05


abline(v=246283/1000,col="red",lwd=2)
abline(v=185075/1000,col="red",lty=2)
abline(v=525848/1000,col="red",lty=2)

abline(v=494166/1000,col="dark green",lwd=2)
abline(v=377094/1000,col="dark green",lty=2)
abline(v=1310036/1000,col="dark green",lty=2)
mtext("Spawners (000s)",1,outer=T,line=1)
mtext("Recruits (000s)",2,outer=T,line=1)

dev.off()

#Hecate Strait-Fjords


Hec.cu <- subset(sx.rs, CU=="Hecate_Strait-Fjords_odd")

alpha <- pk.post[,3]; beta <- pk.post[,14]

max(Hec.cu$Esc, na.rm=T)

spw <- seq(0,max(Hec.cu$Esc, na.rm=T),10000)
spw
recruits <- array(NA,dim=c(1,length(spw),10000))

for(j in 1:10000){
  iter <- sample(length(alpha),1)
  recruits[,,j] <- spw*exp(alpha[iter]-beta[iter]* spw )
  
}


rec.m <- apply(recruits,c(1,2),median)/1000
rec.u <- apply(recruits,c(1,2),quantile,probs=0.975)/1000
rec.l <- apply(recruits,c(1,2),quantile,probs=0.025)/1000

pdf("PKHecFjordsOddSR.pdf", width=6, height=4)

par(mfrow=c(1,1),mar=c(1.5,2.5,1,1.5),oma=c(2,2,2,1),new=F)
plot(Hec.cu$Esc/1000,Hec.cu$Rec/1000,bty='l',xlab="",ylab="",xaxt="n",yaxt="n")
axis(1)
axis(2,las=2)
mtext("Pink - Hecate Strait-Fjords Odd",3,cex=0.85,line=1)

points(spw/1000, rec.m,type="l" )

polygon(c(spw/1000,rev(spw/1000)),c(rec.u,rev(rec.l)),col=grey(0.9),border=NA)
points(spw/1000, rec.m,type="l",lwd=2 )
points(Hec.cu$Esc/1000,Hec.cu$Rec/1000)
949715	0	2.41632E+13	2582148	1179970.25	4.83264E+13


abline(v=949715/1000,col="red",lwd=2)
abline(v=0/1000,col="red",lty=2)
#abline(v=525848/1000,col="red",lty=2)

abline(v=2582148/1000,col="dark green",lwd=2)
abline(v=1179970/1000,col="dark green",lty=2)
#abline(v=1310036/1000,col="dark green",lty=2)
mtext("Spawners (000s)",1,outer=T,line=1)
mtext("Recruits (000s)",2,outer=T,line=1)

dev.off()


# mussel

Mussel.cu <- subset(sx.rs, CU=="Mussel-Kynoch")

max(Mussel.cu$Esc, na.rm=T)

spw <- seq(0,max(Mussel.cu$Esc, na.rm=T),100)
spw
recruits <- array(NA,dim=c(1,length(spw),10000))

for(j in 1:10000){
  iter <- sample(length(alpha),1)
  recruits[,,j] <- spw*exp(alpha[iter]-beta[iter]* spw )
  
}


jpeg("COMusselKynockPercentile.jpg", width=6, height=4, units="in",res=200)

par(mfrow=c(1,1),mar=c(1.5,2.5,1,1.5),oma=c(2,2,2,1),new=F)
plot(Mussel.cu$BY,Mussel.cu$Esc/1000,bty='l', col="dark blue", xlab="",ylab="",xaxt="n",yaxt="n", pch=16)
axis(1)
axis(2,las=2)
mtext("Coho - Mussel-Kynoch",3,cex=0.85,line=1)

points(Mussel.cu$BY,Mussel.cu$Esc/1000, col="dark blue", type="l")

abline(h=2880/1000,col="red",lwd=2, lty=2)
abline(h=8600/1000,col="green",lty=2, lwd=2)

dev.off()
