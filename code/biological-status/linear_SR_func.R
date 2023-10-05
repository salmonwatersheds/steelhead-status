
#Function which computes linear-regresion estimates of parameters and MSY parameters and plots
"LinReg"=function(Nyrs,LNRS,S, R,StNames){
	
	Nstocks=dim(S)[2]

	a=vector(length=Nstocks);b=a	#;tau=a
	ngrows=4;ngcol=trunc(Nstocks/ngrows)+1
	par(mfcol=c(ngrows,ngcol),mai=c(.3,.3,.3,.2),omi=c(0.1,0.1,0.1,0.1))
	

	for (i in 1:Nstocks) {

		reg=lm(LNRS[1:Nyrs[i],i]~S[1:Nyrs[i],i])
		a[i]=as.double(reg$coefficients[1]);b[i]=as.double(abs(reg$coefficients[2]))	#;tau[i]=as.double(1/sd(reg$residuals)^2)

		#Compute production parameters
		b1=a[i]/b[i]
		Prod=round(exp(a[i]),digits=1);Smsy=round(b1*(0.5-0.07*a[i]),digits=0)
		Smax=round(b1/a[i],digits=0);Uopt=round(0.5*a[i]-0.07*a[i]^2,digits=2)
		print(c(a[i],b[i]));print(c("Prod=",Prod));print(c("Smsy=",Smsy));print(c("Smax=",Smax));print(c("Uopt=",Uopt))

		#Plot basic SR curves and data
		#plot(S[,i], LNRS[1:Nyrs[i],i],bty='l',xlab="Spawners",xlim=c(0,max(S)),ylim=c(0,max(LNRS)),ylab="Ln(R/S)",main=StNames[i])
		#abline(reg,lty=1,lwd=2)

		Sx=seq(0,max(S[1:Nyrs[i],i]),by=100)
		pR=Sx*exp(a[i]-b[i]*Sx)
		plot(S[1:Nyrs[i],i],R[1:Nyrs[i],i],bty='l',xlab="",ylab="",xlim=c(0,max(S[1:Nyrs[i],i])),ylim=c(0,max(R[1:Nyrs[i],i])),main=StNames[i])
		lines(pR~Sx,lty=1,lwd=2)
		abline(a=0,b=1,lty=2)
	}
	return(c(a,b))
}
