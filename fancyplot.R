#plot bullshit
par(mfrow=c(1,1))



n=1:5;
space=.1

labs=c('n=25','n=20','n=15','n=10','n=5');
#n=5
plot(HPDinterval(as.mcmc(P1.1.5)),rep(n[5]+space,2),xlab="",ylab="",xaxt='n',yaxt='n',type='l',ylim=c(.5,5.5),xlim=c(0,1),frame.plot=F,
     lwd=2)
lines(HPDinterval(as.mcmc(P2.1.5)),rep(n[5]-space,2),lwd=2,col=4)
points(c(mean(P1.1.5),mean(P2.1.5)),(n[5]+c(space,-space)),pch=16,cex=1,col=c(1,4))
#n=10
lines(HPDinterval(as.mcmc(P1.1.10)),rep(n[4]+space,2),lwd=2)
lines(HPDinterval(as.mcmc(P2.1.10)),rep(n[4]-space,2),lwd=2,col=4)
points(c(mean(P1.1.10),mean(P2.1.10)),(n[4]+c(space,-space)),pch=16,cex=1,col=c(1,4))
#n=15
lines(HPDinterval(as.mcmc(P1.1.15)),rep(n[3]+space,2),lwd=2)
lines(HPDinterval(as.mcmc(P2.1.15)),rep(n[3]-space,2),lwd=2,col=4)
points(c(mean(P1.1.15),mean(P2.1.15)),(n[3]+c(space,-space)),pch=16,cex=1,col=c(1,4))
#n=20
lines(HPDinterval(as.mcmc(P1.1.20)),rep(n[2]+space,2),lwd=2)
lines(HPDinterval(as.mcmc(P2.1.20)),rep(n[2]-space,2),lwd=2,col=4)
points(c(mean(P1.1.20),mean(P2.1.20)),(n[2]+c(space,-space)),pch=16,cex=1,col=c(1,4))
#n=25

lines(HPDinterval(as.mcmc(P1.1.25)),rep(n[1]+space,2),lwd=2)
lines(HPDinterval(as.mcmc(P2.1.25)),rep(n[1]-space,2),lwd=2,col=4)
points(c(mean(P1.1.25),mean(P2.1.25)),(n[1]+c(space,-space)),pch=16,cex=1,col=c(1,4))
axis(2,at=n,labels=labs,las=2,tick=F)
axis(1,at=seq(0,1,.1),labels=seq(0,1,.1),tick=F)
abline(h=.5)

par(cex=.8)
legend("bottomright",legend=c(expression(paste(gamma[1]," (constant)"),paste(gamma[2]," (functional)"))),
       col=c(1,4), lty=1, cex=1.2,inset=c(.1,.1))
