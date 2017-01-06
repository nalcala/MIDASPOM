#################################################
### script to plot posterior distribution     ###
### under in situ die-off hypothesis          ###
###          (output from MIDASPOM)           ###
###                 N. Alcala                 ###
#################################################

print("Running Rscript to plot the posterior density of the parameters under in situ die-off hypothesis")

#recover arguments of script
args <- commandArgs(trailingOnly = TRUE)

#parameters
postfiledieoff = args[1] #file with posterior distribution of parameters under in situ die-off hypothesis
outfile        = args[2] #file with posterior distribution
Kmin  = as.numeric(args[3])
Kmax  = as.numeric(args[4])
stepK = as.numeric(args[5])
Kl = 10**seq(log(Kmin,10),log(Kmax,10),length.out=stepK)

#load likelihood
postdieoff = scan(postfiledieoff) #matrix(scan(postfiledieoff),1,stepK,byrow=T)
#compute posterior
postK = postdieoff/sum(postdieoff)/(log(Kmax,10)-log(Kmin,10))*(stepK-1)

#point estimates
Kdest         = Kl[postK>=(0.95*max(postK))][1] #value for which the posterior distribution is 95% of its maximum value

ma = max(postK)
#plot
coco = rgb(1,seq(1,0,-0.1),seq(1,0,-0.1),1)
pdf(outfile,h=1*3.5,w=1*3.5)
par(mfrow=c(1,1),family="Times",las=1)
plot(0.0001,-1,type="l",log="x",xlab=expression(italic(K[D])),ylab="density",col=2,ylim=c(0,ma),lwd=2,xlim=c(Kmin,Kmax),main="Hypothesis 1: in situ die-off")
polygon(c(Kmin,Kmin,Kmax,Kmax),c(0,1,1,0)/(log(Kmax,10)-log(Kmin,10)),lwd=2,col=rgb(0,0,0,0.4),border=NA)
polygon(c(Kl,rev(Kl)), c(postK,rep(0,stepK) ),col=rgb(1,0,0,0.5),border=NA ) 
abline(v = Kdest,col=2,lty=2,lwd=2)
mtext(format(Kdest,digits=3),side=1,at=Kdest,col=2,line=0.5)
dev.off()

print(paste("plotted to file",outfile))
