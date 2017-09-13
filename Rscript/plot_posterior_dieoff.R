#!/usr/bin/env Rscript

#################################################
### script to plot posterior distribution     ###
### under in situ die-off hypothesis          ###
###          (output from MIDASPOM)           ###
###                 N. Alcala                 ###
#################################################

print("Running Rscript to plot the posterior density of the parameters under in situ die-off hypothesis")

#recover arguments of script
library("optparse")

option_list = list(
  make_option(c("-o", "--out"), type="character", default="out", help="output directory name [default= %default]", metavar="character"),
  make_option(c("-p", "--pattern"), type="character", default="dieoff", help="pattern for input file names [default= %default]", metavar="character"),
  make_option(c("-k", "--Kmin"), type="numeric", default=0.1, help="minimum population size [default= %default]", metavar="number"),
  make_option(c("-K", "--Kmax"), type="numeric", default=100, help="maximum population size [default= %default]", metavar="number"),
  make_option(c("-t", "--tmin"), type="numeric", default=1, help="minimum year [default= %default]", metavar="number"),
  make_option(c("-T", "--tmax"), type="numeric", default=1, help="maximum year [default= %default]", metavar="number"),
  make_option(c("-n", "--segments"), type="numeric", default=10, help="number of segments [default= %default]", metavar="number")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#parameters
postfilespattern = opt$pattern #file with posterior distribution of parameters under in situ die-off hypothesis
outfile        = opt$out #file with posterior distribution
Kmin  = as.numeric(opt$Kmin)
Kmax  = as.numeric(opt$Kmax)
tmin  = as.numeric(opt$tmin)
tmax  = as.numeric(opt$tmax)
n     = as.numeric(opt$segments)

#load likelihood
postfiles  = grep(postfilespattern, list.files(),value = T,perl = T)
postdieoff = sapply(postfiles, scan) #matrix(scan(postfiledieoff),1,stepK,byrow=T)
if(is.array( postdieoff )){
  stepK = nrow(postdieoff)
  stept = ncol(postdieoff)
}else{
  stepK = length(postdieoff)
  stept=1
}
print(stepK)
Kl = 10**seq(log(Kmin,10),log(Kmax,10),length.out=stepK)
tl = seq(tmin,tmax,length.out=stept)
print(tl)

#compute posterior
jpost = postdieoff/sum(postdieoff)/(log(Kmax,10)-log(Kmin,10))*(stepK-1)
postK = rowSums(jpost)
postT = colSums(jpost)/sum(jpost)/(tmax-tmin)*(stept-1)

#point estimates
Kdest         = Kl[postK>=(0.95*max(postK))][1] #value for which the posterior distribution is 95% of its maximum value
tdest         = tl[postT>=(0.95*max(postT))][1] #value for which the posterior distribution is 95% of its maximum value

ma  = max(postK)
mat = max(postT)
#plot
coco = rgb(1,seq(1,0,-0.1),seq(1,0,-0.1),1)
svg(outfile,h=1*3.5,w=2*3.5)
par(mfrow=c(1,2),family="Times",las=1)
plot(0.0001,-1,type="l",log="x",xlab=expression(italic(K[D])),ylab="Density",col=2,ylim=c(0,ma),lwd=2,xlim=c(Kmin,Kmax),main="Hypothesis 1: in situ die-off")
polygon(c(Kmin,Kmin,Kmax,Kmax),c(0,1,1,0)/(log(Kmax,10)-log(Kmin,10)),lwd=2,col=rgb(0,0,0,0.4),border=NA)
polygon(c(Kl,rev(Kl)), c(postK,rep(0,stepK) ),col=rgb(1,0,0,0.5),border=NA ) 
abline(v = Kdest,col=2,lty=2,lwd=2)
abline(v = 1,col=1,lty=2,lwd=2)
mtext(format(Kdest,digits=3),side=3,at=Kdest,col=2,line=0.5)

plot(0.0001,-1,type="l",xlab=expression(italic(t[D])),ylab="Density",col=2,ylim=c(0,mat),lwd=2,xlim=c(tmin,tmax),main="Hypothesis 1: in situ die-off")
polygon(c(tmin,tmin,tmax,tmax),c(0,1/(tmax-tmin),1/(tmax-tmin),0),lwd=2,col=rgb(0,0,0,0.4),border=NA)
polygon(c(tl,rev(tl)), c(postT,rep(0,stept) ),col=rgb(1,0,0,0.5),border=NA ) 
abline(v = tdest,col=2,lty=2,lwd=2)
mtext(format(tdest,digits=3),side=3,at=tdest,col=2,line=0.5)
dev.off()

print(paste("plotted to file",outfile))
print(mean(postdieoff[51,]))
print(mean(postdieoff))

 
tl[which(postT==max(postT))]

print( paste("Kd est=, CI=[", Kl[sort(which(cumsum(postK/sum(postK))<=0.05),decreasing=T)][1],",", Kl[sort(which(cumsum(postK/sum(postK))<=1),decreasing=T)][1],"]") )
print( paste("td CI: [", tl[sort(which(cumsum(postT/sum(postT))<=0.05),decreasing=T)][1],",", tl[sort(which(cumsum(postT/sum(postT))<=1),decreasing=T)][1],"]") )

print(paste("Bayes Factor: ",log( mean(postdieoff[51,])/mean(postdieoff) ,10) ,sep="") )
print(paste("AIC0: ",-2*log(max(postdieoff[51,])/2**n) ,sep="") )
print(paste("AIC1: ",2*2-2*log(max(postdieoff)/2**n) ,sep="") )
