#!/usr/bin/env Rscript

#################################################
### script to plot posterior distribution     ###
###    under habitat loss hypothesis          ###
###          (output from MIDASPOM)           ###
###                 N. Alcala                 ###
#################################################

print("Running Rscript to test hypotheses and plot the posterior density of the parameters")

#recover arguments of script
library("optparse")

option_list = list(
  make_option(c("-o", "--out"), type="character", default="out", help="output directory name [default= %default]", metavar="character"),
  make_option(c("-p", "--pattern"), type="character", default="dieoff", help="pattern for input file names [default= %default]", metavar="character"),
  make_option(c("-k", "--Kmin"), type="numeric", default=0.1, help="minimum population size [default= %default]", metavar="number"),
  make_option(c("-K", "--Kmax"), type="numeric", default=100, help="maximum population size [default= %default]", metavar="number"),
  make_option(c("-d", "--dmin"), type="numeric", default=0.1, help="minimum source population distance [default= %default]", metavar="number"),
  make_option(c("-D", "--dmax"), type="numeric", default=100, help="maximum source population distance [default= %default]", metavar="number"),
  make_option(c("-t", "--tmin"), type="numeric", default=1, help="minimum year [default= %default]", metavar="number"),
  make_option(c("-T", "--tmax"), type="numeric", default=1, help="maximum year [default= %default]", metavar="number"),
  make_option(c("-n", "--segments"), type="numeric", default=10, help="number of segments [default= %default]", metavar="number"),
  make_option(c("-l", "--lH0"), type="numeric", default=1, help="likelihood of null hypothesis [default= %default]", metavar="number")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#parameters
postfilespattern = opt$pattern#args[1] #file with posterior distribution of parameters under in habitat loss hypothesis
outfile          = opt$out # args[2] #file with posterior distribution
Kmin  = as.numeric( opt$Kmin) #args[3])
Kmax  = as.numeric( opt$Kmax) #args[4])
dmin  = as.numeric( opt$dmin) #args[5])
dmax  = as.numeric( opt$dmax) #args[6])
tmin  = as.numeric( opt$tmin) #args[7])
tmax  = as.numeric( opt$tmax) #args[8])
lH0   = as.numeric( opt$lH0) #args[9])
n = as.numeric( opt$segments) #10
print("Arguments read")

#load joint posterior density
postfiles  = grep(postfilespattern, list.files(),value = T,perl = T)
jpostlossl = lapply(postfiles, function(x) as.matrix(read.table(x)) ) #matrix(scan(postfiledieoff),1,stepK,byrow=T)
print("Read posterior densities")
#print(jpostlossl)
stepK = nrow(jpostlossl[[1]])
stepd = ncol(jpostlossl[[1]])
stept =length(jpostlossl)
print(stept)

Kl    = 10**seq(log(Kmin,10),log(Kmax,10),length.out=stepK)
dl    = seq(dmin,dmax,length.out=stepd)
tl = seq(tmin,tmax,length.out=stept)

jpostloss = matrix(0,stepK,stepd)
for(i in 1:length(jpostlossl)) jpostloss = jpostloss + jpostlossl[[i]]
jpostloss = jpostloss/sum(jpostloss)/(log(Kmax,10)-log(Kmin,10))*(stepK-1)/(dmax-dmin)*(stepd-1)
postT = sapply(jpostlossl, sum)
postT = postT/sum(postT)/(tmax-tmin)*(stept-1)

#point estimates
maxpostdloss  = max(colSums(jpostloss))
maxpostKloss  = max(rowSums(jpostloss))
dlest         = dl[colSums(jpostloss)>=(0.95*maxpostdloss)][1]
Klest         = Kl[rowSums(jpostloss)>=(0.95*maxpostKloss)][1]
tdest         = tl[postT>=(0.95*max(postT))][1] #value for which the posterior distribution is 95% of its maximum value

#plot
mat = max(postT)

coco = rgb(1,seq(1,0,-0.1),seq(1,0,-0.1),1)
svg(outfile,h=1*3.5,w=2*3.5)
par(mfrow=c(1,2),family="Times",las=1)
image(log(Kl,10),dl,jpostloss,col=coco,xlab=expression(italic(K[L])),ylab=expression(italic(d[L])) ,main="Hypothesis 2: source population loss",useRaster=T,axes=F)
axis(1,at = log(c(0.1,0.5,1,5,10,50,100),10), labels =  c(0.1,0.5,1,5,10,50,100) )
axis(2,at=seq(200,4000,400))
abline(v=log(Klest,10),col=2,lty=2,lwd=2)
abline(h=dlest,col=2,lty=2,lwd=2)
mtext(format(Klest,digits=3),side=3,at=log(Klest,10),col=2,line=0.5)
mtext(format(dlest,digits=3),side=4,at=dlest,col=2,line=0.5)

plot(0.0001,-1,type="l",xlab=expression(italic(t[L])),ylab="Density",col=2,ylim=c(0,mat),lwd=2,xlim=c(tmin,tmax),main="Hypothesis 2: source population loss")
polygon(c(tmin,tmin,tmax,tmax),c(0,1/(tmax-tmin),1/(tmax-tmin),0),lwd=2,col=rgb(0,0,0,0.4),border=NA)
polygon(c(tl,rev(tl)), c(postT,rep(0,stept) ),col=rgb(1,0,0,0.5),border=NA ) 
abline(v = tdest,col=2,lty=2,lwd=2)
mtext(format(tdest,digits=3),side=3,at=tdest,col=2,line=0.5)
dev.off()

dlest = dl[(which.max(jpostloss)-1)%/%151 +1]
Klest = Kl[(which.max(jpostloss)-1)%%151 +1]

postd = colSums(jpostloss)
postK = rowSums(jpostloss)
print(postd)
#print(postK)
#print(postT)
print(dl)
print( paste("Kl est=, CI=[", Kl[sort(which(cumsum(postK/sum(postK))<=0.05),decreasing=T)][1],",", Kl[sort(which(cumsum(postK/sum(postK))<=1),decreasing=T)][1],"]") )
print( paste("dl est=, CI=[", dl[sort(which(cumsum(postd/sum(postd))<=0),decreasing=T)][1],",", dl[sort(which(cumsum(postd/sum(postd))<=0.95),decreasing=T)][1],"]") )
print( paste("tl CI: [", tl[sort(which(cumsum(postT/sum(postT))<=0.05),decreasing=T)][1],",", tl[sort(which(cumsum(postT/sum(postT))<=1),decreasing=T)][1],"]") )


print(paste("Bayes Factor: ",log( lH0/mean(unlist(jpostlossl)) ,10) ,sep="") )
print(paste("AIC0: ",-2*log(lH0/2**n) ,sep="") )
print(paste("AIC2: ",3*2-2*log(max(unlist(jpostlossl))/2**n) ,sep="") )

print(paste("plotted to file",outfile))
