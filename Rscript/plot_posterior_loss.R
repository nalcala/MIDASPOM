#################################################
### script to plot posterior distribution     ###
###    under habitat loss hypothesis          ###
###          (output from MIDASPOM)           ###
###                 N. Alcala                 ###
#################################################

print("Running Rscript to test hypotheses and plot the posterior density of the parameters")

#recover arguments of script
args <- commandArgs(trailingOnly = TRUE)

#parameters
postfileloss   = args[1] #file with posterior distribution of parameters under in habitat loss hypothesis
outfile        = args[2] #file with posterior distribution
Kmin  = as.numeric(args[3])
Kmax  = as.numeric(args[4])
dmin  = as.numeric(args[5])
dmax  = as.numeric(args[6])
stepK = as.numeric(args[7])
stepd = as.numeric(args[8])
Kl    = 10**seq(log(Kmin,10),log(Kmax,10),length.out=stepK)
dl    = seq(dmin,dmax,length.out=stepd)

#load joint posterior density
jpostloss   = matrix(scan(postfileloss),stepK,stepd,byrow=T)

#point estimates
maxpostKloss  = max(colSums(jpostloss))
maxpostdloss  = max(rowSums(jpostloss))
dlest         = dl[colSums(jpostloss)>=(0.95*maxpostKloss)][1]
Klest         = Kl[rowSums(jpostloss)>=(0.95*maxpostdloss)][1]

#plot
coco = rgb(1,seq(1,0,-0.1),seq(1,0,-0.1),1)
pdf(outfile,h=1*3.5,w=1*3.5)
par(mfrow=c(1,1),family="Times",las=1)
image(log(Kl,10),dl,jpostloss,col=coco,xlab=expression(log[10]~italic(K[L])),ylab=expression(italic(d[L])) ,main="Hypothesis 2: habitat loss",useRaster=T)
abline(v=log(Klest,10),col=2,lty=2,lwd=2)
abline(h=dlest,col=2,lty=2,lwd=2)
mtext(format(log(Klest,10),digits=3),side=1,at=log(Klest,10),col=2,line=0.5)
mtext(format(dlest,digits=3),side=2,at=dlest,col=2,line=0.5)

dev.off()

print(paste("plotted to file",outfile))
