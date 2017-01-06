#################################################
### script to plot the posterior distribution ###
###          (output from MIDASPOM)           ###
###                 N. Alcala                 ###
#################################################

print("Running Rscript to plot the posterior density")

#recover arguments of script
args <- commandArgs(trailingOnly = TRUE)

#parameters
postfile = args[1] #file with posterior distribution
outfile  = args[2] #file with posterior distribution
l = as.numeric(args[3])
u = as.numeric(args[4])
step = as.numeric(args[5])
el = seq(l,u,length.out=step)
cl = seq(l,u,length.out=step)

#load joint posterior density
jpost = as.matrix(read.table(postfile)) #matrix(scan(postfile),step,step,byrow=T)
jpost = jpost/sum(jpost)/0.01/0.01

#point estimates
ML   = which.max(jpost)
eest = el[(ML-1)%%step +1]
cest = cl[(ML-1)%/%step +1]

#compute marginal posterior densities
epost = rowMeans(jpost)
cpost = colMeans(jpost)

qel = which((cumsum(epost)>=(0.025/0.01) )&(cumsum(epost)<(0.975/0.01)))
qcl = which((cumsum(cpost)>=(0.025/0.01))&(cumsum(cpost)<(0.975/0.01)))

#plot
pmax = max(c(epost,cpost))
coco = rgb(1,seq(1,0,-0.1),seq(1,0,-0.1),1)
pdf(outfile,h=1*2.5,w=3*2.5)
svg("FigS1DC.svg",h=1*2.5,w=3*2.5)
par(mfrow=c(1,3),family="Times",las=1)
plot( -1, -1,type="l",xlim=c(l,u),ylim=c(0,8), xlab="extinction parameter",ylab="density",col=2,lwd=3,las=1,main=" ")
polygon(c(l,l,u,u),c(0,1/(u-l),1/(u-l),0),lwd=2,col=rgb(0,0,0,0.3),border=NA)
polygon(c(el,rev(el)), c(epost,rep(0,step) ),col=rgb(1,0,0,0.3),border=NA ) 
polygon(c(el[qel],rev(el[qel])), c(epost[qel],rep(0,step)[qel] ),col=rgb(1,0,0,0.5),border=NA ) 
abline(v=el[qel[1]],col=2,lwd=3,lty=2)
abline(v=el[rev(qel)[1]],col=2,lwd=3,lty=2)
#mtext(eest,side=1,at=eest,col=2,line=0.5)

plot( -1, -1,type="l",xlim=c(l,u),ylim=c(0,8), xlab="colonization parameter",ylab="density",col=2,lwd=3,las=1)
polygon(c(l,l,u,u),c(0,1/(u-l),1/(u-l),0),lwd=2,col=rgb(0,0,0,0.3),border=NA)
polygon(c(cl,rev(cl)), c(cpost,rep(0,step) ),col=rgb(1,0,0,0.5),border=NA ) 
polygon(c(cl[qcl],rev(cl[qcl])), c(cpost[qcl],rep(0,step)[qcl] ),col=rgb(1,0,0,0.5),border=NA ) 
abline(v=cl[qcl[1]],col=2,lwd=3,lty=2)
abline(v=cl[rev(qcl)[1]],col=2,lwd=3,lty=2)

image(el,cl,jpost,col=coco,asp=1,xlab="extinction parameter",ylab="colonization parameter",useRaster = T)
abline(v=eest,lty=2,col=2)
abline(h=cest,lty=2,col=2)
mtext(eest,side=1,at=eest,col=2,line=0.5)
mtext(cest,side=2,at=cest,col=2,line=0.5)
dev.off()

print(paste("plotted to file",outfile))
