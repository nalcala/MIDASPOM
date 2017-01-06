#################################################
### script to perform hypothesis testing      ###
###          (output from MIDASPOM)           ###
###                 N. Alcala                 ###
#################################################

print("Running Rscript to test hypotheses and plot AIC and Bayes factors")

#recover arguments of script
args <- commandArgs(trailingOnly = TRUE)

#parameters
postfiledieoff = args[1] #file with posterior distribution of parameters under in situ die-off hypothesis
postfileloss   = args[2] #file with posterior distribution of parameters under in habitat loss hypothesis
outfile        = args[3] #file with posterior distribution
n      = as.numeric(args[4])
Kmin   = as.numeric(args[5])
Kmax   = as.numeric(args[6])

#load joint posterior density
postdieoff = scan(postfiledieoff)
jpostloss  = as.matrix(read.table(postfileloss))
stepKl     = nrow(jpostloss)
stepd      = ncol(jpostloss)

stepKd = length(postdieoff)
Kdl    = 10**seq(log(Kmin,10),log(Kmax,10),length.out=stepKd)

# find likelihood of null model or approximate it
if( sum(Kdl==1)>0 ){ # null model is present
  idnull = which(Kdl==1) 
  postnull = postdieoff[idnull] # retrieve likelihood of null
}else{  # null model is not present
  idnull   = which(Kdl >1)[1] 
  postnull = (postdieoff[idnull]-postdieoff[idnull-1])/(Kdl[idnull]-Kdl[idnull-1])*(1-Kdl[idnull-1]) + postdieoff[idnull-1] #approximate likelihood using linear regression
}

# compute AIC
AIC0 = 2*1-2*log(postnull/2**n)
AIC1 = 2*2-2*log(max(postdieoff)/2**n)
AIC2 = 2*3-2*log( max(jpostloss)/2**n) 

# compute Bayes factor
lK01 = log( postnull/(sum(postdieoff)/stepKl) ,10)
lK02 = log( postnull/(sum(jpostloss)/stepKl/stepd) ,10)
lK12 = log( sum(postdieoff)/(sum(jpostloss)/stepd) ,10)

#plot
coco = rgb(1,seq(1,0,-0.1),seq(1,0,-0.1),1)
pdf(outfile,h=1*3.5,w=2.5*3.5)
par(mfrow=c(1,2),family="Times",las=1)
barplot(c(AIC0,AIC1,AIC2),main="AIC",names.arg=c("H0","H1","H2"))
barplot(c(lK01,lK02,lK12),main="Log Bayes factor",names.arg=c("H0 vs H1","H0 vs H2","H1 vs H2"))
dev.off()

print(paste("plotted to file",outfile))
