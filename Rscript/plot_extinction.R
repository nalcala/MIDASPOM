##################################################
### script to plot the future extinction proba ###
###          (output from MIDASPOM)            ###
###                 N. Alcala                  ###
##################################################

print("Running Rscript to plot the future probability of extinction")

#recover arguments of script
args <- commandArgs(trailingOnly = TRUE)

#parameters
prextfile = args[1] #file with probability of extinction as a function of time
outfile   = args[2] #output file name
tmax = as.numeric(args[3])
nrep = as.numeric(args[4])

#load joint posterior density
pext  = scan(prextfile)

#plot
t = 1:tmax
pdf(outfile,h=3.5,w=3.5)
par(mfrow=c(1,1),family="Times",las=1)
plot(t, pext/nrep ,type="l",ylim=c(0,1),xlab="year",ylab="Probability of extinction",lwd=2)
dev.off()

print(paste("plotted to file",outfile))
