## simulate trajectories ##
## simulate trajectories ##

simtrajobs <- function(e,c,p,n,p0,J,M,A,t = 17,nr = 1){
  pl  = array(0,c(nr,t,n))
  for(j in 1:nr){
      pEi = e/A
      ptab     = matrix(0,t,n)
      ptab[1,] = rbinom(n,1,p0)
      if(sum(p0==-1)>0) ptab[1,p0==-1] = rbinom(sum(p0==-1),size=1,prob = mean(p0[p0!=-1]))
      pt       = ptab[1,]
      #print(ptab)
      for( i in 1:(t-1)){
        pti = pt*rbinom(n,size=1,prob=1-pEi) #extinction phase
        pCi = sapply(c*M%*%matrix(pti,n,1), function(x) min(x,1) ) #colonization phase
        pt  = pti + (1-pti)*rbinom(n,size=1,prob=pCi)
        ptab[i+1,] = pt
        #print(ptab)
      }
      print(ptab)
      
      pl[j,,] = t(sapply(1:t, function(i) rbinom(n,J[i,],p*ptab[i,]) ))
  }
  return(pl)
}

# FUNCTION
set.seed(1)
simsim <- function(Jt,et,ct,pt,p0t,a= 1/400,nrep=100,tmax=20,suff){
  nt  = ncol(Jt)
  
  dt = sapply(1:nt, function(i){
    if(i<nt) end = 1:(nt-i)
    else end = c()
    c((i-1):0,end) } )*200 # distance matrix in meters
  At = rep(1,nt) #area of each patch
  Mmetat = exp(-a*dt)
  diag(Mmetat)=0

  plt  = simtrajobs(et,ct,pt,nt,p0t,Jt,Mmetat,At,t=tmax,nr = nrep)
  rest = sapply(1:nrep, function(i) rowMeans(plt[i,,]>0) )

  for(k in 1:nrep) write(t(plt[k,,]),file = paste("Y_",suff,k,".txt",sep=""),ncol=nt)
}

impMC = read.table("imputeMC_patches.txt")
impDC = read.table("imputeDC_patches.txt")

simsim(as.matrix(JMC),et=0.11,ct=0.12,pt=0.75,p0t=as.numeric(impMC[1,]),tmax=20,suff="simMC" )
simsim(as.matrix(JDC),et=0.41,ct=0.50,pt=0.75,p0t=as.numeric(impDC[1,]),tmax=19,suff="simDC" )

#### small test
Jt  = matrix( rbinom(nt*20,2,0.8) ,20, nt)
et  = 0.2
ct  = 0.1
pt  = 0.75
nt  = 5
p0t = c(1,1,1,1,1)
a   = 1/400
nrep=100

dt = sapply(1:nt, function(i){
  if(i<nt) end = 1:(nt-i)
  else end = c()
  c((i-1):0,end) } )*200 # distance matrix in meters
At = rep(1,nt) #area of each patch
Mmetat = exp(-a*dt)
diag(Mmetat)=0

plt  = simtrajobs(et,ct,pt,nt,p0t,Jt,Mmetat,At,t = 20,nr = nrep)
rest = sapply(1:nrep, function(i) rowMeans(plt[i,,]>0) )

for(k in 1:nrep) write(t(plt[k,,]),file = paste("test/Y_e02c01p075_",k,".txt",sep=""),ncol=nt)
write(t(Jt),file = paste("test/Jt.txt",sep=""),ncol=nt)



## MC 
piobsMC = matrix(scan("SI/run_CRLF/obsMC.txt"),ncol=11,byrow=T)
eMC = 0.18
cMC = 0.16  
nMC = 11
p0MC = piobsMC[1,]
a=1/400
nrep=1000

dMC = sapply(1:nMC, function(i){
  if(i<nMC) end = 1:(nMC-i)
  else end = c()
  c((i-1):0,end) } )*200 # distance matrix in meters
AMC = rep(1,nMC) #area of each patch
MmetaMC = exp(-a*dMC)
diag(MmetaMC)=0

plMC  = simtraj(eMC,cMC,nMC,piobsMC,MmetaMC,AMC,t = 19,nr = nrep)
resMC = sapply(1:nrep, function(i) colMeans(plMC[i,,]) )

pdf("trajDC.pdf",h=3.5,w=4)
par(las=1,family="Times")
plot(1998:2015, rowMeans(resDC) ,type="l",ylim=c(0,1),lwd=2,col=2,xlab="years",ylab="occupancy",main="Deer Creek")
q005DC = apply(resDC,1,quantile,0.05)
q095DC = apply(resDC,1,quantile,0.95)
polygon(c(1998:2015,2015:1998) , c(q005DC,rev(q095DC)) ,col=rgb(1,0,0,0.5),border=NA)
#lines(1998:2015,obsDC,lwd=2)
dev.off()

for(k in 1:nrep) write(t(plMC[k,,]),file = paste("plMLMCmiss",k,".txt",sep=""),ncol=nMC)

plMCmissl  = lapply(seq(0,0.5,0.1), function(x) simtraj(eMC,cMC,nMC,piobsMC,MmetaMC,AMC,t = 19,nr = nrep,miss = x) )

for(l in 1:length(plMCmissl) ){ for(k in 1:nrep) write(t(plMCmissl[[l]][k,,]),file = paste("plMLMCmiss",l-1,"_",k,".txt",sep=""),ncol=nMC)}

## DC 
piobsDC = matrix(scan("SI/run_CRLF/obsDC.txt"),ncol=10,byrow=T)
eDC = 0.51
cDC = 0.63  
nDC = 10
p0DC = piobsDC[1,]
a=1/400
nrep=1000

dDC = sapply(1:nDC, function(i){
  if(i<nDC) end = 1:(nDC-i)
  else end = c()
  c((i-1):0,end) } )*200 # distance matrix in meters
ADC = rep(1,nDC) #area of each patch
MmetaDC = exp(-a*dDC)
diag(MmetaDC)=0

plDC  = simtraj(eDC,cDC,nDC,piobsDC,MmetaDC,ADC,t = 18,nr = nrep)
resDC = sapply(1:nrep, function(i) colMeans(plDC[i,,]) )

pdf("trajDC.pdf",h=3.5,w=4)
par(las=1,family="Times")
plot(1998:2015, rowMeans(resDC) ,type="l",ylim=c(0,1),lwd=2,col=2,xlab="years",ylab="occupancy",main="Deer Creek")
q005DC = apply(resDC,1,quantile,0.05)
q095DC = apply(resDC,1,quantile,0.95)
polygon(c(1998:2015,2015:1998) , c(q005DC,rev(q095DC)) ,col=rgb(1,0,0,0.5),border=NA)
#lines(1998:2015,obsDC,lwd=2)
dev.off()

for(k in 1:nrep) write(t(plDC[k,,]),file = paste("plMLDCmiss",k,".txt",sep=""),ncol=nDC)

plDCmissl  = lapply(seq(0,0.5,0.1), function(x) simtraj(eDC,cDC,nDC,piobsDC,MmetaDC,ADC,t = 18,nr = nrep,miss = x) )

for(l in 1:length(plDCmissl) ){ for(k in 1:nrep) write(t(plDCmissl[[l]][k,,]),file = paste("plMLDCmiss",l-1,"_",k,".txt",sep=""),ncol=nDC)}

nrep2 = 1000
plDCdetlb  = lapply(detpc, function(x) simtraj(eDC,cDC,nDC,piobsDC,MmetaDC,ADC,t = 18,nr = nrep2,miss = -1) )
plMCdetlb  = lapply(detpc, function(x) simtraj(eMC,cMC,nMC,piobsMC,MmetaMC,AMC,t = 19,nr = nrep2,miss = -1) )
plSFCdetlb = lapply(detpc, function(x) simtraj(eSFC,cSFC,nSFC,piobsSFC,MmetaSFC,ASFC,t = 19,nr = nrep2,miss = 0) )

##
plSFCdetl  = plSFCdetlb
# put some 0s instead of 1s
for(i in 1:length(plSFCdetl)){
  for(j in 1:nrep2){
    for(k in nrow(plSFCdetl[[i]][j,,]):1){
      plSFCdetl[[i]][j,k,][plSFCdetl[[i]][j,k,]==1] = rbinom(sum(plSFCdetl[[i]][j,k,]==1),1,detpc[i]/100) 
    }
  }
}

##
plMCdetl  = plMCdetlb
# put some 0s instead of 1s
for(i in 1:length(plMCdetl)){
  for(j in 1:nrep2){
    for(k in nrow(plMCdetl[[i]][j,,]):1){
      plMCdetl[[i]][j,k,][plMCdetl[[i]][j,k,]==1] = rbinom(sum(plMCdetl[[i]][j,k,]==1),1,detpc[i]/100) 
    }
  }
}
plDCdetl  = plDCdetlb
# put some 0s instead of 1s
for(i in 1:length(plDCdetl)){
  for(j in 1:nrep2){
    for(k in nrow(plDCdetl[[i]][j,,]):1){
      plDCdetl[[i]][j,k,][plDCdetl[[i]][j,k,]==1] = rbinom(sum(plDCdetl[[i]][j,k,]==1),1,detpc[i]/100) 
    }
  }
}

zombiesMC  = sapply(1:9, function(j) sapply(1:nrep2, function(i) sum(((rowSums(plMCdetl[[j]][i,-19,]==0)==nMC) & (rowSums(plMCdetl[[j]][i,-1,]==1)>0))>0  )) ) 
zombiesDC  = sapply(1:9, function(j) sapply(1:nrep2, function(i) sum(((rowSums(plDCdetl[[j]][i,-18,]==0)==nDC) & (rowSums(plDCdetl[[j]][i,-1,]==1)>0))>0  )) ) 
zombiesSFC = sapply(1:9, function(j) sapply(1:nrep2, function(i) sum(((rowSums(plSFCdetl[[j]][i,-19,]==0)==nSFC) & (rowSums(plSFCdetl[[j]][i,-1,]==1)>0))>0  )) ) 

j=9
i=1
zombies2MC  = sapply(1:9, function(j) sapply(1:nrep2, function(i) sum( (rowSums(( (plMCdetl[[j]][i,-19,]==1))&( (plMCdetl[[j]][i,-1,]==1)|(plMCdetl[[j]][i,-1,]==-1) ) )==0)&(rowSums( plMCdetl[[j]][i,-1,]==0 )!=nMC ) ) ) )
zombies2DC  = sapply(1:9, function(j) sapply(1:nrep2, function(i) sum( (rowSums(( (plDCdetl[[j]][i,-18,]==1)|(plDCdetl[[j]][i,-18,]==-1) )&( (plDCdetl[[j]][i,-1,]==1)|(plDCdetl[[j]][i,-1,]==-1) ) )==0)&(rowSums( plDCdetl[[j]][i,-1,]==0 )!=nDC ) ) ) )
zombies2SFC = sapply(1:9, function(j) sapply(1:nrep2, function(i) sum( (rowSums(( (plSFCdetl[[j]][i,-19,]==1)|(plSFCdetl[[j]][i,-19,]==-1) )&( (plSFCdetl[[j]][i,-1,]==1)|(plSFCdetl[[j]][i,-1,]==-1) ) )==0)&(rowSums( plSFCdetl[[j]][i,-1,]==0 )!=nSFC ) ) ) )


svg("prop_null.svg",h=2.5*3,w=2.5)
par(mfrow=c(3,1),family="Times",las=1)
plot((detpc), colMeans(zombiesMC>0)*100 ,col=1,xlim=c(100,60),ylim=c(0,25),xlab="Probability of detection",pch=16,type="o",ylab="Proportion of simulations with a null likelihood (%)")
plot((detpc), colMeans(zombiesDC>0)*100 ,col=1,xlim=c(100,60),ylim=c(0,25),xlab="Probability of detection",pch=16,type="o",ylab="Proportion of simulations with a null likelihood (%)")
plot((detpc), colMeans(zombiesSFC>0)*100 ,col=1,xlim=c(100,60),ylim=c(0,25),xlab="Probability of detection",pch=16,type="o",ylab="Proportion of simulations with a null likelihood (%)")
dev.off()

plot((detpc), colMeans(zombies2MC>0) ,col=1,xlim=c(100,60),ylim=c(0,1),xlab="Probability of detection",pch=16,type="o")
points((detpc), colMeans(zombies2DC>0) ,col=2,pch=16,type="o")
points((detpc), colMeans(zombies2SFC>0) ,col=3,pch=16,type="o")

# SFC
piobsSFC = matrix(scan("obsSFC2.txt"),ncol=26,byrow=T)
eSFC = 0.51
cSFC = 0.29  
nSFC = 26
p0SFC = piobsSFC[1,]
a=1/400
nrep=100

dSFC = sapply(1:nSFC, function(i){
  if(i<nSFC) end = 1:(nSFC-i)
  else end = c()
  c((i-1):0,end) } )*250 # distance matrix in meters
ASFC = rep(1,nSFC) #area of each patch
MmetaSFC = exp(-a*dSFC)
diag(MmetaSFC)=0

plSFC  = simtraj(eSFC,cSFC,nSFC,piobsSFC,MmetaSFC,ASFC,t = 12,nr = nrep)
resSFC = sapply(1:nrep, function(i) rowMeans(plSFC[i,,]) )

for(k in 1:nrep) write(t(plSFC[k,,]),file = paste("plMLSFCmiss",k,".txt",sep=""),ncol=nSFC)

plSFCmissl  = lapply(seq(0,0.5,0.1), function(x) simtraj(eSFC,cSFC,nSFC,piobsSFC,MmetaSFC,ASFC,t = 12,nr = nrep,miss = x) )

for(l in 1:length(plSFCmissl) ){ for(k in 1:nrep) write(t(plSFCmissl[[l]][k,,]),file = paste("plMLSFCmiss",l-1,"_",k,".txt",sep=""),ncol=nSFC)}

detpc     = seq(100,60,-5)

errorsSFC = sapply(1:9, function(j) sapply(1:100, function(i) sum( (rowSums(plSFCdetl[[j]][i,-12,]*plSFCdetl[[j]][i,-1,])==0)&(rowSums(plSFCdetl[[j]][i,-1,])!=0) )  ) )
## check percentage of impossible


barplot( colMeans( errorsSFC>0 ) )

for(l in 1:length(plSFCdetl) ){ for(k in 1:nrep) write(t(plSFCdetl[[l]][k,,]),file = paste("plMLSFCdet",detpc[l],"_",k,".txt",sep=""),ncol=nSFC)}

############ simul likelihood #####################
#function for H1
simtrajH1 <- function(e,c,n,p0,M,K,ts,td,nr){
  pl  = array(0,c(nr,n))
  for(j in 1:nr){
    pEi = e/K
    if(pEi>1) pEi = 1
    pt       = p0
    for( i in 1:(ts-1)){
      pti = pt*rbinom(n,size=1,prob=1-pEi) #extinction phase
      pCi = sapply(c*M%*%matrix(pti,n,1)*K, function(x) min(x,1) ) #colonization phase
      pt  = pti + (1-pti)*rbinom(n,size=1,prob=pCi)
      if(sum(pt)==0){
        pl[j,] = matrix(0,1,n)
        break
      }
    }
    pEi = e
    if(pEi>1) pEi = 1
    for( i in 1:(td)){
      pti = pt*rbinom(n,size=1,prob=1-pEi) #extinction phase
      pCi = sapply(c*M%*%matrix(pti,n,1), function(x) min(x,1) ) #colonization phase
      pt  = pti + (1-pti)*rbinom(n,size=1,prob=pCi)
      if(sum(pt)==0){
        pl[j,] = matrix(0,1,n)
        break
      }
    }
    pl[j,]=pt
  }
  return(pl)
}

simtrajH2 <- function(e,c,n,p0,Mm,K,dd,a,ts,td,nr){
  M = cbind(Mm, exp(-dd*1:n/a)*K )
  M = rbind(M, c(exp(-dd*1:n/a)*K,0) )
  pl  = array(0,c(nr,n))
  for(j in 1:nr){
    pEi = e
    pt  = p0
    for( i in 1:(ts-1)){
      pti = pt*rbinom(n,size=1,prob=1-pEi) #extinction phase
      pCi = sapply(c*M%*%matrix(c(pti,1),n+1,1), function(x) min(x,1) )[-(n+1)] #colonization phase
      pt  = pti + (1-pti)*rbinom(n,size=1,prob=pCi)
      if(sum(pt)==0){
        pl[j,] = matrix(0,1,n)
        break
      }
    }
    pEi = e
    if(pEi>1) pEi = 1
    for( i in 1:(td)){
      pti = pt*rbinom(n,size=1,prob=1-pEi) #extinction phase
      pCi = sapply(c*Mm%*%matrix(pti,n,1), function(x) min(x,1) ) #colonization phase
      pt  = pti + (1-pti)*rbinom(n,size=1,prob=pCi)
      if(sum(pt)==0){
        pl[j,] = matrix(0,1,n)
        break
      }
    }
    pl[j,]=pt
  }
  return(pl)
}

#####
# test MC
pendMC = obsMC[,1]
ppendMC = sum(pendMC)

ns=1000
simpMC = matrix(-1,length(Kl),ns)
simpMC2 = matrix(-1,length(Kl),ns)
for(i in 1:length(Kl)){
  if(i%%100==0) print(i)
  p0MC = rbinom(nMC,size = 1,0.5)
  plMC  = rowSums(simtrajH1(eMC,cMC,nMC,p0MC,Mmeta,Kl[i],20,26,ns))
  plMC2  = rowSums(simtrajH1(eMC,cMC,nMC,p0MC,Mmeta,Kl[i],20,1,ns))
  simpMC[i,] = plMC
  simpMC2[i,] = plMC2
}
lsimpMC = rowMeans(simpMC==ppendMC)
lsimpMC2 = sapply(1:200, function(i) mean(lsimpMC[((i-1)*5+1):(i*5)] ) )

plot(Kl, lsimpMC ,type="l",log="x")
lines(Kl,rowMeans(simpMC)/max(rowMeans(simpMC))*mean(lsimpMC2[150:200]),col=4 )
lines(Kl,rowMeans(simpMC2)/max(rowMeans(simpMC2))*mean(lsimpMC2[150:200]),col=5 )

lines(1:200/2,lsimpMC2,col=2)
lines(Kl,LikMCdis/max(LikMCdis)*0.065,col=3)


####
eSFC = 0.51
cSFC = 0.29
nSFC = 26
pendSFC = obsSFC[,1]
ppendSFC = sum(pendSFC)

ns=1000
simpSFC = matrix(-1,length(Kl),ns)
for(i in 1:length(Kl)){
  if(i%%100==0) print(i)
  p0SFC = rep(1,nSFC)#rbinom(nSFC,size = 1,0.5)
  plSFC  = rowSums(simtrajH1(eSFC,cSFC,nSFC,p0SFC,MmetaSFC,Kl[i],20,1,ns))
  simpSFC[i,] = plSFC
}
lsimpSFC = rowMeans(simpSFC==ppendSFC)
lsimpSFC2 = sapply(1:200, function(i) mean(lsimpSFC[((i-1)*5+1):(i*5)] ) )

plot(Kl, rowMeans(simpSFC) ,log="x",type="l")
#lines(Kl,rowMeans(simpMC)/max(rowMeans(simpMC))*mean(lsimpMC2[150:200]),col=4 )

old1 = simpSFC

plot(Kl, lsimpSFC ,type="l",log="x")
lines(1:200/2,lsimpSFC2,col=2)
#lines(Kl,LikSFCdis/max(LikSFCdis)*mean(lsimpSFC2[150:200]),col=3)

#loss 
ns=1000
simpSFC = array(-1,c(length(dl),length(Kl),ns) )
for(j in 1:length(dl)){
  print(j)
  for(i in 1:length(Kl)){
    if(i%%100==0) print(i)
    p0SFC = rep(1,nSFC)#rbinom(nSFC,size = 1,0.5)
    plSFC  = rowSums(simtrajH2(e=eSFC,c=cSFC,n=nSFC,p0 = p0SFC,Mm = MmetaSFC,K = Kl[i],dd = dl[j],a = 400, ts = 20,td = 1,nr = ns))
    simpSFC[j,i,] = plSFC
  }
}
jlsimpSFC = matrix(0,length(dl),length(Kl) )
for(j in 1:length(dl)){
  for(i in 1:length(Kl)){
    jlsimpSFC[j,i] = mean(simpSFC[j,i,])
  }
}

lsimpSFC = rowMeans(simpSFC==ppendSFC)
lsimpSFC2 = sapply(1:200, function(i) mean(lsimpSFC[((i-1)*5+1):(i*5)] ) )

# MC to compute matrix
simtrajH11 <- function(e,c,n,p0,M,K,nr){
  pl  = array(0,c(nr,n))
  for(j in 1:nr){
    pEi = e/K
    if(pEi>1) pEi = 1
    pti = p0*rbinom(n,size=1,prob=1-pEi) #extinction phase
    pCi = sapply(c*M%*%matrix(pti,n,1)*K, function(x) min(x,1) ) #colonization phase
    pl[j,]=pti + (1-pti)*rbinom(n,size=1,prob=pCi)
  }
  return(pl)
}

simtrajH21 <- function(e,c,n,p0,Mm,K,dd,a,nr){
  M = cbind(Mm, exp(-dd*1:n/a)*K )
  M = rbind(M, c(exp(-dd*1:n/a)*K,0) )
  pl  = array(0,c(nr,n))
  for(j in 1:nr){
    pEi = e
    pti = p0*rbinom(n,size=1,prob=1-pEi) #extinction phase
    pCi = sapply(c*M%*%matrix(c(pti,1),n+1,1), function(x) min(x,1) )[-(n+1)] #colonization phase
    pl[j,]= pti + (1-pti)*rbinom(n,size=1,prob=pCi)
  }
  return(pl)
}

## compute relative contribution of states:
ns2 = 1000

listp0MC = vector("list",length = nMC+1)
for( i in 1:nMC ){
  for(j in 1:ns2){
    p0MC = rep(0,nMC)
    p0MC[sample(1:nMC,i)] = 1
    ll = simtrajH11(e=eMC,c=cMC,n=nMC,p0 = p0MC,M = Mmeta,K = 1,nr = 1) 
    ll = simtrajH11(e=eMC,c=cMC,n=nMC,p0 = ll,M = Mmeta,K = 1,nr = 1) 
    listp0MC[[sum(ll)+1]] = rbind(listp0MC[[sum(ll)+1]], ll)
  }
}
##


ns=10000
P0 = matrix(0,nMC+1,nMC+1)
P0[1,1] = 1
for(i in 1:nMC){
  for(j in 1:ns){
    p0MC = listp0MC[[i+1]][sample(1:nrow(listp0MC[[i+1]]),1),]
    ll = rowSums( simtrajH11(e=eMC,c=cMC,n=nMC,p0 = p0MC,M = Mmeta,K = 1,nr = 1) )
    for(l in ll) P0[i+1,l+1] = P0[i+1,l+1]+1/ns
  }
}
Pt01 = P0%^%26
Pt02 = P0%^%38

#P0old = P0

# MC 2 

ns1=100
MCMC = rep(-1,length(Kl))
resMCl = vector("list", length(Kl))
resMCl2 = vector("list", nMC+1) # regroup res into classes of different # occupied segments
for(i in 1:length(Kl)){
  if(i%%100==0) print(i)
  for(j in 1:ns1){
    p0MC = rbinom(nMC,size = 1, 0.5) #sample(c(rep(1,k),rep(0,nMC-k))) #listp0MC[[k+1]][sample(1:nrow(listp0MC[[k+1]]),1),]
    l = simtrajH11(e=eMC,c=cMC,n=nMC,p0 = p0MC,M = Mmeta,K = Kl[i],nr = 1)
    for(tt in 1:19) l = simtrajH11(e=eMC,c=cMC,n=nMC,p0 = l,M = Mmeta,K = Kl[i],nr = 1)
    resMCl[[i]] = rbind(resMCl[[i]],l)
    resMCl2[[sum(l)+1]] = rbind( resMCl2[[sum(l)+1]], l)
  }
}

ns2=10000
lMC0 = rep(0,nMC+1)
for(i in 1:nMC){
  print(i)
  for(j in 1:ns2){
    p0MC = resMCl2[[i+1]][sample( 1:nrow(resMCl2[[i+1]]) ,1),]
    if(sum(p0MC) > 0 ){ 
      l = simtrajH11(e=eMC,c=cMC,n=nMC,p0 = p0MC,M = Mmeta,K = 1,nr = 1)
      for(tt in 1:25) l = simtrajH11(e=eMC,c=cMC,n=nMC,p0 = l,M = Mmeta,K = 1,nr = 1)
      #resMCl2[[i]] = rbind(resMCl2[[i]],l)
      ll = sum(l)
    }
    else ll = 0
    #print(c(i,ll))
    lMC0[i+1] = lMC0[i+1] + 1*(ll==9)
  }
}

ns2=10000
lMC02 = rep(0,nMC+1)
for(i in 1:nMC){
  print(i)
  for(j in 1:ns2){
    p0MC = resMCl2[[i+1]][sample( 1:nrow(resMCl2[[i+1]]) ,1),]
    if(sum(p0MC) > 0 ){ 
      l = simtrajH11(e=eMC,c=cMC,n=nMC,p0 = p0MC,M = Mmeta,K = 1,nr = 1)
      for(tt in 1:37) l = simtrajH11(e=eMC,c=cMC,n=nMC,p0 = l,M = Mmeta,K = 1,nr = 1)
      ll = sum(l)
    }
    else ll = 0
    lMC02[i+1] = lMC02[i+1] + 1*(ll==9)
  }
}


la = sapply(1:1000, function(j) sum( sapply(0:nMC, function(i) mean(rowSums( resMCl[[j]] )==i ) )*lMC0/ns2 ) )

restest = c()
for(i in 1:ns2){
  ltest = simtrajH11(e=eMC,c=cMC,n=nMC,p0 = rep(1,nMC),M = Mmeta,K = 1,nr = 1)
  for(tt in 1:37) ltest = rbind(ltest, simtrajH11(e=eMC,c=cMC,n=nMC,p0 = ltest[tt,],M = Mmeta,K = 1,nr = 1) )
  restest = rbind(restest,rowSums( ltest ))
}

## H2 
ns1=100
MCMCloss = rep(-1,length(Kl))
resMClossl = vector("list", length(dl))
for(i in 1:length(dl)) resMClossl[[i]] = vector("list", length(Kl))
#resMClossl2 = vector("list", nMC+1) # regroup res into classes of different # occupied segments
for(o in 1:length(dl)){
  print(o)
  for(i in 1:length(Kl)){
    if(i%%100==0) print(i)
    for(j in 1:ns1){
      p0MC = rbinom(nMC,size = 1, 0.5) #sample(c(rep(1,k),rep(0,nMC-k))) #listp0MC[[k+1]][sample(1:nrow(listp0MC[[k+1]]),1),]
      l = simtrajH21(e=eMC,c=cMC,n=nMC,p0 = p0MC,Mm = Mmeta,K = Kl[i],dd=dl[o],a = 400,nr = 1)
      for(tt in 1:19) l = simtrajH21(e=eMC,c=cMC,n=nMC,p0 = l,Mm = Mmeta,K = Kl[i],dd=dl[o],a = 400,nr = 1)
      resMClossl[[o]][[i]] = rbind(resMClossl[[o]][[i]],l)
      #resMClossl2[[sum(l)+1]] = rbind( resMClossl2[[sum(l)+1]], l)
    }
  }
}

laloss = matrix(0,length(dl),length(Kl))
for(o in 1:length(dl)){
  for(i in 1:length(Kl)){
    laloss[o,i] = sum( sapply(0:nMC, function(j) mean(rowSums( resMClossl[[o]][[i]] )==j ) )*lMC02/ns2 )
  }
}

sapply(1:200, function(x) mean(rowSums(resMCl[[x]])==9 ) )

##

ns=100
MCMC = rep(-1,length(Kl))
resMCl = vector("list", nMC)
for(i in 1:length(Kl)){
    if(i%%100==0) print(i)
    P11 = matrix(0,nMC+1,nMC+1)
    P11[1,1] = 1
    for(k in 1:nMC){
      for(j in 1:ns){
        p0MC = sample(c(rep(1,k),rep(0,nMC-k))) #listp0MC[[k+1]][sample(1:nrow(listp0MC[[k+1]]),1),]
        l = simtrajH11(e=eMC,c=cMC,n=nMC,p0 = p0MC,M = Mmeta,K = Kl[i],nr = 1)
        for(tt in 1:19) l = simtrajH11(e=eMC,c=cMC,n=nMC,p0 = l,M = Mmeta,K = Kl[i],nr = 1)
        resMCl[[k]] = rbind(resMCl[[k]],l)
        ll = sum(l)
        P11[k+1,ll+1] = P11[k+1,ll+1]+1/ns
      }
    }
    #Pt11 = P11%^%20
    Pf = P11%*%Pt01
    
    MCMC[i] = mean( Pf[,10] )
}

plot( Kl, MCMC ,log="x")
lines(Kl,LikMCdis/max(LikMCdis)*max(MCMC),col=3)


ns=100
MCMCloss = matrix(-1,length(dl),length(Kl))
for(j in 1:length(dl)){
  print(j)
  for(i in 1:length(Kl)){
    if(i%%100==0) print(i)
    P11 = matrix(0,nMC+1,nMC+1)
    P11[1,1] = 1
    for(k in 1:nMC){
      p0MC = listp0MC[[k+1]][sample(1:nrow(listp0MC[[k+1]]),1),]
      for(o in 1:ns){
        ll = sum( simtrajH21(e=eMC,c=cMC,n=nMC,p0 = p0MC,Mm = Mmeta,K = Kl[i],dd=dl[j],a = 400,nr = 1) )
        P11[k+1,ll+1] = P11[k+1,ll+1]+1/ns
      }
    }
  Pt11 = P11%^%20
  Pf = Pt11%*%Pt02
  MCMCloss[j,i] = mean( Pf[,10] )
  }
}

image(log(Kl,10),dl,t(MCMCloss),col=c(0,coco) ,xlab="",ylab="" ,axes=F)

### SFC 
ns=100000
P0SFC = matrix(0,nSFC+1,nSFC+1)
P0SFC[1,1] = 1
for(i in 1:nSFC){
  for(j in 1:ns){
    p0SFC = sample(c(rep(1,i),rep(0,nSFC-i))  ) #listp0SFC[[i+1]][sample(1:nrow(listp0SFC[[i+1]]),1),]
    l = sum( simtrajH11(e=eSFC,c=cSFC,n=nSFC,p0 = p0SFC,M = MmetaSFC,K = 1,nr = 1) )
    P0SFC[i+1,l+1] = P0SFC[i+1,l+1]+1/ns
  }
}
Pt01SFC = P0SFC%^%26
Pt02SFC = P0SFC%^%38

ns=1000
MCSFC = matrix(-1,1,length(Kl))
for(i in 1:length(Kl)){
    if(i%%100==0) print(i)
    P11 = matrix(0,nSFC+1,nSFC+1)
    P11[1,1] = 1
    for(k in 1:nSFC){
      for(ii in 1:ns){
        p0SFC = rep(0,nSFC)
        p0SFC[sample(1:nSFC,k)] = 1
        l = sum( simtrajH11(e=eSFC,c=cSFC,n=nSFC,p0 = p0SFC,M = MmetaSFC,K = Kl[i],nr = 1) )
        P11[k+1,l+1] = P11[k+1,l+1]+1/ns
      }
    }
    Pt11 = P11%^%20
    Pf = Pt11%*%Pt01SFC
    MCSFC[i] = mean( Pf[,12] )
}


ns=1000
MCSFCloss = matrix(-1,length(dl),length(Kl))
for(j in 1:length(dl)){
  print(j)
  for(i in 1:length(Kl)){
    if(i%%100==0) print(i)
    P11 = matrix(0,nSFC+1,nSFC+1)
    P11[1,1] = 1
    for(k in 1:nSFC){
      for(ii in 1:ns){
        p0SFC = rep(0,nSFC)
        p0SFC[sample(1:nSFC,k)] = 1
        l = sum( simtrajH21(e=eSFC,c=cSFC,n=nSFC,p0 = p0SFC,Mm = MmetaSFC,K = Kl[i],dd=dl[j],a = 400,nr = 1) )
        P11[k+1,l+1] = P11[k+1,l+1]+1/ns
      }
    }
    Pt11 = P11%^%20
    Pf = Pt11%*%Pt02SFC
    MCSFCloss[j,i] = mean( Pf[,12] )
  }
}

png("KlossestSFC.png",width=250,height=250,bg=rgb(0,0,0,0))#h=3.5,w=3.5*2)
par(mfrow=c(1,1),family="Times",las=1)
image(log(Kl,10),dl,t(MCSFCloss)/max(MCSFCloss),col=c(0,coco) ,xlab="",ylab="" ,axes=F)
axis(2,at=seq(200,2000,400),labels=rep("",5))
axis(1,at=seq(-1,2,1),labels=rep("",4))#labels=c("0.1","1","10","100"))
dev.off()


MCSFC2 = sapply(0:99, function(x) mean(MCSFC[(x*10+1):((x+1)*10)]))
pdf("Kdisest2.pdf",h=3,w=3*3)
par(mfrow=c(1,3),family="Times",las=1)
Kl = seq(0.1,100,0.1)
plot(0.0001,-1,type="l",log="x",xlab=expression(italic(K[dis])),ylab="probability distribution",col=2,ylim=c(0,0.6),lwd=2,xlim=c(0.1,100))
polygon(c(Kl,rev(Kl)), c(LikMCdis[1,]/maMC/1.885,rep(0,1000) ),col=rgb(1,0,0,0.5),border=NA ) 
polygon(c(0.1,0.1,100,100),c(0,1,1,0)/(log(100,10)-log(0.1,10)),lwd=2,col=rgb(0,0,0,0.4),border=NA)
abline(v = Kl[idmaMC],col=2,lty=2,lwd=2)
plot(0.0001,-1,type="l",log="x",xlab=expression(italic(K[dis])),ylab="",col=2,ylim=c(0,0.6),lwd=2,xlim=c(0.1,100))
polygon(c(Kl,rev(Kl)), c(LikDCdis[1,]/maDC/1.885,rep(0,1000) ),col=rgb(1,0,0,0.5),border=NA ) 
polygon(c(0.1,0.1,100,100),c(0,1,1,0)/(log(100,10)-log(0.1,10)),lwd=2,col=rgb(0,0,0,0.4),border=NA)
abline(v = Kl[idmaDC],col=2,lty=2,lwd=2)

plot(0.0001,-1,type="l",log="x",xlab=expression(italic(K[dis])),ylab="",col=2,ylim=c(0,4),lwd=2,xlim=c(0.1,100))
polygon(c(0:99+0.1,rev(0:99+0.1)), c(MCSFC2/max(MCSFC2)/0.26,rep(0,100) ),col=rgb(1,0,0,0.5),border=NA ) 
polygon(c(0.1,0.1,100,100),c(0,1,1,0)/(log(100,10)-log(0.1,10)),lwd=2,col=rgb(0,0,0,0.4),border=NA)
abline(v = Kl[950],col=2,lty=2,lwd=2)
dev.off()
