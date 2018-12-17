library("pracma")
Sigma=read.table("Sigmany.txt")
Sigma=as.matrix(Sigma)
R=read.table("Cny.txt")
R=as.matrix(R)
X=read.table("D.txt")
X=as.matrix(X)
effekt=read.table("effekt.txt")
effekt=as.matrix(effekt)
y=c(effekt/400,400)
invsqrtSigma=solve(sqrt(Sigma))
ytilde=invsqrtSigma%*%y
Xtilde=invsqrtSigma%*%X
XTX=t(Xtilde)%*%Xtilde
Nulm=matrix(rep(0,676),26,26)
A=matrix(c(rbind(XTX,R),rbind(t(R),Nulm)),63,63)
Ainv=solve(A)
Beta=Ainv%*%c(t(Xtilde)%*%ytilde,rep(0,26))
Beta[1:37]#Her er det første estimat for strømstyrker og spændinger i hele netværket:
for(i in 1:50)
{
  y=c(effekt/c(Beta[21:22],Beta[25],Beta[27],Beta[29],Beta[31],Beta[33],Beta[35:37]),400)
  Sigma=diag((y*0.01)^2)
  invsqrtSigma=solve(sqrt(Sigma))
  ytilde=invsqrtSigma%*%y
  Xtilde=invsqrtSigma%*%X
  XTX=t(Xtilde)%*%Xtilde
  A=matrix(c(rbind(XTX,R),rbind(t(R),Nulm)),63,63)
  Ainv=solve(A)
  Beta=Ainv%*%c(t(Xtilde)%*%ytilde,rep(0,26))
}
y=rnorm(37,mean=Beta[1:37],sd=abs(Beta[1:37])*0.01)
Sigma=diag((c(rep(6.3,18),rep(400,19))*0.01)^2)
invsqrtSigma=solve(sqrt(Sigma))
X=diag(rep(1,37))
Xtilde=invsqrtSigma%*%X
XTX=t(Xtilde)%*%Xtilde
Nulm=matrix(rep(0,676),26,26)
A=matrix(c(rbind(XTX,R),rbind(t(R),Nulm)),63,63)
Ainv=solve(A)
P=Xtilde%*%Ainv[1:37,1:37]%*%t(Xtilde)
lev=c()
for(i in 1:37)
{
  lev=c(lev,P[i,i])
}
levord=sort(lev)
levordindex=sort(lev,index.return=TRUE)[2]
levordindex=as.vector(levordindex[[1]])
flag=0
v=rep(1,length(lev))
for(i in 1:37)
{
  flag=0
  
for(j in 1:37)
{
  if(flag==0)
  {
    remove=levordindex[i]
    Sigmat=Sigma[-remove,-remove]
    Xt=X[-remove,]
    vt=v[-remove]
    yt=y[-remove]
    invsqrtSigmat=solve(sqrt(Sigmat))
    Xtildet=invsqrtSigmat%*%Xt
    XTXt=t(Xtildet)%*%Xtildet
    A=matrix(c(rbind(XTXt,R),rbind(t(R),Nulm)),63,63)
    Ainv=solve(A)
    P=Xtildet%*%Ainv[1:37,1:37]%*%t(Xtildet)
    lev=c()
    for(k in 1:length(vt))
    {
     lev=c(lev,P[k,k])
    }
    if(max(lev)<0.9999)
    {
      flag=1
      X=Xt
      Sigma=Sigmat
      v=vt
      y=yt
    }
  }
}
}#Ovenfor har vi fjernet såmange målinger som muligt, uden at nogen af leveragene bliver 1.
#Indfører målefejl, og tester om den kan detekteres med standardiserede residualer:
y[15]=450
invsqrtSigma=solve(sqrt(Sigma))
ytilde=invsqrtSigma%*%y
Xtilde=invsqrtSigma%*%X
XTX=t(Xtilde)%*%Xtilde
A=matrix(c(rbind(XTX,R),rbind(t(R),Nulm)),63,63)
Ainv=solve(A)
P=Xtilde%*%Ainv[1:37,1:37]%*%t(Xtilde)
v=c()
for(i in 1:length(y))
{
  v=c(v,(ytilde[i]-(P%*%ytilde)[i])/sqrt(1-P[i,i]))
}
v
#Vi fjerner den 15. måling, da det standardiserede residual er meget højt for denne måling. Så finder vi de standardiserede residualer igen:
y=y[-15]
X=X[-15,]
Sigma=Sigma[-15,-15]
invsqrtSigma=solve(sqrt(Sigma))
ytilde=invsqrtSigma%*%y
Xtilde=invsqrtSigma%*%X
XTX=t(Xtilde)%*%Xtilde
A=matrix(c(rbind(XTX,R),rbind(t(R),Nulm)),63,63)
Ainv=solve(A)
P=Xtilde%*%Ainv[1:37,1:37]%*%t(Xtilde)
v=c()
for(i in 1:length(y))
{
  v=c(v,(ytilde[i]-(P%*%ytilde)[i])/sqrt(1-P[i,i]))
}
v