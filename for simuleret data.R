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
Sigma=diag((y*0.01)^2)
X=diag(1,37)
invsqrtSigma=solve(sqrt(Sigma))
ytilde=invsqrtSigma%*%y
Xtilde=invsqrtSigma%*%X
XTX=t(Xtilde)%*%Xtilde
A=matrix(c(rbind(XTX,R),rbind(t(R),Nulm)),63,63)
Ainv=solve(A)
Beta2=Ainv%*%c(t(Xtilde)%*%ytilde,rep(0,26))
y
Beta[1:37]
konf=matrix(rep(0,74),37,2)
for(i in 1:37)
{
  konf[i,1]=Beta2[i]-1.959963984540*sqrt(Ainv[i,i])
  konf[i,2]=Beta2[i]+1.959963984540*sqrt(Ainv[i,i])
}
konf#Konfidensintervaller for spændingerne:
#determinationskoefficient udregnes:
P=Xtilde%*%Ainv[1:37,1:37]%*%t(Xtilde)
mu<-P%*%ytilde
V<-c()
for (i in 1:length(ytilde)){
  V<-c(V,(ytilde[i]-mu[i])^2)
}
sum(V)
V2<-c()
for (i in 1:length(ytilde)){
  V2<-c(V2,(ytilde[i]-(sum(ytilde)/length(ytilde)))^2)
}
R2<-1-(sum(V)/sum(V2))
R2
#Vi fjerner 5 målinger ved leverage metoden
varians<-c(rep((0.01*0.0155)^2,18),rep((0.01*1)^2,19))
varians
Sigma1=diag(varians)
invsqrtSigma1=solve(sqrt(Sigma1))
Xtilde=invsqrtSigma1%*%X
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
v=rep(1,37)
for(i in 1:5)
{
  if(min(lev)<0.9999999999)
  {
    remove=match(min(lev),lev)
    Sigma=Sigma[-remove,-remove]
    Sigma1=Sigma1[-remove,-remove]
    y=y[-remove]
    X=X[-remove,]
    v=v[-remove]
    invsqrtSigma1=solve(sqrt(Sigma1))
    Xtilde=invsqrtSigma1%*%X
    XTX=t(Xtilde)%*%Xtilde
    A=matrix(c(rbind(XTX,R),rbind(t(R),Nulm)),63,63)
    Ainv=solve(A)
    P=Xtilde%*%Ainv[1:37,1:37]%*%t(Xtilde)
    lev=c()
    for(j in 1:length(v))
    {
      lev=c(lev,P[j,j])
    }
    print(lev)
  }
}
invsqrtSigma=solve(sqrt(Sigma))
ytilde=invsqrtSigma%*%y
Xtilde=invsqrtSigma%*%X
XTX=t(Xtilde)%*%Xtilde
A=matrix(c(rbind(XTX,R),rbind(t(R),Nulm)),63,63)
Ainv=solve(A)
Beta3=Ainv%*%c(t(Xtilde)%*%ytilde,rep(0,26))
Beta[1:37]
konf=matrix(rep(0,74),37,2)
for(i in 1:37)
{
  konf[i,1]=Beta3[i]-1.959963984540*sqrt(Ainv[i,i])
  konf[i,2]=Beta3[i]+1.959963984540*sqrt(Ainv[i,i])
}
konf