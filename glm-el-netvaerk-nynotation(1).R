library("pracma")
Sigma=read.table("Sigmany.txt")
Sigma=as.matrix(Sigma)
#Sigma=Sigma/400^2
#Sigma[11,11]=Sigma[11,11]*400^2
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
  ytilde=invsqrtSigma%*%y
  Beta=Ainv%*%c(t(Xtilde)%*%ytilde,rep(0,26))
}
Beta[1:37]#Her har vi divideret med de estimerede spændinger for at finde strømstyrken. (Istedet for bare at dividere dem alle med 400.) Vi itererer til konvergens.
konf=matrix(rep(0,74),37,2)
for(i in 1:37)
{
  konf[i,1]=Beta[i]-1.959963984540*sqrt(Ainv[i,i])
  konf[i,2]=Beta[i]+1.959963984540*sqrt(Ainv[i,i])
}
konf#Konfidensintervaller for spændingerne:
z=c()
T=360
for(i in 19:37)
{
  z=c(z,(Beta[i]-T)/sqrt(Ainv[i,i])) #Wald's test
}
z
pnorm(z)# p-værdier udregnes. De giver alle tal større end 0.05, så vi kan for hver knude sige med 95% sikkerhed, at den knude ikke har en spænding under 360V.
# simdata=rnorm(37,mean=Beta[1:37],sd=sqrt(Beta[1:37]^2)*0.01)
# Sigma=diag((simdata*0.01)^2)
# D=diag(rep(1,37))
# invsqrtSigma=solve(sqrt(Sigma))
# dtilde=invsqrtSigma%*%simdata
# Dtilde=invsqrtSigma%*%D
# DTD=t(Dtilde)%*%Dtilde
# Nulm=matrix(rep(0,676),26,26)
# A=matrix(c(rbind(DTD,C),rbind(t(C),Nulm)),63,63)
# Ainv=solve(A)
# Beta=Ainv%*%c(t(Dtilde)%*%dtilde,rep(0,26))
# Beta[1:37]
# P=Dtilde%*%Ainv[1:37,1:37]%*%t(Dtilde)
# plot((dtilde-P%*%dtilde)[19:37])#Residualplot i forhold til indeks
# plot((sqrt(Sigma)%*%(P%*%dtilde))[1:18],(dtilde-P%*%dtilde)[1:18])#residualplot ift fittede værdier
# P%*%dtilde
# P%*%simdata
# P
# (dtilde-P%*%dtilde)#residualer på den ene måde
# dtilde-Dtilde%*%Beta[1:37]#residualer på den anden måde
# Beta[1:18]
# (sqrt(Sigma)%*%(P%*%dtilde))[1:18]
# lev=c()
# for(i in 1:37)
# {
#   lev=c(lev,P[i,i])
# }
# remove=match(min(lev),lev)
simdata=rnorm(37,mean=Beta[1:37],sd=sqrt(Beta[1:37]^2)*0.01)
Sigma=diag((Beta[1:37]*0.01)^2)
X=diag(rep(1,37))
invsqrtSigma=solve(sqrt(Sigma))
ytilde=invsqrtSigma%*%simdata
Xtilde=invsqrtSigma%*%X
XTX=t(Xtilde)%*%Xtilde
Nulm=matrix(rep(0,676),26,26)
A=matrix(c(rbind(XTX,R),rbind(t(R),Nulm)),63,63)
Ainv=solve(A)
Beta=Ainv%*%c(t(Xtilde)%*%ytilde,rep(0,26))
Beta[1:37]
P=Xtilde%*%Ainv[1:37,1:37]%*%t(Xtilde)
plot((ytilde-P%*%ytilde)[19:37])#Residualplot i forhold til indeks
plot((sqrt(Sigma)%*%(P%*%ytilde))[1:18],(ytilde-P%*%ytilde)[1:18])#residualplot ift fittede værdier
P%*%ytilde
P%*%simdata
P
(ytilde-P%*%ytilde)#residualer på den ene måde
ytilde-Xtilde%*%Beta[1:37]#residualer på den anden måde
Beta[1:18]
(sqrt(Sigma)%*%(P%*%ytilde))[1:18]
lev=c()
for(i in 1:37)
{
  lev=c(lev,P[i,i])
}
y=simdata
for(i in 1:37)
{
  if(min(lev)<0.9999999999)
  {
    remove=match(min(lev),lev)
    Sigma=Sigma[-remove,-remove]
    X=X[-remove,]
    y=y[-remove]
    invsqrtSigma=solve(sqrt(Sigma))
    Xtilde=invsqrtSigma%*%X
    XTX=t(Xtilde)%*%Xtilde
    A=matrix(c(rbind(XTX,R),rbind(t(R),Nulm)),63,63)
    Ainv=solve(A)
    P=Xtilde%*%Ainv[1:37,1:37]%*%t(Xtilde)
    lev=c()
    for(j in 1:length(y))
    {
      lev=c(lev,P[j,j])
    }
  }
}
X
# c(Beta[2:3],Beta[6],Beta[8],Beta[10],Beta[12],Beta[14],Beta[16:19])
# d=rnorm(11,mean=c(Beta[2:3],Beta[6],Beta[8],Beta[10],Beta[12],Beta[14],Beta[16:19]),sd=sqrt((c(Beta[2:3],Beta[6],Beta[8],Beta[10],Beta[12],Beta[14],Beta[16:19])*0.01)^2))
# Sigma=diag((d*0.01)^2)
# invsqrtSigma=solve(sqrt(Sigma))
# dtilde=invsqrtSigma%*%d
# Dtilde=invsqrtSigma%*%D
# DTD=t(Dtilde)%*%Dtilde
# Nulm=matrix(rep(0,676),26,26)
# A=matrix(c(rbind(DTD,C),rbind(t(C),Nulm)),63,63)
# Ainv=solve(A)
# Beta=Ainv%*%c(t(Dtilde)%*%dtilde,rep(0,26))
# dtilde-Dtilde%*%Beta[1:37]
# dtilde-Dtilde%*%Ainv[1:37,1:37]%*%t(Dtilde)%*%dtilde
# plot(1:11,dtilde-Dtilde%*%Ainv[1:37,1:37]%*%t(Dtilde)%*%dtilde)