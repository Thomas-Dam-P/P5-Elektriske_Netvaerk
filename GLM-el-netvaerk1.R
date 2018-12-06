library("pracma")
#Sigma=read.table("Sigmany400.txt")
#Sigma=as.matrix(Sigma)
#Sigma=Sigma/400^2
#Sigma[11,11]=Sigma[11,11]*400^2
C=read.table("Cny.txt")
C=as.matrix(C)
D=read.table("D.txt")
D=as.matrix(D)
effekt=read.table("effekt.txt")
effekt=as.matrix(effekt)
d=c(effekt/400,400)
Sigma=diag((d*0.01)^2)
invsqrtSigma=solve(sqrt(Sigma))
dtilde=invsqrtSigma%*%d
Dtilde=invsqrtSigma%*%D
DTD=t(Dtilde)%*%Dtilde
Nulm=matrix(rep(0,676),26,26)
A=matrix(c(rbind(DTD,C),rbind(t(C),Nulm)),63,63)
Ainv=solve(A)
Beta=Ainv%*%c(t(Dtilde)%*%dtilde,rep(0,26))
Beta[1:37]#Her er det første estimat for strømstyrker og spændinger i hele netværket:
for(i in 1:500)
{
  d=c((effekt)/c(Beta[21:22],Beta[25],Beta[27],Beta[29],Beta[31],Beta[33],Beta[35:37]),400)
  Sigma=diag((d*0.01)^2)
  invsqrtSigma=solve(sqrt(Sigma))
  dtilde=invsqrtSigma%*%d
  Dtilde=invsqrtSigma%*%D
  DTD=t(Dtilde)%*%Dtilde
  A=matrix(c(rbind(DTD,C),rbind(t(C),Nulm)),63,63)
  Ainv=solve(A)
  Beta=Ainv%*%c(t(Dtilde)%*%dtilde,rep(0,26))
}
Beta[1:37]#Her har vi divideret med de estimerede spændinger for at finde strømstyrken. (Istedet for bare at dividere dem alle med 400.) Vi itererer til konvergens.
konf=matrix(rep(0,74),37,2)
for(i in 1:37)
{
  konf[i,1]=Beta[i]-1.959963984540*sqrt(Ainv[i,i])
  konf[i,2]=Beta[i]+1.959963984540*sqrt(Ainv[i,i])
}
konf#Konfidensintervaller for parametrene:
z=c()#H_0i: spændingen i knude i er mindre end T volt:
T=390
for(i in 19:37)
{
  z=c(z,(Beta[i]-T)/sqrt(Ainv[i,i])) #Wald's test
}
z
1-pnorm(z)# p-værdier udregnes. hvis p-værdien er under 0.05, så kan vi med 95% sikkerhed sige, at den knude ikke har spænding lavere end T (forkaste H_0i).


simdata=rnorm(37,mean=Beta[1:37],sd=abs(Beta[1:37])*0.01)
Sigma=diag((Beta[1:37]*0.01)^2)#diag(37)
D=diag(rep(1,37))
invsqrtSigma=solve(sqrt(Sigma))
dtilde=invsqrtSigma%*%simdata
Dtilde=invsqrtSigma%*%D
DTD=t(Dtilde)%*%Dtilde
Nulm=matrix(rep(0,676),26,26)
A=matrix(c(rbind(DTD,C),rbind(t(C),Nulm)),63,63)
Ainv=solve(A)
Beta=Ainv%*%c(t(Dtilde)%*%dtilde,rep(0,26))
Beta[1:37]
P=Dtilde%*%Ainv[1:37,1:37]%*%t(Dtilde)
plot((dtilde-P%*%dtilde)[19:37])#Residualplot i forhold til indeks
plot((sqrt(Sigma)%*%(P%*%dtilde))[1:18],(dtilde-P%*%dtilde)[1:18])#residualplot ift fittede værdier
# P%*%dtilde
P%*%simdata
(dtilde-P%*%dtilde)#residualer på den ene måde
dtilde-Dtilde%*%Beta[1:37]#residualer på den anden måde
# Beta[1:18]
# (sqrt(Sigma)%*%(P%*%dtilde))[1:18]
lev=c()
for(i in 1:37)
{
  lev=c(lev,P[i,i])
}
d=simdata
for(i in 1:37)
{
  if(min(lev)<0.9999999999)
  {
    remove=match(min(lev),lev)
    Sigma=Sigma[-remove,-remove]
    D=D[-remove,]
    d=d[-remove]
    invsqrtSigma=solve(sqrt(Sigma))
    Dtilde=invsqrtSigma%*%D
    DTD=t(Dtilde)%*%Dtilde
    A=matrix(c(rbind(DTD,C),rbind(t(C),Nulm)),63,63)
    Ainv=solve(A)
    P=Dtilde%*%Ainv[1:37,1:37]%*%t(Dtilde)
    lev=c()
    for(j in 1:length(d))
    {
      lev=c(lev,P[j,j])
    }
    print(lev)
  }
}
D #D efter vi har fjernet målinger med leverage metoden
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