library("pracma")
#Sigma=read.table("Sigmany400.txt")
#Sigma=as.matrix(Sigma)
#Sigma=Sigma/400^2
#Sigma[11,11]=Sigma[11,11]*400^2
C=read.table("Cny.txt")
C=as.matrix(C)
D=read.table("D.txt")
D=as.matrix(D)
Dg=D
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
for(i in 1:50)
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
konfg=matrix(rep(0,74),37,2)
for(i in 1:37)
{
  konfg[i,1]=Beta[i]-1.959963984540*sqrt(Ainv[i,i])
  konfg[i,2]=Beta[i]+1.959963984540*sqrt(Ainv[i,i])
}
konfg#Konfidensintervaller for parametrene:


z=c()#H_0i: spændingen i knude i er mindre end T volt:
T=360
for(i in 19:37)
{
  z=c(z,(Beta[i]-T)/sqrt(Ainv[i,i])) #Wald's test
}
z
1-pnorm(z)

#T-test for 440 W: 
z=c()
T=440
for(i in 19:37)
{
  z=c(z,(Beta[i]-T)/sqrt(Ainv[i,i])) #Wald's test
}
z
pnorm(z)
# p-værdier udregnes. hvis p-værdien er under 0.05, så kan vi med 95% sikkerhed sige, at den knude ikke har spænding lavere end T (forkaste H_0i).
#Nedenfor piller vi målinger med mindst leverage fra, indtil man ikke kan undvære flere
Sigma=diag((Beta[1:37]*0.01)^2)#diag(37)
D=diag(rep(1,37))
invsqrtSigma=solve(sqrt(Sigma))
Dtilde=invsqrtSigma%*%D
DTD=t(Dtilde)%*%Dtilde
Nulm=matrix(rep(0,676),26,26)
A=matrix(c(rbind(DTD,C),rbind(t(C),Nulm)),63,63)
Ainv=solve(A)
P=Dtilde%*%Ainv[1:37,1:37]%*%t(Dtilde)
lev=c()
for(i in 1:37)
{
  lev=c(lev,P[i,i])
}
v=rep(1,37)
for(i in 1:37)
{
  if(min(lev)<0.9999999999)
  {
    remove=match(min(lev),lev)
    Sigma=Sigma[-remove,-remove]
    D=D[-remove,]
    v=v[-remove]
    invsqrtSigma=solve(sqrt(Sigma))
    Dtilde=invsqrtSigma%*%D
    DTD=t(Dtilde)%*%Dtilde
    A=matrix(c(rbind(DTD,C),rbind(t(C),Nulm)),63,63)
    Ainv=solve(A)
    P=Dtilde%*%Ainv[1:37,1:37]%*%t(Dtilde)
    lev=c()
    for(j in 1:length(v))
    {
      lev=c(lev,P[j,j])
    }
    print(lev)
  }
}
#Determinationskoefficient for modellen:
mu<-P%*%dtilde
 
  V<-c()
for (i in 1:11){
  V<-c(V,(dtilde[i]-mu[i])^2)
}
sum(V)
V2<-c()
for (i in 1:11){
  V2<-c(V2,(dtilde[i]-(sum(dtilde)/11))^2)
}
R2<-1-(sum(V)/sum(V2))
R2
#Her estimeres beta ud fra den nye D og Sigma:
  d=c(Beta[2:3],Beta[6],Beta[8],Beta[10:12],Beta[14],Beta[16],Beta[18],Beta[25])
  Sigma=diag((d*0.01)^2)
  invsqrtSigma=solve(sqrt(Sigma))
  dtilde=invsqrtSigma%*%d
  Dtilde=invsqrtSigma%*%D
  DTD=t(Dtilde)%*%Dtilde
  A=matrix(c(rbind(DTD,C),rbind(t(C),Nulm)),63,63)
  Ainv=solve(A)
  Beta=Ainv%*%c(t(Dtilde)%*%dtilde,rep(0,26))
#Dette giver samme estimat som før.
#Nu findes konfidens intervaller med den nye D og Sigma:
konf=matrix(rep(0,74),37,2)
for(i in 1:37)
{
  konf[i,1]=Beta[i]-1.959963984540*sqrt(Ainv[i,i])
  konf[i,2]=Beta[i]+1.959963984540*sqrt(Ainv[i,i])
}
konf#Konfidensintervaller for parametrene:
konfg
konfg-konf#forskel i konfidensintervallet fra før og det nye.

#leverage kun afhængigt af netværk, med anden varians: 
varians<-c(rep((0.01*0.0155)^2,18),rep((0.01*1)^2,19))
varians
Sigma=diag(varians)
invsqrtSigma=solve(sqrt(Sigma))
D=diag(rep(1,37))
Dtilde=invsqrtSigma%*%D
DTD=t(Dtilde)%*%Dtilde
Nulm=matrix(rep(0,676),26,26)
A=matrix(c(rbind(DTD,C),rbind(t(C),Nulm)),63,63)
Ainv=solve(A)
P=Dtilde%*%Ainv[1:37,1:37]%*%t(Dtilde)
lev=c()
for(i in 1:37)
{
  lev=c(lev,P[i,i])
}
v=rep(1,37)
for(i in 1:37)
{
  if(min(lev)<0.9999999999)
  {
    remove=match(min(lev),lev)
    Sigma=Sigma[-remove,-remove]
    D=D[-remove,]
    v=v[-remove]
    invsqrtSigma=solve(sqrt(Sigma))
    Dtilde=invsqrtSigma%*%D
    DTD=t(Dtilde)%*%Dtilde
    A=matrix(c(rbind(DTD,C),rbind(t(C),Nulm)),63,63)
    Ainv=solve(A)
    P=Dtilde%*%Ainv[1:37,1:37]%*%t(Dtilde)
    lev=c()
    for(j in 1:length(v))
    {
      lev=c(lev,P[j,j])
    }
    print(lev)
  }
}
D
d=c(Beta[2:3],Beta[5],Beta[8],Beta[10],Beta[12],Beta[14],Beta[16:18],Beta[36])
Sigma=diag((d*0.01)^2)
invsqrtSigma=solve(sqrt(Sigma))
dtilde=invsqrtSigma%*%d
Dtilde=invsqrtSigma%*%D
DTD=t(Dtilde)%*%Dtilde
A=matrix(c(rbind(DTD,C),rbind(t(C),Nulm)),63,63)
Ainv=solve(A)
Beta=Ainv%*%c(t(Dtilde)%*%dtilde,rep(0,26))
#Dette giver samme estimat som før.
#Nu findes konfidens intervaller med den nye D og Sigma:
konf=matrix(rep(0,74),37,2)
for(i in 1:37)
{
  konf[i,1]=Beta[i]-1.959963984540*sqrt(Ainv[i,i])
  konf[i,2]=Beta[i]+1.959963984540*sqrt(Ainv[i,i])
}
konfg-konf
konfg
#D efter vi har fjernet målinger med leverage metoden
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