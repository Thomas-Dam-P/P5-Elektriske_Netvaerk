library("pracma")
Sigma=read.table("Sigmany.txt")
Sigma=as.matrix(Sigma)
#Sigma=Sigma/400^2
#Sigma[11,11]=Sigma[11,11]*400^2
C=read.table("Cny.txt")
C=as.matrix(C)
D=read.table("D.txt")
D=as.matrix(D)
effekt=read.table("effekt.txt")
effekt=as.matrix(effekt)
d=c(effekt/400,400)
invsqrtSigma=solve(sqrt(Sigma))
dtilde=invsqrtSigma%*%d
Dtilde=invsqrtSigma%*%D
DTD=t(Dtilde)%*%Dtilde
Nulm=matrix(rep(0,676),26,26)
A=matrix(c(rbind(DTD,C),rbind(t(C),Nulm)),63,63)
Ainv=solve(A)
Beta=Ainv%*%c(t(Dtilde)%*%dtilde,rep(0,26))
Beta[1:37]
for(i in 1:50)
{
  d=c(effekt/c(Beta[21:22],Beta[25],Beta[27],Beta[29],Beta[31],Beta[33],Beta[35:37]),400)
  dtilde=invsqrtSigma%*%d
  Beta=Ainv%*%c(t(Dtilde)%*%dtilde,rep(0,26))
}
Beta[1:37]
konf=matrix(rep(0,74),37,2)
for(i in 1:37)
{
  konf[i,1]=Beta[i]-1.959963984540*sqrt(Ainv[i,i])
  konf[i,2]=Beta[i]+1.959963984540*sqrt(Ainv[i,i])
}
konf
z=c()
T=360
for(i in 19:37)
{
  z=c(z,(Beta[i]-T)/sqrt(Ainv[i,i])) #Wald's test
}
z
pnorm(z)# De giver alle et, da det er meget usandsynligt, at den rigtige værdi ligger under 360V.