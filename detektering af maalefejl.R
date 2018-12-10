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
mu=P%*%ytilde
V=c()
for (i in 1:length(ytilde)){
  V=c(V,(ytilde[i]-mu[i])^2)
}
sum(V)
V2=c()
for (i in 1:length(ytilde)){
  V2=c(V2,(ytilde[i]-(sum(ytilde)/length(ytilde)))^2)
}
R2<-1-(sum(V)/sum(V2))
R2
#Vi plotter residualer:
res0=ytilde-P%*%ytilde
plot((ytilde-P%*%ytilde)[1:18])#strømstyrker
plot((ytilde-P%*%ytilde)[19:37])#spændinger
#Indfører målefejl, og tester om den kan detekteres med standardiserede residualer:
y[29]=440
ytilde=invsqrtSigma%*%y
P=Xtilde%*%Ainv[1:37,1:37]%*%t(Xtilde)
v=c()
for(i in 1:length(y))
{
  v=c(v,(ytilde[i]-(P%*%ytilde)[i])/sqrt(1-P[i,i]))
}
v
#For spændinger kan vi sige nøjagtig hvilken knude, der er målefejl ved. For strømstyrker kan vi kun sige hvilken del af træet fejlen er i.


# #Da der er store residualer i 1,2,3 indgang prøver vi at fjerne dem én efter én, og tjekker, om det giver bedre residualer:
# #Fjerner den førstemåling:
# Xny=X[-1,]
# Sigmany=Sigma[-1,-1]
# yny=y[-1]
# invsqrtSigmany=solve(sqrt(Sigmany))
# ytildeny=invsqrtSigmany%*%yny
# Xtildeny=invsqrtSigmany%*%Xny
# XTXny=t(Xtildeny)%*%Xtildeny
# Any=matrix(c(rbind(XTXny,R),rbind(t(R),Nulm)),63,63)
# Ainvny=solve(Any)
# Betany=Ainvny%*%c(t(Xtildeny)%*%ytildeny,rep(0,26))
# Pny=Xtildeny%*%Ainvny[1:37,1:37]%*%t(Xtildeny)
# v=c()
# for(i in 1:length(yny))
# {
#   v=c(v,(ytildeny[i]-(Pny%*%ytildeny)[i])/sqrt(1-Pny[i,i]))
# }
# v
# #Stadig store residualer, så vi fjerner den anden måling i stedet:
# Xny=X[-2,]
# Sigmany=Sigma[-2,-2]
# yny=y[-2]
# invsqrtSigmany=solve(sqrt(Sigmany))
# ytildeny=invsqrtSigmany%*%yny
# Xtildeny=invsqrtSigmany%*%Xny
# XTXny=t(Xtildeny)%*%Xtildeny
# Any=matrix(c(rbind(XTXny,R),rbind(t(R),Nulm)),63,63)
# Ainvny=solve(Any)
# Betany=Ainvny%*%c(t(Xtildeny)%*%ytildeny,rep(0,26))
# Pny=Xtildeny%*%Ainvny[1:37,1:37]%*%t(Xtildeny)
# v=c()
# for(i in 1:length(yny))
# {
#   v=c(v,(ytildeny[i]-(Pny%*%ytildeny)[i])/sqrt(1-Pny[i,i]))
# }
# v