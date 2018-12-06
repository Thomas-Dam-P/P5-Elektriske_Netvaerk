#Determinationskoefficienten beregnes for model 1: 
effekt=read.table("effekt.txt")
effekt=as.matrix(effekt)
d=c(effekt/400,400)
Sigma=diag((d*0.01)^2)
invsqrtSigma=solve(sqrt(Sigma))
dtilde=invsqrtSigma%*%d
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
#Den justerede determinationskoefficient beregnes for model 1 ud fra følgende relation mellem R2 og R2adj
#R2adj=1-((n-1)/(n-k))(1-R2):
R2adj<-1-((11-1)/(11-37))*(1-R2)
R2adj
