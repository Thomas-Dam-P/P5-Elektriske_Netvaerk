Sigma=diag((y*0.01)^2)
X=diag(1,37)
invsqrtSigma=solve(sqrt(Sigma))
ytilde=invsqrtSigma%*%y
Xtilde=invsqrtSigma%*%X
XTX=t(Xtilde)%*%Xtilde
A=matrix(c(rbind(XTX,R),rbind(t(R),Nulm)),63,63)
Ainv=solve(A)
Beta3=Ainv%*%c(t(Xtilde)%*%ytilde,rep(0,26))
konf=matrix(rep(0,74),37,2)
for(i in 1:37)
{
  konf[i,1]=Beta3[i]-1.959963984540*sqrt(Ainv[i,i])
  konf[i,2]=Beta3[i]+1.959963984540*sqrt(Ainv[i,i])
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
#Vi fjerner 5 målinger ved leverage metoden
varians=c(rep((0.01*0.0155)^2,18),rep((0.01*1)^2,19))
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
y10=y
for(i in 1:10)
{
  if(min(lev)<0.9999999999)
  {
    remove=match(min(lev),lev)
    Sigma=Sigma[-remove,-remove]
    Sigma1=Sigma1[-remove,-remove]
    y10=y10[-remove]
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
ytilde=invsqrtSigma%*%y10
Xtilde=invsqrtSigma%*%X
XTX=t(Xtilde)%*%Xtilde
A=matrix(c(rbind(XTX,R),rbind(t(R),Nulm)),63,63)
Ainv=solve(A)
Beta10=Ainv%*%c(t(Xtilde)%*%ytilde,rep(0,26))
Beta10[1:37]#Estimat med 5 fjernede målinger.
konf10=matrix(rep(0,74),37,2)
#konfidensinterval med 5 fjernede målinger:
for(i in 1:37)
{
  konf10[i,1]=Beta10[i]-1.959963984540*sqrt(Ainv[i,i])
  konf10[i,2]=Beta10[i]+1.959963984540*sqrt(Ainv[i,i])
}
konf10
konf
konf-konf10
#Vi finder R^2 med 5 fjernede målinger:
P=Xtilde%*%Ainv[1:37,1:37]%*%t(Xtilde)
mu=P%*%ytilde
V=c()
for (i in 1:length(ytilde)){
  V=c(V,(ytilde[i]-mu[i])^2)
}
V2=c()
for (i in 1:length(ytilde)){
  V2=c(V2,(ytilde[i]-(sum(ytilde)/length(ytilde)))^2)
}
R210=1-(sum(V)/sum(V2))
R210adj=1-(((37-1)/(37-length(ytilde))))*(1-R210)