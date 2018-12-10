C=read.table("Cny.txt")
C=as.matrix(C)
Nulm=matrix(rep(0,676),26,26)
#Function til at estimerer beta:
betaestimat<-function(d,D,Sigma)
{invsqrtSigma=solve(sqrt(Sigma))
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
Beta[1:37]}
y=c((effekt)/c(Beta[21:22],Beta[25],Beta[27],Beta[29],Beta[31],Beta[33],Beta[35:37]),400)
D=read.table("D.txt")
D=as.matrix(D)
Sigma=diag((d*0.01)^2)
betaestimat(y,D,Sigma)
#Funktion til bestemmelse af konfidensintervaller: 
konfidens<-function(Beta){
  konfg=matrix(rep(0,74),37,2)
  for(i in 1:37)
  {
    konfg[i,1]=Beta[i]-1.959963984540*sqrt(Ainv[i,i])
    konfg[i,2]=Beta[i]+1.959963984540*sqrt(Ainv[i,i])
  }
  konfg
}
konfidens(Beta)
#Funktion til levarage: 
levarage<-function(D,Sigma){
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
  }}  
}

