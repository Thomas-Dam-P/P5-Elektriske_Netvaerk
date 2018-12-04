#data model 2: 
y<-cars$dist
x1<-cars$speed
x2<-x1^2
data<-matrix(c(y,x1,x2),50,3)
colnames(data)<-c("y-værdier","x-værdier","x^2-værdier")
print(data)
#første søjle i designmatrix: 
k<-c(rep.int(1,50))
#designmatrix og transponeret designmatrix: 
X<-matrix(c(k,x1,x2),50,3) 
tX<-t(X)
print(X)
print(tX)
#Hatmatrix og hatmatrice for null modellen:
P1<-X%*%solve(tX%*%X)%*%tX
P0<-k%*%solve(t(k)%*%k)%*%t(k)
#devians beregnes: 
devians0<-t(c(P1%*%y-P0%*%y))%*%c(P1%*%y-P0%*%y)
devians1<-t(c(y-P1%*%y))%*%c(y-P1%*%y)
devians2<-t(c(y-P0%*%y))%*%c(y-P0%*%y) 
#den justerede determinationskoefficient for model 1 og model 2 beregnes: 
a<-devians1/(50-3)
b<-devians2/(50-1)
c<-a/b
R2adj2<-round(1-c,3)
print(paste("Den justerede determinationskoefficient er da",R2adj2))
ff
