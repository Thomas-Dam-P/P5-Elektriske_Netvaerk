#Determinationskoeficient eksempel med udgangspunkt i datasættet cars:
plot.default(cars,ylab ="y",xlab="x")
#data:
y<-cars$dist
x1<-cars$speed
data<-matrix(c(y,x1),50,2)
colnames(data)<-c("y-værdier","x-værdier")
print(data)
#første søjle i designmatrix: 
k<-c(rep.int(1,50))
#Designmatrix og transponeret designmatrix:
X<-matrix(c(k,x1),50,2)
tX<-t(X)
print(X)
print(tX)
#Hatmatrix og hatmatrice for null modellen:
P1<-X%*%solve(tX%*%X)%*%tX
P0<-k%*%solve(t(k)%*%k)%*%t(k)
#devians beregnes:
devians0<-t(c(y-P0%*%y))%*%c(y-P0%*%y) 
devians1<-t(c(P1%*%y-P0%*%y))%*%c(P1%*%y-P0%*%y)
devians2<-t(c(y-P1%*%y))%*%c(y-P1%*%y)
#beregning af determinationskoefficient: 
R2<-round(devians1/devians0,3)
print(paste("Determinationskoefficent er da",R2))
#beregning af den justerede determinationskoefficient: 
a<-devians2/(50-2)
b<-devians0/(50-1)
c<-a/b
R2adj1<-round(1-c,3)
print(paste("Den justerede determinationskoefficent er da",R2adj1)) 