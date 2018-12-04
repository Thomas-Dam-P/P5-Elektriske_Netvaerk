#determinationskoeficent eksempel med udgangspunkt i datasættet cars:
plot(cars)
k<-c(rep.int(1,50))
x1<-cars$speed
#Designmatrix:
X<-matrix(c(k,x1),50,2)
tX<-t(X)
#Hatmatrix og hatmatrice for null modellen:
P1<-X%*%solve(tX%*%X)%*%tX
P0<-x1%*%solve(t(x1)%*%x1)%*%t(x1)
#devarians beregnes: 
devarians1<-t(c(y-P1%*%y))%*%c(y-P1%*%y)
devarians0<-t(c(y-P0%*%y))%*%c(y-P0%*%y)
#beregning af determinationskoefficent: 
R2<-1-(devarians1/devarians0)
#giver ikke den rigtige R2..
#y-værdier: 
y<-cars$dist
#estimation af beta: 
beta<-solve(tX%*%X)%*%tX%*%y
#test for rigtigt beta ved benyttese af lm funktionen: 
betatest<-lm(dist~speed,cars)
summary(betatest)
#R2=0,65
#R2adj=0,64
#ny model 
model1R2adj<-1-((devarians1/50-2)/(devarians0/50-1))
summary(betatest)
betatest2<-lm(dist~speed+I(speed^2),cars)
summary(betatest2)
#R2adj =0,65
