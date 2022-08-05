library(tidyverse)
library(scatterplot3d)
rm(list=ls())

#Simulate uniformly from unit sphere

n<-10000

x<-matrix(rnorm(n*3),n,3)
x<-t(apply(x,1,function(a){a/sqrt(sum(a^2))}))

theta<-acos(x[,3])

phi<-atan(x[,2]/x[,1])

phi[(x[,2]<0)&(x[,1]<0)]<-phi[(x[,2]<0)&(x[,1]<0)]+pi
phi[(x[,2]>0)&(x[,1]<0)]<-phi[(x[,2]>0)&(x[,1]<0)]-pi

y<-cbind(theta*cos(phi),theta*sin(phi))
ycap<-y[(theta<1),]

#Checks
cc<-rep('blue',n)
cc[theta<1]<-'red'
plot(y[,1],y[,2],col=cc,asp=1,pch=20)
plot(ycap[,1],ycap[,2],col=cc[theta<1],asp=1,pch=20)
scatterplot3d(x[,1],x[,2],x[,3],color=cc,asp=1,pch=20)


caparea<-2*pi*(1-cos(1/pi))
areaball_local<-pi

meandetH1<-caparea/areaball_local

#Next Steps

#1) Treat y as an output embedding and x as an input. Run the LearnMetric algorithm to get a H at each y point
#2) Find Average root determinant of H in the spherical cap
