###############################################
Distance decay
###############################################
rm(list=ls())
data<- read.csv2("C:/Fungal ecology/Complete_data_distance_decay.csv", header=T)
head(data)
attach(data)

data_1<-data[,-c(1:8)]
head(data_1)
attach(data_1)

distance_dist.0<-read.table("C:/Fungal ecology/Distance_km.txt", header=T)
head(distance_dist.0)
attach(distance_dist.0)
distance_dist.1<-vegdist(distance_dist.0,method="euclidean")
distance_dist.2<-log(distance_dist.1)
max(distance_dist.2)
distance_dist.3<-7.751263-(distance_dist.2)

divers<-vegdist(wisconsin(sqrt(data_1)), method="jaccard")
plot(distance_dist.3,divers,pch=".", cex=3,xlim=c(0,13),ylim=c(0,1.2))
abline(lm(divers~distance_dist.2), col="black", lwd=2)
summary(lm(divers~distance_dist.2))


###############################################
Distance decay_without singletons
###############################################
rm(list=ls())
data<- read.table("C:/Fungal ecology/2 agosto/tudo.txt", header=T)
head(data)
attach(data)

data_1<-data[,-c(1:10)]
head(data_1)
attach(data_1)

distance_dist.0<-read.table("C:/Fungal ecology//2 agosto/Distancia_km.txt", header=T)
head(distance_dist.0)
attach(distance_dist.0)
distance_dist.1<-vegdist(distance_dist.0,method="euclidean")
distance_dist.2<-log(distance_dist.1)
max(distance_dist.2)
distance_dist.3<-(distance_dist.2)

divers<-vegdist(wisconsin(sqrt(data_1)), method="jaccard")
plot(distance_dist.3,divers,pch=".", cex=3,xlim=c(0,13),ylim=c(0,1.2))
abline(lm(divers~distance_dist.2), col="black", lwd=2)
summary(lm(divers~distance_dist.2))



