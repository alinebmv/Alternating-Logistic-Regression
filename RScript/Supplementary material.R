############################################################################################
############################################################################################
Manuscript informations
############################################################################################
############################################################################################
Article title: Fungal endophytes associated with three South American Myrtae (Myrtaceae) exhibit preferences in the colonization at leaf level
Journal name: Fungal Biology
Author names: Aline B. M. VAZ, Andre G.F.C DA COSTA, Lucelia V.V. RAAD, Aristoteles GoES-NETO

############################################################################################
############################################################################################
Analysis
############################################################################################
############################################################################################

rm(list=ls())
data<- read.csv2(file = "data.csv", header=T)
attach(data)
require(geepack)

############################################################################################
############################################################################################
1. Fungi taxonomic grouping
############################################################################################
############################################################################################
dim(data)
sordariomycetes<-D1+D2+D3+D4+D5+D6+Dp+Ds+Dh+Gree+Am1+C1+C2+C3+C4+C5+C6+C7+C8+Cep+X1+X2+X3+X4+An1+An2+An3
sordariomycetes<- ifelse(sordariomycetes>=1,1,0)

dothideomycetes<-Gm+Dcan+Pb+Myc+Cs+Cc+Pa+Li+Li
dothideomycetes<- ifelse(dothideomycetes>=1,1,0)

xylariales<-X1+X2+X3+X4+An1+An2+An3
xylariales<- ifelse(xylariales>=1,1,0)

capnodiales<-Pb+Myc+Cs+Cc
capnodiales<- ifelse(capnodiales>=1,1,0)

xylaria<-X1+X2+X3+X4
xylaria<-ifelse(xylaria>=1,1,0)

data$sordariomycetes<- sordariomycetes
data$xylariales<- xylariales
data$xylaria<- xylaria
data$dothideomycetes<- dothideomycetes
data$capnodiales<- capnodiales

############################################################################################
############################################################################################
2.Preparing the dependence structure (matrix design)_Sordariomycetes
############################################################################################
############################################################################################
n = as.vector(table(as.factor(data$Pl)))
last <- cumsum(n);
first <- last - n + 1;
Design <- NULL
for ( i in 1:length(n) )
{
n.i <- n[i]
id.i <- data$Pl[first[i]]
ind.i <- data$Ind[ first[i]:last[i] ]
lf.i <- data$lf[ first[i]:last[i] ]
fr.i <- data$Fr[ first[i]:last[i] ]
dis.i <- data$dis[ first[i]:last[i] ]
resp.i<- data$sordariomycetes[ first[i]:last[i] ]

l <- 1
if (n.i == 1) {z.i <- cbind(NA,NA,NA, NA, NA, NA)}
else
{
## Note: ch2(m) = m(m - 1)/2
id.i <- rep(id.i, choose(n.i, 2))
z.i1 <- rep(NA, choose(n.i, 2) )
z.i2 <- rep(NA, choose(n.i, 2) )
z.i3 <- rep(NA, choose(n.i, 2) )
z.i4 <- rep(NA, choose(n.i, 2) )
z.i5 <- rep(NA, choose(n.i, 2) )

for( j in seq(1, n.i - 1) )
{
for( k in seq(j+1, n.i) )
{
z.i1[l] <- ifelse(ind.i[j]- ind.i[k]==0,1,0) 
z.i2[l] <- ifelse(lf.i[j]- lf.i[k]==0,1,0)
z.i3[l] <- ifelse(fr.i[j]== fr.i[k],1,0)
z.i4[l] <- abs(dis.i[j] -  dis.i[k])
z.i5[l] <- ifelse(resp.i[j]&resp.i[k]==1,1,0)

l <- l+1
}
}
z.i <- cbind(id.i, z.i1, z.i2, z.i3, z.i4, z.i5)
}
Design <- rbind(Design, z.i)
}

Design<- data.frame(Design)
DesignF<- Design
DesignF_sordariomycetes<-DesignF
colnames(DesignF_sordariomycetes)<- c("Collection_site", "Host_tree_individual", "Leaf", "Leaf_fragment", "Distance", "sordariomycetes")
z0a_sordariomycetes<- cbind(ifelse(DesignF[,1]==1000,0,1), DesignF[,2], DesignF[,3], DesignF[,5])
head(z0a_sordariomycetes)

############################################################################################
############################################################################################
3.Preparing the dependence structure (matrix design)_Xylariales
############################################################################################
############################################################################################
n = as.vector(table(as.factor(data$Pl)))
last <- cumsum(n);
first <- last - n + 1;
Design <- NULL
for ( i in 1:length(n) )
{
n.i <- n[i]
id.i <- data$Pl[first[i]]
ind.i <- data$Ind[ first[i]:last[i] ]
lf.i <- data$lf[ first[i]:last[i] ]
fr.i <- data$Fr[ first[i]:last[i] ]
dis.i <- data$dis[ first[i]:last[i] ]
resp.i<- data$xylariales[ first[i]:last[i] ]

l <- 1
if (n.i == 1) {z.i <- cbind(NA,NA,NA, NA, NA, NA)}
else
{
## Note: ch2(m) = m(m - 1)/2
id.i <- rep(id.i, choose(n.i, 2))
z.i1 <- rep(NA, choose(n.i, 2) )
z.i2 <- rep(NA, choose(n.i, 2) )
z.i3 <- rep(NA, choose(n.i, 2) )
z.i4 <- rep(NA, choose(n.i, 2) )
z.i5 <- rep(NA, choose(n.i, 2) )

for( j in seq(1, n.i - 1) )
{
for( k in seq(j+1, n.i) )
{
z.i1[l] <- ifelse(ind.i[j]- ind.i[k]==0,1,0) 
z.i2[l] <- ifelse(lf.i[j]- lf.i[k]==0,1,0)
z.i3[l] <- ifelse(fr.i[j]== fr.i[k],1,0)
z.i4[l] <- abs(dis.i[j] -  dis.i[k])
z.i5[l] <- ifelse(resp.i[j]&resp.i[k]==1,1,0)

l <- l+1
}
}
z.i <- cbind(id.i, z.i1, z.i2, z.i3, z.i4, z.i5)
}
Design <- rbind(Design, z.i)
}
Design<- data.frame(Design)
DesignF<- Design

DesignF_xylariales<-DesignF
colnames(DesignF_xylariales)<- c("Collection_site", "Host_tree_individual", "Leaf", "Leaf_fragment", "Distance", "xylariales")
z0a_xylariales<- cbind(ifelse(DesignF[,1]==1000,0,1), DesignF[,2], DesignF[,3], DesignF[,5])
head(z0a_xylariales)

############################################################################################
############################################################################################
4.Preparing the dependence structure (matrix design)_Xylaria
############################################################################################
############################################################################################
n = as.vector(table(as.factor(data$Pl)))
last <- cumsum(n);
first <- last - n + 1;
Design <- NULL
for ( i in 1:length(n) )
{
n.i <- n[i]
id.i <- data$Pl[first[i]]
ind.i <- data$Ind[ first[i]:last[i] ]
lf.i <- data$lf[ first[i]:last[i] ]
fr.i <- data$Fr[ first[i]:last[i] ]
dis.i <- data$dis[ first[i]:last[i] ]
resp.i<- data$xylaria[ first[i]:last[i] ]

l <- 1
if (n.i == 1) {z.i <- cbind(NA,NA,NA, NA, NA, NA)}
else
{
## Note: ch2(m) = m(m - 1)/2
id.i <- rep(id.i, choose(n.i, 2))
z.i1 <- rep(NA, choose(n.i, 2) )
z.i2 <- rep(NA, choose(n.i, 2) )
z.i3 <- rep(NA, choose(n.i, 2) )
z.i4 <- rep(NA, choose(n.i, 2) )
z.i5 <- rep(NA, choose(n.i, 2) )

for( j in seq(1, n.i - 1) )
{
for( k in seq(j+1, n.i) )
{
z.i1[l] <- ifelse(ind.i[j]- ind.i[k]==0,1,0) 
z.i2[l] <- ifelse(lf.i[j]- lf.i[k]==0,1,0)
z.i3[l] <- ifelse(fr.i[j]== fr.i[k],1,0)
z.i4[l] <- abs(dis.i[j] -  dis.i[k])
z.i5[l] <- ifelse(resp.i[j]&resp.i[k]==1,1,0)

l <- l+1
}
}
z.i <- cbind(id.i, z.i1, z.i2, z.i3, z.i4, z.i5)
}
Design <- rbind(Design, z.i)
}

Design<- data.frame(Design)
DesignF<- Design

DesignF_xylaria<-DesignF
colnames(DesignF_xylaria)<- c("Collection_site", "Host_tree_individual", "Leaf", "Leaf_fragment", "Distance", "xylaria")
z0a_xylaria<- cbind(ifelse(DesignF[,1]==1000,0,1), DesignF[,2], DesignF[,3], DesignF[,5])
head(z0a_xylaria)

############################################################################################
############################################################################################
5.Preparing the dependence structure (matrix design)_dothideomycetes
############################################################################################
############################################################################################
n = as.vector(table(as.factor(data$Pl)))
last <- cumsum(n);
first <- last - n + 1;
Design <- NULL
for ( i in 1:length(n) )
{
n.i <- n[i]
id.i <- data$Pl[first[i]]
ind.i <- data$Ind[ first[i]:last[i] ]
lf.i <- data$lf[ first[i]:last[i] ]
fr.i <- data$Fr[ first[i]:last[i] ]
dis.i <- data$dis[ first[i]:last[i] ]
resp.i<- data$dothideomycetes[ first[i]:last[i] ]

l <- 1
if (n.i == 1) {z.i <- cbind(NA,NA,NA, NA, NA, NA)}
else
{
## Note: ch2(m) = m(m - 1)/2
id.i <- rep(id.i, choose(n.i, 2))
z.i1 <- rep(NA, choose(n.i, 2) )
z.i2 <- rep(NA, choose(n.i, 2) )
z.i3 <- rep(NA, choose(n.i, 2) )
z.i4 <- rep(NA, choose(n.i, 2) )
z.i5 <- rep(NA, choose(n.i, 2) )

for( j in seq(1, n.i - 1) )
{
for( k in seq(j+1, n.i) )
{
z.i1[l] <- ifelse(ind.i[j]- ind.i[k]==0,1,0) 
z.i2[l] <- ifelse(lf.i[j]- lf.i[k]==0,1,0)
z.i3[l] <- ifelse(fr.i[j]== fr.i[k],1,0)
z.i4[l] <- abs(dis.i[j] -  dis.i[k])
z.i5[l] <- ifelse(resp.i[j]&resp.i[k]==1,1,0)

l <- l+1
}
}
z.i <- cbind(id.i, z.i1, z.i2, z.i3, z.i4, z.i5)
}
Design <- rbind(Design, z.i)
}

Design<- data.frame(Design)
DesignF<- Design
DesignF_dothideomycetes<-DesignF
colnames(DesignF_dothideomycetes)<- c("Collection_site", "Host_tree_individual", "Leaf", "Leaf_fragment", "Distance", "dothideomycetes")
z0a_dothideomycetes<- cbind(ifelse(DesignF[,1]==1000,0,1), DesignF[,2], DesignF[,3], DesignF[,5])
head(z0a_dothideomycetes)

############################################################################################
############################################################################################
6.Preparing the dependence structure (matrix design)_Capnodiales
############################################################################################
############################################################################################
n = as.vector(table(as.factor(data$Pl)))
last <- cumsum(n);
first <- last - n + 1;
Design <- NULL
for ( i in 1:length(n) )
{
n.i <- n[i]
id.i <- data$Pl[first[i]]
ind.i <- data$Ind[ first[i]:last[i] ]
lf.i <- data$lf[ first[i]:last[i] ]
fr.i <- data$Fr[ first[i]:last[i] ]
dis.i <- data$dis[ first[i]:last[i] ]
resp.i<- data$capnodiales[ first[i]:last[i] ]

l <- 1
if (n.i == 1) {z.i <- cbind(NA,NA,NA, NA, NA, NA)}
else
{
## Note: ch2(m) = m(m - 1)/2
id.i <- rep(id.i, choose(n.i, 2))
z.i1 <- rep(NA, choose(n.i, 2) )
z.i2 <- rep(NA, choose(n.i, 2) )
z.i3 <- rep(NA, choose(n.i, 2) )
z.i4 <- rep(NA, choose(n.i, 2) )
z.i5 <- rep(NA, choose(n.i, 2) )

for( j in seq(1, n.i - 1) )
{
for( k in seq(j+1, n.i) )
{
z.i1[l] <- ifelse(ind.i[j]- ind.i[k]==0,1,0) 
z.i2[l] <- ifelse(lf.i[j]- lf.i[k]==0,1,0)
z.i3[l] <- ifelse(fr.i[j]== fr.i[k],1,0)
z.i4[l] <- abs(dis.i[j] -  dis.i[k])
z.i5[l] <- ifelse(resp.i[j]&resp.i[k]==1,1,0)

l <- l+1
}
}
z.i <- cbind(id.i, z.i1, z.i2, z.i3, z.i4, z.i5)
}
Design <- rbind(Design, z.i)
}

Design<- data.frame(Design)
DesignF<- Design
DesignF_capnodiales<-DesignF
colnames(DesignF_capnodiales)<- c("Collection_site", "Host_tree_individual", "Leaf", "Leaf_fragment", "Distance", "capnodiales")
z0a_capnodiales<- cbind(ifelse(DesignF[,1]==1000,0,1), DesignF[,2], DesignF[,3], DesignF[,5])
head(z0a_capnodiales)

############################################################################################
############################################################################################
7. Principal components analysis of environmental variables
############################################################################################
############################################################################################

fit1<- prcomp(cbind(ele, wp, temp),scale=T)
fit1
summary(fit1)
PCA1<- fit1$x[,1]

############################################################################################
############################################################################################
8. Considering the C leaf fragment on the intercept
############################################################################################
############################################################################################

Fr_T<- ifelse(Fr=="C", "0C", ifelse(Fr=="A","A",ifelse(Fr=="B","B",ifelse(Fr=="D","D",ifelse(Fr=="E","E","F")))))

############################################################################################
############################################################################################
9. Mean structure graphs
############################################################################################
############################################################################################

### Principal component analysis_Figure 2

sordariomycetes_PCA<-tapply(sordariomycetes, PCA1, mean)
xylariales_PCA<-tapply(xylariales, PCA1, mean)
xylaria_PCA<-tapply(xylaria, PCA1, mean)
dothideomycetes_PCA<-tapply(dothideomycetes, PCA1, mean)
capnodiales_PCA<-tapply(capnodiales, PCA1, mean)

tiff(filename = "figure2.tiff" , width = 600, height =600, units="px")
par(mfrow=c(1,1))
plot(names(table(PCA1)), 100*sordariomycetes_PCA, type="b", pch=2, lwd=2, lty=1, col="black", cex=0.5, ylim=c(0, 35), ylab="Porcentage (%)", xlab="PC1")
lines(names(table(PCA1)),100*xylariales_PCA, type="b", pch=1, lwd=2, lty=1 ,col= "red",cex=0.5)
lines(names(table(PCA1)),100*xylaria_PCA, type="b", pch=3, lwd=2, lty=1, col="blue", cex=0.5)
lines(names(table(PCA1)),100*dothideomycetes_PCA,  type="b", pch=2,lwd=2, lty=1, cex=0.5, col="orange", ylim=c(0, 35), ylab="Porcentage (%)", xlab="PC1")
lines(names(table(PCA1)),100*capnodiales_PCA,  type="b", pch=1,lwd=2, lty=1, col="green", cex=0.5)
legend("topleft", c("Sordariomycete","Xylariales","Xylaria","Dothideomycetes","Capnodiales"),lty=1, col=c("black","red","blue","orange","green"),lwd=2, bty="n")
dev.off()

### Leaf fragment_Figure 3

sordariomycetes_Fr<-tapply(sordariomycetes, Fr, mean)
xylariales_Fr<-tapply(xylariales, Fr, mean)
xylaria_Fr<-tapply(xylaria, Fr, mean)
dothideomycetes_Fr<-tapply(dothideomycetes, Fr, mean)
capnodiales_Fr<-tapply(capnodiales, Fr, mean)

y<-cbind(sordariomycetes_Fr,xylariales_Fr,xylaria_Fr,dothideomycetes_Fr,capnodiales_Fr)

tiff(filename = "Figure3.tiff" , width = 600, height = 600, units="px")
par(mfrow=c(1,1))
barplot(t(y*100), beside=T, ylim=c(0, 35),legend=F,col=c("black","red","blue","orange","green"),ylab="Porcentage (%)", xlab="Leaf fragment")
legend("topleft", c("Sordariomycetes","Xylariales","Xylaria","Dothideomycetes","Capnodiales"), fill=c("black","red","blue","orange","green"),bty="n", cex=0.9)
dev.off()

############################################################################################
############################################################################################
10. Dependence structure graphs_for each fungal taxonomic group
############################################################################################
############################################################################################

### Individual host and leaf_ Figure 4 and 5

sordariomycetes_Ht<-tapply(DesignF_sordariomycetes$sordariomycetes, DesignF_sordariomycetes$Host_tree_individual, mean)
xylariales_Ht<-tapply(DesignF_xylariales$xylariales, DesignF_xylariales$Host_tree_individual, mean)
xylaria_Ht<-tapply(DesignF_xylaria$xylaria, DesignF_xylaria$Host_tree_individual, mean)
dothideomycetes_Ht<-tapply(DesignF_dothideomycetes$dothideomycetes, DesignF_dothideomycetes$Host_tree_individual, mean)
capnodiales_Ht<-tapply(DesignF_capnodiales$capnodiales, DesignF_capnodiales$Host_tree_individual, mean)

sordariomycetes_L<-tapply(DesignF_sordariomycetes$sordariomycetes, DesignF_sordariomycetes$Leaf, mean)
xylariales_L<-tapply(DesignF_xylariales$xylariales, DesignF_xylariales$Leaf, mean)
xylaria_L<-tapply(DesignF_xylaria$xylaria, DesignF_xylaria$Leaf, mean)
dothideomycetes_L<-tapply(DesignF_dothideomycetes$dothideomycetes, DesignF_dothideomycetes$Leaf, mean)
capnodiales_L<-tapply(DesignF_capnodiales$capnodiales, DesignF_capnodiales$Leaf, mean)

w<-cbind(sordariomycetes_Ht,xylariales_Ht,xylaria_Ht,dothideomycetes_Ht,capnodiales_Ht)
z<-cbind(sordariomycetes_L,xylariales_L,xylaria_L,dothideomycetes_L,capnodiales_L)

tiff(filename = "Figure4.tiff" , width = 600, height = 600, units="px")
par(mfrow=c(1,2))
barplot(t(w), beside=T, ylim=c(0, 0.08), col=c("black","red","blue","orange","green"),ylab=expression(Prob(Y[j]==1,Y[k]==1)), names=c("No","Yes"), xlab="Individual host tree",legend=F)
legend("topleft", fill=c("black","red","blue","orange","green"), c("Sordariomycetes","Xylariales","Xylaria","Dothideomycetes","Capnodiales"), bty="n", cex=0.9)
barplot(t(z), beside=T, ylim=c(0, 0.08), col=c("black","red","blue","orange","green"),ylab=expression(Prob(Y[j]==1,Y[k]==1)), names=c("No","Yes"),xlab="Leaf fragment", legend=F)
legend("topleft", fill=c("black","red","blue","orange","green"), c("Sordariomycetes","Xylariales","Xylaria","Dothideomycetes","Capnodiales"), bty="n", cex=0.9)
dev.off()

### Distance_ Figure 6

sordariomycetes_Dis<-tapply(DesignF_sordariomycetes$sordariomycetes, DesignF_sordariomycetes$Distance, mean)
xylariales_Dis<-tapply(DesignF_xylariales$xylariales, DesignF_xylariales$Distance, mean)
xylaria_Dis<-tapply(DesignF_xylaria$xylaria, DesignF_xylaria$Distance, mean)
dothideomycetes_Dis<-tapply(DesignF_dothideomycetes$dothideomycetes, DesignF_dothideomycetes$Distance, mean)
capnodiales_Dis<-tapply(DesignF_capnodiales$capnodiales, DesignF_capnodiales$Distance, mean)

tiff(filename = "Figure5.tiff" , width = 600, height = 600, units="px")
par(mfrow=c(1,1))
plot(sort(unique(DesignF_sordariomycetes$Dis)), sordariomycetes_Dis, pch=1, lwd=2, lty=1, cex=0, col="black", ylim=c(0,0.08), ylab=expression(Prob(Y[j]==1,Y[k]==1)), xlab="Distance (m)")
lines(lowess(sort(unique(DesignF_sordariomycetes$Dis)), sordariomycetes_Dis, 3/3), lwd=2, lty=1,col="black")
lines(lowess(sort(unique(DesignF_xylariales$Dis)), xylariales_Dis, 3/3), lwd=2, lty=1,col="red")
lines(lowess(sort(unique(DesignF_xylaria$Dis)), xylaria_Dis, 3/3), lwd=2, lty=1,col="blue")
lines(lowess(sort(unique(DesignF_dothideomycetes$Dis)), dothideomycetes_Dis, 3/3), lwd=2, lty=1,col="orange")
lines(lowess(sort(unique(DesignF_capnodiales$Dis)), capnodiales_Dis, 3/3), lwd=2, lty=1, col="green")
legend("topleft", col=c("black","red","blue","orange","green"), lty=1, lwd=2, c("Sordariomycetes","Xylariales","Xylaria","Dothideomycetes","Capnodiales"), bty="n", cex=0.9)
dev.off()

############################################################################################
############################################################################################
11. Alternating Logistic regression
############################################################################################
############################################################################################

data1<- data.frame(xylaria,sordariomycetes,dothideomycetes,xylariales,capnodiales,dis,PCA1, Fr_T)

ordgee0a = ordgee(ordered(sordariomycetes)~PCA1+Fr_T,id=Pl,trace=T,mean.link="logit",corstr=("userdefined"),z=z0a,data=data1)
ordgee1a = ordgee(ordered(xylaria)~PCA1+Fr_T,id=Pl,trace=T,mean.link="logit",corstr=("userdefined"),z=z0a,data=data1)
ordgee2a = ordgee(ordered(xylariales)~PCA1+Fr_T,id=Pl,trace=T, maxit = 35, mean.link="logit",corstr=("userdefined"),z=z0a,data=data1)
ordgee3a = ordgee(ordered(capnodiales)~PCA1+Fr_T,id=Pl,trace=T,mean.link="logit",corstr=("userdefined"),z=z0a,data=data1)
ordgee4a = ordgee(ordered(dothideomycetes)~PCA1+Fr_T,id=Pl,trace=T,mean.link="logit",corstr=("userdefined"),z=z0a,data=data1)

summary(ordgee0a)
summary(ordgee1a)
summary(ordgee2a)
summary(ordgee3a)
summary(ordgee4a)

############################################################################################
############################################################################################
12. Alternating Logistic regression_Inferior and Superior limits 
############################################################################################
############################################################################################

#### Confidence interval - 95% for (Odds ratio) the mean structure

cbind(exp(summary(ordgee0a)$mean[,1]-1.96*summary(ordgee0a)$mean[,2]), exp(summary(ordgee0a)$mean[,1]+1.96*summary(ordgee0a)$mean[,2]))
cbind(exp(summary(ordgee1a)$mean[,1]-1.96*summary(ordgee1a)$mean[,2]), exp(summary(ordgee1a)$mean[,1]+1.96*summary(ordgee1a)$mean[,2]))
cbind(exp(summary(ordgee2a)$mean[,1]-1.96*summary(ordgee2a)$mean[,2]), exp(summary(ordgee2a)$mean[,1]+1.96*summary(ordgee2a)$mean[,2]))
cbind(exp(summary(ordgee3a)$mean[,1]-1.96*summary(ordgee3a)$mean[,2]), exp(summary(ordgee3a)$mean[,1]+1.96*summary(ordgee3a)$mean[,2]))
cbind(exp(summary(ordgee4a)$mean[,1]-1.96*summary(ordgee4a)$mean[,2]), exp(summary(ordgee4a)$mean[,1]+1.96*summary(ordgee4a)$mean[,2]))

#### Confidence interval - 95% for (Odds ratio) the dependence structure

cbind(exp(summary(ordgee0a)$correlation[,1]-1.96*summary(ordgee0a)$correlation[,2]), exp(summary(ordgee0a)$correlation[,1]+1.96*summary(ordgee0a)$correlation[,2]))
cbind(exp(summary(ordgee1a)$correlation[,1]-1.96*summary(ordgee1a)$correlation[,2]), exp(summary(ordgee1a)$correlation[,1]+1.96*summary(ordgee1a)$correlation[,2]))
cbind(exp(summary(ordgee2a)$correlation[,1]-1.96*summary(ordgee2a)$correlation[,2]), exp(summary(ordgee2a)$correlation[,1]+1.96*summary(ordgee2a)$correlation[,2]))
cbind(exp(summary(ordgee3a)$correlation[,1]-1.96*summary(ordgee3a)$correlation[,2]), exp(summary(ordgee3a)$correlation[,1]+1.96*summary(ordgee3a)$correlation[,2]))
cbind(exp(summary(ordgee4a)$correlation[,1]-1.96*summary(ordgee4a)$correlation[,2]), exp(summary(ordgee4a)$correlation[,1]+1.96*summary(ordgee4a)$correlation[,2]))









