ll=load('manzello_angsd_dec4.RData')
head(meta)
me3=subset(me2,clones!="black")
me4=subset(me2,clones=="black")
me3$cn=as.numeric(as.factor(me3$clones))
me4$cn=seq(max(me3$cn)+1,max(me3$cn)+nrow(me4),1)
me5=rbind(me3,me4)
me5$cn=factor(me5$cn)
levels(me5$cn)
me5$May2016Condition[!(me5$May2016Condition %in% c("0","1","2","3"))]=NA

#----------------------
# are there more clones inshore?

me2$inshore=0
me2$inshore[grep("I",me2$Location)]=1

ci=nrow(subset(me2,me2$inshore==1 & me2$clones!="black")) # 81
co=nrow(subset(me2,me2$inshore==0 & me2$clones!="black")) # 14
ui=nrow(subset(me2,me2$inshore==1 & me2$clones=="black")) # 13
uo=nrow(subset(me2,me2$inshore==0 & me2$clones=="black")) # 72

clonesites=matrix(c(ci,ui,co,uo),nrow=2,ncol=2)
dimnames(clonesites)=list(c("clones","unique"),c("inshore","offshore"))
fisher.test(clonesites)
# data:  clonesites
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
 # 13.21693 79.41045
# sample estimates:
# odds ratio 
  # 31.07495 

#-----------------------

me6$fracZ[is.na(me6$fracZ)]=1e-5
me6$Location=factor(me6$Location,levels=c("LKI1","LKI2","UKI1","UKI2","UKI3","LKO1","LKO2","UKO1","UKO2","UKO3"))
# centered log10-fracZ
me6$lfz=log(me6$fracZ,10)
me6$lfz=me6$lfz-mean(me6$lfz)

me6$bleach=1
me6$bleach[me6$Condition..P.PALE..PB.PARTIALLY.BLEACHED..H.HEALTHY..B.BLEACHED.=="PB"]=2
me6$bleach[me6$Condition..P.PALE..PB.PARTIALLY.BLEACHED..H.HEALTHY..B.BLEACHED.=="P"]=3
me6$bleach[me6$Condition..P.PALE..PB.PARTIALLY.BLEACHED..H.HEALTHY..B.BLEACHED.=="H"]=4
me6[,c("Condition..P.PALE..PB.PARTIALLY.BLEACHED..H.HEALTHY..B.BLEACHED.","bleach")]
me6$bleach=as.factor(me6$bleach)
table(me6$bleach)
me6$inshore=0
me6$inshore[grep("O",me6$Location)]=1
me6$isclone=1
me6$isclone[me6$clones=="black"]=0
head(me6,20)

save(me6,file='~/Dropbox/Documents/manzello/angsd/clones_meta_me6.RData')
write.csv(me6,quote=F, row.names=F,file='~/Dropbox/Documents/manzello/angsd/clones_meta_me6.csv')

#----------------------------
# is bleaching different across locations? (LRT)

library(lme4)
lr=lmer(lfz~Location+(1|cn),me6)
lrc=lmer(lfz~isclone+Location+(1|cn),me6)
lr0=lmer(lfz~(1|cn),me6)
anova(lr,lrc,lr0)
# lr0: lfz ~ (1 | cn)
# lr: lfz ~ Location + (1 | cn)
# lrc: lfz ~ isclone + Location + (1 | cn)
    # Df    AIC    BIC  logLik deviance   Chisq Chi Df Pr(>Chisq)    
# lr0  3 319.58 329.15 -156.79   313.58                              
# lr  12 290.93 329.24 -133.46   266.93 46.6507      9  4.558e-07 ***
# lrc 13 292.92 334.43 -133.46   266.92  0.0018      1      0.966    

library(ordinal)
cl=clmm2(location = bleach ~ Location, random = cn, data = subset(me6,Location!="UKI2"),Hess=T)
clc=clmm2(location = bleach ~ Location +isclone, random = cn, data = subset(me6,Location!="UKI2"),Hess=T)
cl0=clmm2(location = bleach ~1, random = cn, data = subset(me6,Location!="UKI2"),Hess=T)
anova(cl,clc,cl0)
# Response: bleach
                     # Model Resid. df -2logLik   Test    Df   LR stat.      Pr(Chi)
# 1                  1 |  |        156 376.0614                                     
# 2           Location |  |        148 327.9381 1 vs 2     8 48.1233209 9.357423e-08
# 3 Location + isclone |  |        147 327.2943 2 vs 3     1  0.6437501 4.223561e-01

#---------------------
# is bleaching higher or lower inshore? (accounting for clone effects)

# LRT
library(lme4)
lr=lmer(lfz~inshore+(1|cn),me6)
lrc=lmer(lfz~isclone+inshore+(1|cn),me6)
lr0=lmer(lfz~(1|cn),me6)
anova(lr,lrc,lr0)
# lr0  3 319.58 329.15 -156.79   313.58                           
# lr   4 317.16 329.93 -154.58   309.15 4.4206      1    0.03551 *
# lrc  5 318.96 334.92 -154.48   308.96 0.1967      1    0.65742  
lr=lmer(lfz~isclone+(1|cn),me5)
lr0=lmer(lfz~(1|cn),me5)
anova(lr,lr0)

mcb=MCMCglmm(lfz~inshore,random=~cn,data=me5,prior=prior,nitt=55000,thin=50,burnin=5000,pr=T)
summary(mcb)
# inshore      -0.33861 -0.66457 -0.03487     1000 0.044 *

library(ordinal)
cl=clmm2(location = bleach ~ inshore, random = cn, data = subset(me6,Location!="UKI2"),Hess=T)
clc=clmm2(location = bleach ~ inshore+isclone, random = cn, data = subset(me6,Location!="UKI2"),Hess=T)
cl0=clmm2(location = bleach ~1, random = cn, data = subset(me6,Location!="UKI2"),Hess=T)
anova(cl,clc,cl0)
# Response: bleach
                    # Model Resid. df -2logLik   Test    Df   LR stat.    Pr(Chi)
# 1                 1 |  |        156 376.0614                                   
# 2           inshore |  |        155 372.0100 1 vs 2     1 4.05143081 0.04413392
# 3 inshore + isclone |  |        154 371.9199 2 vs 3     1 0.09001768 0.76415468


# ordinal model for inshore-offshore bleaching
prioro = list(R = list(V = 1, nu=0.002,fix=1),G=list(G1=list(V=1, nu = 0.002)))
mcbo=MCMCglmm(bleach~inshore,random=~cn,data=me5,prior=prioro,family="ordinal",nitt=155000,thin=100,burnin=10000,pr=T)
summary(mcbo)
# inshore        -2.126   -4.225   -0.228   1450.0 0.0331 *  

# plotting bleaching bar chart
par(las=2,mgp=c(2.1,1,0))
plot(as.factor(5-as.numeric(bleach))~Location, me6,xlab="",ylab="bleaching",col=c("darkorange4","darkorange3","burlywood1","floralwhite"))
me6$May2016Condition=factor(me6$May2016Condition, levels=unique(me6$May2016Condition))
plot(May2016Condition~Location, me6,xlab="",ylab="recovery",col=c("darkorange4","darkorange3","burlywood1"))

sols=data.frame(summary(mc)$solutions)
sols$location=sub("Location","",row.names(sols))
sols$location=factor(sols$location,levels=c("LKI1","LKI2","UKI1","UKI2","UKI3","LKO1","LKO2","UKO1","UKO2","UKO3"))
quartz()
library(ggplot2)
ggplot(sols,aes(location,post.mean))+geom_point()+geom_errorbar(aes(ymin=l.95..CI,ymax=u.95..CI),lwd=0.4,width=0.4)+ylab("log10(zooxanthellae)")+theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust = 1))+ylim(min(sols$l.95..CI),max(sols$u.95..CI)+0)

#------------------------------
# broad-sense heritability of bleaching tolerance

library(MCMCglmm)
ll=load('~/Dropbox/Documents/manzello/angsd/clones_meta_me6.RData')
ll=load('~/Dropbox/Documents/manzello/colony_distances.RData')
me6=me6[!is.na(me6$fracD),]
goods=intersect(colnames(mapd),me6$Tag.Number)
mapd=mapd[goods,goods]

dd=matrix(0,ncol=nrow(me6),nrow=nrow(me6))
for (c1 in 1:(nrow(me6)-1)) {
	for (c2 in (c1+1):nrow(me6)) {
			dd[c1,c2]=round(abs(asin(sqrt(me6$fracD[c1]))-asin(sqrt(me6$fracD[c2]))),3)
	}
}
dd=dd+t(dd)
plot(hclust(as.dist(dd),"ave"))

head(me6)
dc=matrix(0,ncol=nrow(me6),nrow=nrow(me6))
for (c1 in 1:(nrow(me6)-1)) {
	for (c2 in (c1+1):nrow(me6)) {
			dc[c1,c2]=as.numeric(me6$cn[c1]==me6$cn[c2])
	}
}
dc=dc+t(dc)
plot(hclust(as.dist(dc),"ave"))

dimnames(dc)=list(me6$Tag.Number,me6$Tag.Number)
dimnames(dd)=list(me6$Tag.Number,me6$Tag.Number)

data(bird.families)
str(bird.families)
?inverseA
1-mapd

head(me6)

me5=me6


# With fixed effect of location, based on reads
prior = list(R = list(V = 1, nu = 0.002),G=list(G1=list(V=1, nu = 0.002,alpha.mu=0,alpha.V=1000)))
mc=MCMCglmm(lfz~0+Location,random=~cn,data=me6,prior=prior,nitt=155000,thin=50,burnin=10000,pr=T)
summary(mc$Sol)
bigH=c()
bigH=(mc$VCV[,1])/(mc$VCV[,1]+mc$VCV[,2])
mean(bigH) # 0.71
HPDinterval(bigH) # 0.56-0.82 

cor.test(me6$lfz,asin(sqrt(me6$fracD)))


# ordinal model with fixed residual var (omitting UKI2 as it does not mix well)
prior = list(R = list(V = 1, nu=0.002,fix=1),G=list(G1=list(V=1, nu = 0.002)))
me6.s=subset(me6,Location!="UKI2") 
bb=MCMCglmm(bleach~0+Location,random=~cn,data=me6.s,family="ordinal",prior=prior,nitt=510000,thin=500,burnin=10000,pr=T)
#plot(bb)
summary(bb)

me6.s$ordp=predict(bb,marginal=NULL)
plot(ordp~bleach,me6.s)

bigH=c()
summary(bb$VCV)
colnames(bb$VCV)
bigH=(bb$VCV[,1])/(bb$VCV[,1]+bb$VCV[,2])
mean(bigH) # 0.93 
HPDinterval(bigH) # 0.87-0.98 

solsbb=data.frame(summary(bb)$solutions)
solsbb$location=sub("Location","",row.names(solsbb))
solsbb$location=factor(solsbb$location,levels=c("LKI1","LKI2","UKI1","UKI3","LKO1","LKO2","UKO1","UKO2","UKO3"))

#quartz()
library(ggplot2)
ggplot(solsbb,aes(location,post.mean/log(10)))+geom_point()+geom_errorbar(aes(ymin=l.95..CI/log(10),ymax=u.95..CI/log(10)),lwd=0.4,width=0.4)+ylab("log10 (odds of \nhigher color score )")+theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust = 1))

quartz()
bbx=solsbb[,1]/log(10)
par(las=0)
plot(sols[-4,1]~bbx,mgp=c(2.1,1,0),bty="n",xlab="color",ylab="reads", main="model coefficients",xlim=c(-2,3.2),ylim=c(-1.4,0.8))
abline(v=0,lty=3)
abline(h=0,lty=3)
abline(lm(sols[-4,1]~bbx),lty=3,col="red")
summary(lm(sols[-4,1]~bbx))  # R2=0.8 P=0.0006
cor.test(bbx,sols[-4,1]) # 0.91 ( 0.6265425 0.9814666 )
mtext("R = 0.91, P < 0.001",cex=0.9)

#--------------------
# heritability of fraction of D

hist(asin(sqrt(me5$fracD)))

mc=MCMCglmm(asin(sqrt(fracD))~1,random=~cn,data=me6,nitt=55000,thin=50,burnin=5000)
summary(mc)
bigH=c()
colnames(mc$VCV)
bigH=(mc$VCV[,1])/(mc$VCV[,1]+mc$VCV[,2])
mean(bigH) # 0.77 
HPDinterval(bigH) # 0.67-0.86 

mc=MCMCglmm(asin(sqrt(fracD))~Location,random=~cn,data=me6,nitt=55000,thin=50,burnin=5000)
summary(mc)
bigH=c()
bigH=(mc$VCV[,1])/(mc$VCV[,1]+mc$VCV[,2])
mean(bigH) # 0.73 
HPDinterval(bigH) # 0.61-0.84 

#-------------------------
# heritability in double-response model

# standartizing for double-response model
me5=me6
me5$zstandard=(mean(log(me5$fracZ))-log(me5$fracZ))/sd(log(me5$fracZ))
me5$dstandard=(mean(asin(sqrt(me5$fracD)),na.rm=T)-asin(sqrt(me5$fracD)))/sd(asin(sqrt(me5$fracD)),na.rm=T)
plot(density(me5$dstandard,na.rm=T))

library(MCMCglmm)
prior = list(R = list(V = diag(2), nu = 2-0.998),G=list(G1=list(V=diag(2), nu = 2-0.998)))
mc=MCMCglmm(cbind(dstandard,zstandard)~Location,random=~us(trait):cn,rcov=~idh(trait):units,data=me5,family=c("gaussian","gaussian"),nitt=155000,thin=50,burnin=5000,prior=prior)

summary(mc)
plot(mc)
bigH=c()
summary(mc$VCV)
colnames(mc$VCV)
vcv=matrix(colMeans(mc$VCV[,1:4]),2,2)
colnames(vcv)=colnames(mc$VCV[,1:2])
head(mc$VCV)

HPDinterval((mc$VCV[,1]+mc$VCV[,4])/rowSums(mc$VCV[,c(1,4,5,6)]))
# var1 0.660163 0.8227876
mean((mc$VCV[,1]+mc$VCV[,4])/rowSums(mc$VCV[,c(1,4,5,6)]))
# 0.747

# do we calculate double-trait heritability correctly?
colnames(mc$VCV)
bigH=(mc$VCV[,1]+mc$VCV[,4])/(mc$VCV[,1]+mc$VCV[,4]+mc$VCV[,5]+mc$VCV[,6])
mean(bigH) # 0.73 
HPDinterval(bigH) # 0.63 0.81 
# for D fraction
bigHd=(mc$VCV[,1])/(mc$VCV[,1]+mc$VCV[,5])
mean(bigHd) # 0.74 
HPDinterval(bigHd) # 0.62 0.85 
# for zoox proportion
bigHz=(mc$VCV[,4])/(mc$VCV[,4]+mc$VCV[,6])
mean(bigHz) # 0.71
HPDinterval(bigHz) # 0.57 0.82 

corrs=cov2cor(matrix(apply(mc$VCV[,1:4],2,mean),nrow=2))
corrs

bigH=(mc$VCV[,1]+mc$VCV[,2]+mc$VCV[,3]+mc$VCV[,4])/(mc$VCV[,1]+mc$VCV[,2]+mc$VCV[,3]+mc$VCV[,4]+mc$VCV[,5]+mc$VCV[,6])
mean(bigH) # 0.82 
HPDinterval(bigH) # 0.75 0.88 

colnames(mc$VCV[,1:4])
# fraction of variation explained by clone effect on covariance:
covar=(mc$VCV[,2]+mc$VCV[,3])/(mc$VCV[,1]+mc$VCV[,2]+mc$VCV[,3]+mc$VCV[,4]+mc$VCV[,5]+mc$VCV[,6])
mean(covar) # 0.28 
HPDinterval(covar) # 0.19 0.36 


library(MCMC.qpcr)
mcmc.pval(mc$VCV[,"traitzstandard:traitdstandard.cn"])
# 0.0005697883

covs=summary(mc)$Gcovariances  
                                 # post.mean  l-95% CI  u-95% CI eff.samp
# traitdstandard:traitdstandard.cn 0.8610487 0.5096670 1.2384532     3000
# traitzstandard:traitdstandard.cn 0.4024957 0.1840186 0.6311813     3000
# traitdstandard:traitzstandard.cn 0.4024957 0.1840186 0.6311813     3000
# traitzstandard:traitzstandard.cn 0.6593045 0.4088476 0.9284492     3000
gmat=matrix(covs,nrow=2,ncol=2)
          # [,1]      [,2]
# [1,] 0.8610487 0.4024957
# [2,] 0.4024957 0.6593045
cov2cor(gmat)
          # [,1]      [,2]
# [1,] 1.0000000 0.5342006
# [2,] 0.5342006 1.0000000

cor(me5$zstandard,me5$dstandard,use="complete")
# 0.52
cor.test(me5$zstandard,me5$dstandard,use="complete")
# t = 8.0662, df = 173, p-value = 1.161e-13
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
 # 0.4059184 0.6228245
# sample estimates:
      # cor 
# 0.5227821 

#--------------------
# recovery by site, based on visual scores

library(MCMCglmm)

me6.s$May2016Condition=factor(me6.s$May2016Condition,levels=unique(me6.s$May2016Condition))
prior = list(R = list(V = 1, fix=1),G=list(G1=list(V=1, nu = 0.002)))
bbr=MCMCglmm(May2016Condition~0+Location,random=~cn,data=me6,family="ordinal",prior=prior,nitt=55000,thin=50,burnin=5000)
plot(bbr)
summary(bbr)
mm=abs(mctukey(bb,prob=0.95))
dimnames(mm)=list(sub("Location","",dimnames(mm)[[1]]),sub("Location","",dimnames(mm)[[1]]))
mm[mm>0]=1
plot(hclust(as.dist(mm),method="complete"))

sols=data.frame(summary(bb)$solutions)
sols$location=sub("Location","",row.names(sols))
quartz()

ggplot(sols,aes(location,post.mean))+geom_point()+geom_errorbar(aes(ymin=l.95..CI,ymax=u.95..CI),lwd=0.4,width=0.4)+ylab("bleaching tolerance")+theme_bw()

#-----------------------
# do clones bleach more/less?

# ordinal model with fixed residual var (omitting UKI2 as it does not mix well)
prior = list(R = list(V = 1, nu=0.002,fix=1),G=list(G1=list(V=1, nu = 0.002)))
me6.s=subset(me6,Location!="UKI2") 
names(me6.s)
me6.s$isclonal=as.numeric(me6.s$clones!="black")

# visual scores, effect of location
prior = list(R = list(V = 1, nu=0.002,fix=1),G=list(G1=list(V=1, nu = 0.002)))
bbc=MCMCglmm(bleach~0+Location+isclonal,random=~cn,data=me6.s,family="ordinal",prior=prior,nitt=55000,thin=50,burnin=5000,pr=T)
plot(bbc)
summary(bbc)
# isclonal       -0.9599  -3.9351   1.8357    633.8  0.526    

# visual scores, no random clone effect
prior2 = list(R = list(V = 1, nu=0.002,fix=1))
bbc2=MCMCglmm(bleach~0+Location+isclonal,data=me6.s,family="ordinal",prior=prior2,nitt=55000,thin=50,burnin=5000,pr=T)
plot(bbc2)
summary(bbc2)
# isclonal       -0.2451  -0.9522   0.4938   1000.0  0.470    

# reads
me5$isclonal=as.numeric(me5$clones!="black")
prior = list(R = list(V = 1, nu = 0.002),G=list(G1=list(V=1, nu = 0.002,alpha.mu=0,alpha.V=1000)))
mc=MCMCglmm(lfz~0+Location+isclonal,random=~cn,data=me5,prior=prior,nitt=55000,thin=50,burnin=5000,pr=T)
summary(mc)
# isclonal      0.0073990 -0.3594471  0.3307057   1000.0  0.954    

# reads, no random clone effect
prior2 = list(R = list(V = 1, nu=0.002,fix=1))
me5$isclonal=as.numeric(me5$clones!="black")
mc2=MCMCglmm(lfz~0+Location+isclonal,data=me5,prior=prior2,nitt=55000,thin=50,burnin=5000,pr=T)
summary(mc2)
# isclonal      -0.02059 -0.41641  0.44520   1000.0  0.956    


table(me6.s$May2016Condition)
me6.s$May2016Condition[me6.s$May2016Condition==1]=2
me6.s$May2016Condition=factor(me6.s$May2016Condition,levels=unique(me6.s$May2016Condition))
prior = list(R = list(V = 1, nu=0.002,fix=1),G=list(G1=list(V=1, nu = 0.002)))
bbc=MCMCglmm(May2016Condition~0+Location+isclonal,random=~cn,data=me6.s,family="categorical",prior=prior,nitt=1555000,thin=1500,burnin=15000,pr=T)
plot(bbc)
summary(bbc)
# isclonal          8197   -30659    47873   53.672  0.628    



#--------------------
# bleaching status heritability based on visual scores

library(MCMCglmm)
ll=load('~/Dropbox/Documents/manzello/angsd/clones_meta_me6.RData')


library(ordinal)
cl=clmm2(location = bleach ~ Location, random = cn, data = me6,Hess=T)
summary(cl)

library(MCMCglmm)
prior = list(R = list(V = 1, fix=1),G=list(G1=list(V=1, nu = 0.02)))
bb=MCMCglmm(bleach~0+Location,random=~cn,data=me6,family="ordinal",prior=prior,nitt=510000,thin=500,burnin=10000)
plot(bb)
summary(bb)
str(bb$CP)
mctukey=function(model,prob=0.95) {
	require(MCMCglmm)
	diffs=c()
	for (l in colnames(model$Sol)){
		for (l2 in colnames(model$Sol)){
			if (l==l2) { 
				diffs=append(diffs,0)
				next
			}
			d= model$Sol[,l]-model$Sol[,l2]
			hpd=HPDinterval(d,prob=prob)
			if (prod(hpd)>0) {
				diffs=append(diffs,mean(d)) 
			} else { diffs=append(diffs,0) }
		 }
	}
	m=matrix(diffs,ncol=ncol(model$Sol),nrow=ncol(model$Sol))
	dimnames(m)=list(colnames(model$Sol),colnames(model$Sol))
	return(m)
}

mm=abs(mctukey(bb,prob=0.95))
dimnames(mm)=list(sub("Location","",dimnames(mm)[[1]]),sub("Location","",dimnames(mm)[[1]]))
mm[mm>0]=1
plot(hclust(as.dist(mm),method="complete"))

sols=data.frame(summary(bb)$solutions)
sols$location=sub("Location","",row.names(sols))
sols=subset(sols,location!="UKI2")
sols$location=factor(sols$location,levels=c("LKI1","LKI2","UKI1","UKI3","LKO1","LKO2","UKO1","UKO2","UKO3"))
quartz()
ggplot(sols,aes(location,post.mean))+geom_point()+geom_errorbar(aes(ymin=l.95..CI,ymax=u.95..CI),lwd=0.4,width=0.4)+ylab("bleaching tolerance")+theme_bw()

plot(sols$post.mean~solsbb$post.mean)

bigH=c()
bigH=(bb$VCV[,1])/(bb$VCV[,1]+bb$VCV[,2])
mean(bigH) # 0.78
HPDinterval(bigH) # 0.58-0.91

#-------------------
# sanity check for ordinal model
library(Rmisc)

me6$bl=as.numeric(as.character(me6$bleach))
su=summarySE(me6,measurevar="bl",groupvars="Location")
me6$bleach[me6$Location=="LKI1"]


head(me6)
# do clones recover better than non-clones, overall?
me2=me6
head(me2)
ci=nrow(subset(me2,me2$May2016Condition==3 & me2$clones!="black")) # 89
co=nrow(subset(me2,me2$May2016Condition!=3 & me2$clones!="black")) # 6
ui=nrow(subset(me2,me2$May2016Condition==3 & me2$clones=="black")) # 62
uo=nrow(subset(me2,me2$May2016Condition!=3 & me2$clones=="black")) # 22

clonesites=matrix(c(ci,ui,co,uo),nrow=2,ncol=2)
dimnames(clonesites)=list(c("clones","unique"),c("inshore","offshore"))
fisher.test(clonesites)

# data:  clonesites
# p-value = 0.0003346
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
  # 1.909038 16.663077
# sample estimates:
# odds ratio 
  # 5.215836 

# do inshore reefs bleach less than offshores?
me2i=subset(me2,inshore==0)
ci=nrow(subset(me2i,me2i$May2016Condition==3 & me2i$clones!="black")) # 8
co=nrow(subset(me2i,me2i$May2016Condition!=3 & me2i$clones!="black")) # 6
ui=nrow(subset(me2i,me2i$May2016Condition==3 & me2i$clones=="black")) # 49
uo=nrow(subset(me2i,me2i$May2016Condition!=3 & me2i$clones=="black")) # 22
clonesites=matrix(c(ci,ui,co,uo),nrow=2,ncol=2)
fisher.test(clonesites)


# do inshore clones recover better than non-clones?
me2i=subset(me2,inshore==1)
ci=nrow(subset(me2i,me2i$May2016Condition==3 & me2i$clones!="black")) # 81
co=nrow(subset(me2i,me2i$May2016Condition!=3 & me2i$clones!="black")) # 0
ui=nrow(subset(me2i,me2i$May2016Condition==3 & me2i$clones=="black")) # 13
uo=nrow(subset(me2i,me2i$May2016Condition!=3 & me2i$clones=="black")) # 0
# all fully recovered

# do inshore clones recover better than non-clones?
me2i=subset(me2,inshore==0)
ci=nrow(subset(me2i,me2i$May2016Condition==3 & me2i$clones!="black")) # 8
co=nrow(subset(me2i,me2i$May2016Condition!=3 & me2i$clones!="black")) # 6
ui=nrow(subset(me2i,me2i$May2016Condition==3 & me2i$clones=="black")) # 49
uo=nrow(subset(me2i,me2i$May2016Condition!=3 & me2i$clones=="black")) # 22
clonesites=matrix(c(ci,ui,co,uo),nrow=2,ncol=2)
fisher.test(clonesites)
# data:  clonesites
# p-value = 0.5346
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
 # 0.1604163 2.3757583
# sample estimates:
# odds ratio 
 # 0.6024549 
# nope.

# do inshore corals recover better than offshore?
ci=nrow(subset(me2,me2$May2016Condition==3 & me2$inshore==1)) # 94
co=nrow(subset(me2,me2$May2016Condition!=3 & me2$inshore==1)) # 0
ui=nrow(subset(me2,me2$May2016Condition==3 & me2$inshore==0)) # 57
uo=nrow(subset(me2,me2$May2016Condition!=3 & me2$inshore==0)) # 28
clonesites=matrix(c(ci,ui,co,uo),nrow=2,ncol=2)
fisher.test(clonesites)
# data:  clonesites
# p-value = 5.373e-11
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
 # 11.00616      Inf
# sample estimates:
# odds ratio 
       # Inf 


head(me2)
# do inshore corals bleach less than offshore?
ci=nrow(subset(me2,me2$May2016Condition==3 & me2$inshore==1)) # 94
co=nrow(subset(me2,me2$May2016Condition!=3 & me2$inshore==1)) # 0
ui=nrow(subset(me2,me2$May2016Condition==3 & me2$inshore==0)) # 57
uo=nrow(subset(me2,me2$May2016Condition!=3 & me2$inshore==0)) # 28
clonesites=matrix(c(ci,ui,co,uo),nrow=2,ncol=2)
fisher.test(clonesites)
# data:  clonesites
# p-value = 5.373e-11
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
 # 11.00616      Inf
# sample estimates:
# odds ratio 
       # Inf 
       
names(me6)[14]="visual"
me6$visual=factor(me6$visual,levels=c("B","PB","P","H"))
library(ggplot2)
ggplot(me6,aes(x=visual,y=lfz))+geom_violin(draw_quantiles=0.5)+theme_bw()
plot(density(me6$lfz))
plot(lfz~visual,me6)
summary(lm(lfz~visual,me6[me6$visual!="B",]))
