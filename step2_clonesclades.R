#----------------------- preparing metadata 

ll=load('clones_meta_me6.RData')

# computing standardized log10-transformed zooxanthellae abundance:
me6$fracZ[is.na(me6$fracZ)]=1e-5
me6$lfz=log(me6$fracZ,10)
me6$lfz=me6$lfz-mean(me6$lfz)

# releveling locations
me6$Location=factor(me6$Location,levels=c("LKI1","LKI2","UKI1","UKI2","UKI3","LKO1","LKO2","UKO1","UKO2","UKO3"))

# converting bleaching scores to ordered categories
me6$bleach=1
me6$bleach[me6$Condition..P.PALE..PB.PARTIALLY.BLEACHED..H.HEALTHY..B.BLEACHED.=="PB"]=2
me6$bleach[me6$Condition..P.PALE..PB.PARTIALLY.BLEACHED..H.HEALTHY..B.BLEACHED.=="P"]=3
me6$bleach[me6$Condition..P.PALE..PB.PARTIALLY.BLEACHED..H.HEALTHY..B.BLEACHED.=="H"]=4
me6[,c("Condition..P.PALE..PB.PARTIALLY.BLEACHED..H.HEALTHY..B.BLEACHED.","bleach")]
me6$bleach=as.factor(me6$bleach)
table(me6$bleach)

# inshore vs offshore
me6$inshore=0
me6$inshore[grep("I",me6$Location)]=1

# is a coral member of a conal group (i.e. shares genotype with at least one other sample) 
me6$isclone=1
me6$isclone[me6$clones=="black"]=0

save(me6,file='clones_meta_me6.RData')

# ----------- comparing visual scores and read-based zoox abundances

# releveling visual scores
me6$Condition..P.PALE..PB.PARTIALLY.BLEACHED..H.HEALTHY..B.BLEACHED.=factor(me6$Condition..P.PALE..PB.PARTIALLY.BLEACHED..H.HEALTHY..B.BLEACHED.,levels=c("B","PB","P","H"))

plot(lfz~Condition..P.PALE..PB.PARTIALLY.BLEACHED..H.HEALTHY..B.BLEACHED.,me6,xlab="visual bleaching", ylab="log ( zoox reads )")

# correlation with bleacing categories
ct=cor.test(me6$lfz,as.numeric(me6$bleach))
ct
	# Pearson's product-moment correlation

# data:  me6$lfz and as.numeric(me6$bleach)
# t = 11.833, df = 178, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
 # 0.5728810 0.7381674
# sample estimates:
      # cor 
# 0.6635449 
mtext(side=3,paste(" r = ",round(ct$estimate,2)))

TukeyHSD(aov(lm(lfz~Condition..P.PALE..PB.PARTIALLY.BLEACHED..H.HEALTHY..B.BLEACHED.,me6)))
          # diff         lwr       upr     p adj
# PB-B 0.9676728  0.70544392 1.2299016 0.0000000
# P-B  1.1201358  0.77510976 1.4651619 0.0000000
# H-B  1.2461191  1.01372251 1.4785156 0.0000000
# P-PB 0.1524631 -0.20207545 0.5070016 0.6805449
# H-PB 0.2784463  0.03214824 0.5247444 0.0197584
# H-P  0.1259832 -0.20709598 0.4590625 0.7604281

# lfz does not discriminate Pale from Healthy or 2 Part-Bleached
# it does discriminate Bleached from everyone else and Healthy from Part-Bleached

# ----------------------- Correlation of read-based bleaching and fraction of Durusdinium

cor.test(me6$lfz,asin(sqrt(me6$fracD)))
	# Pearson's product-moment correlation

# data:  me6$lfz and asin(sqrt(me6$fracD))
# t = 8.0662, df = 173, p-value = 1.161e-13
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
 # 0.4059184 0.6228245
# sample estimates:
      # cor 
# 0.5227821 

plot(me6$lfz~asin(sqrt(me6$fracD)),pch=16,col=rgb(0,0,0,alpha=0.3),xlab="fraction of Durusdinium",ylab="log (zoox reads)")

# --------------------- Plotting bleaching bar chart

par(mgp=c(2.1,1,0),mfrow=c(1,2))
plot(as.factor(5-as.numeric(bleach))~Location, me6,xlab="",ylab="bleaching",col=c("darkorange4","darkorange3","burlywood1","floralwhite"))
me6$May2016Condition=factor(me6$May2016Condition, levels=unique(me6$May2016Condition))
plot(May2016Condition~Location, me6,xlab="",ylab="recovery",col=c("darkorange4","darkorange3","burlywood1"))

# --------------------- clones by location bar chart

library(ggplot2)

me6$Location=factor(me6$Location,levels=c("LKI1","LKI2","UKI1","UKI2","UKI3","LKO1","LKO2","UKO1","UKO2","UKO3"))
bysite=data.frame(cbind("reef"=as.character(me6$Location),"clone"=as.character(me6$clones)))
cols=sort(as.character(unique(bysite$clone)))
bysite$clone=as.character(bysite$clone)
bysite$clone[bysite$clone=="black"]="notClonal"
bysite$clone[bysite$clone=="black"]="notClonal"
cols=cols[-1]
bysite$clone=factor(bysite$clone,levels=c(cols,"notClonal"))
cols=append(cols,"grey90")
bysite$reef=factor(bysite$reef,levels=c("LKI1","LKI2","UKI1","UKI2","UKI3","LKO1","LKO2","UKO1","UKO2","UKO3"))
ggplot(bysite,aes(reef))+geom_bar(aes(fill=clone))+scale_fill_manual(values=cols)+theme_bw()

#----------------------- Are there more clones inshore?

# counting how many corals are members of clonal groups inshore vs offshore:
ci=nrow(subset(me6,me6$inshore==0 & me6$clones!="black")) # 81
co=nrow(subset(me6,me6$inshore==1 & me6$clones!="black")) # 14
ui=nrow(subset(me6,me6$inshore==0 & me6$clones=="black")) # 13
uo=nrow(subset(me6,me6$inshore==1 & me6$clones=="black")) # 72

clonesites=matrix(c(ci,ui,co,uo),nrow=2,ncol=2)
dimnames(clonesites)=list(c("clones","unique"),c("inshore","offshore"))
fisher.test(clonesites)
# data:  clonesites
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
 # 0.01259280 0.07566051
# sample estimates:
# odds ratio 
# 0.03218026 


#----------------------- Is bleaching different across locations? 

# ( cn = genet id ; random effect)

# Likelihood Ratio Test, based on standardized log-transformed zoox read counts (lfz)

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

# ordinal model based on visual bleaching scores (omitting UKI2 because there is a single genet at that site)
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

# ==> Bleaching varies significantly by location, but whether a coral is a member of a clonal group does not matter 
# (i.e. clones are not inherently more bleaching-resistant) .

# ---------------- Effect of location on bleaching (based on reads and visual scores)

library(ggplot2)
library(MCMCglmm)

# based on reads

prior = list(R = list(V = 1, nu = 0.002),G=list(G1=list(V=1, nu = 0.002,alpha.mu=0,alpha.V=1000)))
mc=MCMCglmm(lfz~0+Location,random=~cn,data=me6,prior=prior,nitt=55000,thin=50,burnin=5000,pr=T)
solsr=data.frame(summary(mc)$solutions)
solsr$location=sub("Location","",row.names(solsr))
solsr$location=factor(solsr$location,levels=c("LKI1","LKI2","UKI1","UKI2","UKI3","LKO1","LKO2","UKO1","UKO2","UKO3"))
frame()
ggplot(solsr,aes(location,post.mean))+geom_point()+geom_errorbar(aes(ymin=l.95..CI,ymax=u.95..CI),lwd=0.4,width=0.4)+ylab("log10(zooxanthellae)")+theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust = 1))

# based on bleaching scores: ordinal model with fixed residual var 

# (omitting UKI2 - single-genet site - as the model does not mix well)
me6.s=subset(me6,Location!="UKI2") 

prior = list(R = list(V = 1, nu=0.002,fix=1),G=list(G1=list(V=1, nu = 0.002)))
bb=MCMCglmm(bleach~0+Location,random=~cn,data=me6.s,family="ordinal",prior=prior,nitt=55000,thin=50,burnin=5000,pr=T)
summary(bb)
solso=data.frame(summary(bb)$solutions)
solso$location=sub("Location","",row.names(solso))
solso$location=factor(solso$location,levels=c("LKI1","LKI2","UKI1","UKI3","LKO1","LKO2","UKO1","UKO2","UKO3"))
frame()
ggplot(solso,aes(location,post.mean))+geom_point()+geom_errorbar(aes(ymin=l.95..CI,ymax=u.95..CI),lwd=0.4,width=0.4)+ylab("log10(zooxanthellae)")+theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust = 1))
ordp=predict(bb,marginal=NULL)

# similarity of location effects between model based on reads and prediction
quartz()
bbx=solso[,1]/log(10)
par(las=0)
plot(solsr[-4,1]~bbx,mgp=c(2.1,1,0),bty="n",xlab="color",ylab="reads", main="model coefficients",xlim=c(-2,3.2),ylim=c(-1.4,0.8))
abline(v=0,lty=3)
abline(h=0,lty=3)
abline(lm(solsr[-4,1]~bbx),lty=3,col="red")
summary(lm(solsr[-4,1]~bbx))  # R2=0.8 P=0.0006
cor.test(bbx,solsr[-4,1]) # 0.91 ( 0.6265425 0.9814666 )
mtext("R = 0.91, P < 0.001",cex=0.9)

#--------------------- Is bleaching different between inshore an offshore? 

# (accounting for random clone effects, cn)

# LRT, two ways: lme4 and MCMCglmm packages

library(lme4)
lr=lmer(lfz~inshore+(1|cn),me6)
lrc=lmer(lfz~isclone+inshore+(1|cn),me6)
lr0=lmer(lfz~(1|cn),me6)
anova(lr,lrc,lr0)
# lr0  3 319.58 329.15 -156.79   313.58                           
# lr   4 317.16 329.93 -154.58   309.15 4.4206      1    0.03551 *
# lrc  5 318.96 334.92 -154.48   308.96 0.1967      1    0.65742  

library(MCMCglmm)
prior = list(R = list(V = 1, nu = 0.002),G=list(G1=list(V=1, nu = 0.002,alpha.mu=0,alpha.V=1000)))
mcb=MCMCglmm(lfz~inshore,random=~cn,data=me6,prior=prior,nitt=55000,thin=50,burnin=5000,pr=T)
summary(mcb)
# inshore       0.34632  0.02814  0.66964     1000 0.044 * 

# ordinal GLM for inshore-offshore bleaching (based on visual scores) - two ways again

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


prioro = list(R = list(V = 1, nu=0.002,fix=1),G=list(G1=list(V=1, nu = 0.002)))
mcbo=MCMCglmm(bleach~inshore,random=~cn,data=me6,prior=prioro,family="ordinal",nitt=55000,thin=50,burnin=5000,pr=T)
summary(mcbo)
# inshore        2.1517   0.1062   4.4155    763.8 0.0234 *


#-------------------------- Broad-sense heritability of bleaching tolerance

# (proportion of bleaching variation attributable to differences between clones)

library(MCMCglmm)
library(ggplot2)

ll=load('clones_meta_me6.RData')

# Read-based, With fixed effect of location
prior = list(R = list(V = 1, nu = 0.002),G=list(G1=list(V=1, nu = 0.002,alpha.mu=0,alpha.V=1000)))
mc=MCMCglmm(lfz~0+Location,random=~cn,data=me6,prior=prior,nitt=155000,thin=50,burnin=10000,pr=T)
summary(mc$Sol)
bigH=c()
bigH=(mc$VCV[,1])/(mc$VCV[,1]+mc$VCV[,2])
mean(bigH) # 0.73
HPDinterval(bigH) # 0.60-0.84 

#-------------------- Heritability of fraction of D

prior = list(R = list(V = 1, nu = 0.002),G=list(G1=list(V=1, nu = 0.002,alpha.mu=0,alpha.V=1000)))
mcd=MCMCglmm(asin(sqrt(fracD))~Location,random=~cn,data=me6,prior=prior,nitt=155000,thin=50,burnin=10000)
summary(mcd)
bigHd=c()
bigHd=(mcd$VCV[,1])/(mcd$VCV[,1]+mcd$VCV[,2])
mean(bigHd) # 0.73 
HPDinterval(bigHd) # 0.62-0.84 


# ----------- simple tests for obvious things

# do inshore clones recover better than non-clones?
me2i=subset(me2,inshore==0)
ci=nrow(subset(me2i,me2i$May2016Condition==3 & me2i$clones!="black")) # 8
co=nrow(subset(me2i,me2i$May2016Condition!=3 & me2i$clones!="black")) # 6
ui=nrow(subset(me2i,me2i$May2016Condition==3 & me2i$clones=="black")) # 47
uo=nrow(subset(me2i,me2i$May2016Condition!=3 & me2i$clones=="black")) # 18
clonesites=matrix(c(ci,ui,co,uo),nrow=2,ncol=2)
fisher.test(clonesites)
# data:  clonesites
# p-value = 0.3384
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
 # 0.1341525 2.0700639
# sample estimates:
# odds ratio 
   # 0.51536 
# ==> Nope.

# do inshore corals recover better than offshore?
ci=nrow(subset(me6,me6$May2016Condition==3 & me6$inshore==1)) # 55
co=nrow(subset(me6,me6$May2016Condition!=3 & me6$inshore==1)) # 24
ui=nrow(subset(me6,me6$May2016Condition==3 & me6$inshore==0)) # 94
uo=nrow(subset(me6,me6$May2016Condition!=3 & me6$inshore==0)) # 0
clonesites=matrix(c(ci,ui,co,uo),nrow=2,ncol=2)
fisher.test(clonesites)
# data:  clonesites
# p-value = 7.268e-10
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
 # 0.0000000 0.1035352
# sample estimates:
# odds ratio 
         # 0 
# ==> Yes.
