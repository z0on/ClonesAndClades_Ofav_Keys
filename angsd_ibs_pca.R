bams=read.table("bams_noclones")[,1] # list of bam files

#--------------------
# loading individual to population correspondences
i2p=read.table("inds2pops",sep="\t") # 2-column tab-delimited table of individual assignments to populations; must be in the same order as samples in the bam list or vcf file.
row.names(i2p)=i2p[,1]
site=as.character(i2p[,2])

# settign up colors for plotting
palette(rainbow(length(unique(site))))
colors=as.numeric(as.factor(site))
colpops=as.numeric(as.factor(sort(unique(site))))

#-------------
# clustering / PCoA based on identity by state (IBS) based on single read resampling
# (for low and/or uneven coverage)

ma = as.matrix(read.table("ofav_nc.ibsMat"))
hc=hclust(as.dist(ma),"ave")
plot(hc,cex=0.5)  # this shows how similar clones are

# performing PCoA and CAP
library(vegan)
conds=data.frame(cbind(site))
pp0=capscale(ma~1)
pp=capscale(ma~site,conds)

# significance of by-site divergence
adonis(ma~site,conds)

# eigenvectors
plot(pp0$CA$eig) 

axes2plot=c(1,2)  
quartz()
library(adegenet) # for transp()
cmd=pp0  # change to cmd=pp to see constrained ordination (data projection to maximize by-site separation)
plot(cmd,choices=axes2plot,display="sites",type="n") # choices - axes to display
points(cmd,choices=axes2plot,pch=19,col=transp(colors,alpha=0.7))
#ordihull(cmd,choices= axes2plot,groups= conds$grp,draw="polygon",col=1+as.numeric(unique(as.factor(conds$grp))),label=T)
ordispider(cmd,choices= axes2plot,groups=conds$site,col="grey80")
ordiellipse(cmd,choices= axes2plot,groups= conds$site,draw="polygon",col=colpops,label=T)

# unscaled, to identify outliers
plot(cmd$CA$u[,axes2plot],pch=19,col=colors)
ordispider(cmd$CA$u[,axes2plot],groups=conds$site,col="grey80")
ordiellipse(cmd$CA$u[,axes2plot],groups= conds$site,draw="polygon",col=colpops,label=T)
identify(cmd$CA$u[,axes2plot],labels=colnames(ma),n=3,cex=0.7)

