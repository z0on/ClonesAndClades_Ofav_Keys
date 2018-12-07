library(WGCNA) # for coloring clonal groups
library(sparcl) # for ColorDendrogram

ll=load('manzello_metadata_dec4.RData')

# reading list of bam files = order of samples in IBS matrix
bams=read.table("bams",header=F)[,1]
bams=sub("\\.fastq.*","",bams,perl=T)

# aligning metadata with bams list
row.names(meta)=meta$Run
meta=meta[bams,]

# reading IBS matrix based on SNPs with allele frequency >= 0.05:
ma = as.matrix(read.table("ibs05.ibsMat"))
samples=meta$Tag.Number
dimnames(ma)=list(samples,samples)

# plotting hierarchical clustering tree
hc=hclust(as.dist(ma),"ave")
plot(hc,cex=0.6)

# retaining faveolata-only samples
mafa=ma[meta$faveolata,meta$faveolata]
bamfa=bams[meta$faveolata]
meta=meta[meta$faveolata,]

# re-clustering
hc=hclust(as.dist(mafa),"ave")
plot(hc,cex=0.6)

# "artificial clones" (genotyping replicates)
cl1=c("73.GACT","73.TCAC","73.TGTC")
cl2=c("83.GACT","83.TCAC","83.TGTC")
cl3=c("88.GACT","88.TCAC","88.TGTC")
cl4=c("91.GACT","91.TCAC","91.TGTC")
cl5=c("79.GACT","79.TCAC","79.TGTC")
# color artificial clones red
art.clones=rep("black",nrow(mafa))
art.clones[meta$fileName %in% c(cl1,cl2,cl3,cl4)]="red"
ColorDendrogram(hclust(as.dist(mafa),"ave"), y = art.clones,branchlength=0.05)

# cutoff for defining clones
abline(h=0.15,col="red",lty=3)

# sorting samples into clonal groups; singletons go into the same "color" group (for plotting) but different "cn" groups (for GLM modeling later)
cc=cutree(hc,h=0.15)
meta$cn=cc
tc=table(cc)
singletons=as.numeric(names(tc)[tc==1])
cc[cc %in% singletons]= singletons[1]
library(WGCNA) # just for coloring
clones=labels2colors(as.numeric(as.factor(as.numeric(cc))))
# changing some colors manually for better visibility
clones[clones=="white"]="khaki3"
clones[clones=="black"]="goldenrod"
clones[clones=="turquoise"]="black"
meta$cn.color=clones
ColorDendrogram(hclust(as.dist(mafa),"ave"), y = clones, labels = F,branchlength=0.035)

# selecting one coral per clone
sel=c()
for (i in unique(clones)) {
	if (i=="black") { 
		sel=append(sel,subset(bamfa,clones==i)) 
	} else {
		sel=append(sel,sample(subset(bamfa,clones==i),1))
	}
}
sel=bamfa %in% sel
table(sel)

# writing list of bams with no clonal replicated - rerun the "ibs" ANGSD command with this one instead of 'bams' for PCoA or ADMIXTURE
write.table(paste(bamfa[sel],".fastq.bam",sep=""),quote=F, row.names=F,col.names=F,sep="\t",file="bams_noclones")

save(bamfa,meta,mafa,sel,file='manzello_angsd_dec5.RData')

