library("SNFtool")
library("vegan")
source("function_snf.R")
library("reticulate")
library("phyloseq")
library("microbiome")

source_python("sil.py")
use_python("/public/apps/anaconda2/bin/python2")

sps <- load("sps.rda")

rel_sps <- microbiome::transform(sps, "hellinger")
# Hellinger-transformed Bray-Curtis dissimilarity

arc = subset_taxa(rel_sps, k == "k__Archaea")
bac = subset_taxa(rel_sps, k == "k__Bacteria")
fun = subset_taxa(rel_sps, k == "k__Fungi")
vir = subset_taxa(rel_sps, k == "k__Viruses")

a_data<- as.data.frame(t(arc@otu_table))
b_data<- as.data.frame(t(bac@otu_table))
f_data<- as.data.frame(t(fun@otu_table))
v_data<- as.data.frame(t(vir@otu_table))

# a_data<- as.data.frame(t(subset_samples(arc, project != "PRJEB43119")@otu_table))
# b_data<- as.data.frame(t(subset_samples(bac, project != "PRJEB43119")@otu_table))
# f_data<- as.data.frame(t(subset_samples(fun, project != "PRJEB43119")@otu_table))
# v_data<- as.data.frame(t(subset_samples(vir, project != "PRJEB43119")@otu_table))
## PRJEB22863  PRJEB22893  PRJEB43119  PRJEB54704 PRJNA399742 PRJNA541981 PRJNA615114 PRJNA751792 
##         217          25         165          28          39          42          39         335 
## PRJNA762360   SRP115355 
##          94          44 
## a_dsim=vegdist(a_data,method='euclidean',diag=TRUE,upper=TRUE)
## 基于CLR向量的欧氏距离的样品间的比例距离，称为Aitchison距离
a_dsim=vegdist(a_data,method='bray',diag=TRUE,upper=TRUE)
b_dsim=vegdist(b_data,method='bray',diag=TRUE,upper=TRUE)
f_dsim=vegdist(f_data,method='bray',diag=TRUE,upper=TRUE)
v_dsim=vegdist(v_data,method='bray',diag=TRUE,upper=TRUE)

## mantel test
#mantel(a_dsim, b_dsim ,method = "spearman")
## Mantel statistic r: 0.159 
##       Significance: 0.001 
#mantel(a_dsim, f_dsim ,method = "spearman")
## Mantel statistic r: 0.1821 
##       Significance: 0.001
#mantel(a_dsim, v_dsim ,method = "spearman")
## Mantel statistic r: 0.1026 
##       Significance: 0.001 
#mantel(b_dsim, f_dsim ,method = "spearman")
## Mantel statistic r: 0.1642 
##       Significance: 0.001 
#mantel(b_dsim, v_dsim ,method = "spearman")
## Mantel statistic r: 0.1467 
##       Significance: 0.001 
#mantel(f_dsim, v_dsim ,method = "spearman")
## Mantel statistic r: 0.1552 
##       Significance: 0.001 

a_dsim[is.nan(a_dsim)]<-0
b_dsim[is.nan(b_dsim)]<-0
f_dsim[is.nan(f_dsim)]<-0
v_dsim[is.nan(v_dsim)]<-0

W1=(as.matrix(a_dsim)-1)*-1
W2=(as.matrix(b_dsim)-1)*-1
W3=(as.matrix(f_dsim)-1)*-1
W4=(as.matrix(v_dsim)-1)*-1

#Assigning weight values
weight_a=dim(a_data)[2]
weight_b=dim(b_data)[2]
weight_f=dim(f_data)[2]
weight_v=dim(v_data)[2]

weights_snf=c(weight_a,weight_b,weight_f,weight_v)
sil_values=c()
for (i in 2:20){
  W = SNF_weighted_iter(list(W1,W2,W3,W4),i,20,weight = weights_snf)
  z=estimateNumberOfClustersGivenGraph(W)[[1]]
  labels=spectralClustering(W,z)
  sil_values<-c(sil_values,silhouette_score(W,labels))
}
tuned_k<-which.max(sil_values)+1 #since starts from 2
print(paste(tuned_k,sil_values[tuned_k-1],sep = " "))

W = SNF_weighted_iter(list(W1,W2,W3,W4),tuned_k,20,weight = weights_snf)
z=estimateNumberOfClustersGivenGraph(W)[[1]]
labels=spectralClustering(W,z)
print(table(labels))
lab=as.data.frame(labels,row.names = row.names(b_data))

lab$response <- meta(sps)[rownames(lab),]$Response
t <- table(lab)

write.csv(lab,paste("labels",".csv",sep=''))
write.csv(W,paste("matrix",".csv",sep=''))
write.csv(t,paste("matrix",".csv",sep=''))
