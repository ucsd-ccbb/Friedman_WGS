##########################################
#import SNV table of mutations in patients
##########################################
#this is filtered for FLAGS, and has only deleterious SNPs (in both PolyPhen and SIFT):
X<-read.table("menieres.531.gnomadv3.filtered.noFLAGS.toleratedSNPout.INDELs.pvalue.lfdr.0.1.tsv",sep="\t",header=TRUE,stringsAsFactors=FALSE,quote="\"")
#dim(X)
#[1] 74359    140

sampleID<-X$Tumor_Sample_Barcode
nmut<-as.data.frame(table(sampleID))
#median number of deleterious mutations per patient
median(nmut$Freq)
#[1] 88

samples<-sort(unique(sampleID))
ns<-length(samples)
#511 samples


##########################################
#import STRING 11.5 gene-gene network
##########################################
library(igraph)
library(ndexr)
library(mygene)
library(Matrix)

#login to the NDEx server. Use your own username and password.
ndexcon<-ndex_connect("username", "password")
rcx<-ndex_get_network(ndexcon, "11937658-62ea-11ee-aa50-005056ae23aa")
#this ^^ takes a few minutes
#convert network to igraph format:
g2<-toIgraph(rcx)

vsym<-V(g2)$nodeName
nv<-length(V(g2))
#as edge weight is character, transform into numeric
w<-as.numeric(E(g2)$weight)
g2<-set_edge_attr(g2,name="w",value=w)

##########################################
#calculate observed frequency of mutations per gene in cohort
##########################################
#simple SNV list first
simple<-list() #this list contains ENSG
for (s in samples) {
   sym<-unique(X$SYMBOL[X$Tumor_Sample_Barcode==s])
   simple[[s]]<-intersect(sym,vsym)
}
m<-sapply(simple,length)

Y0<-matrix(0,nrow=nv,ncol=ns) #input matrix
for (j in 1:ns) {
   Y0[,j]<- vsym %in% simple[[j]]
}
frequency<-rowSums(Y0) #per gene
names(frequency)<-vsym


##########################################
#calculate the Wprime matrix for network propagation
##########################################
A<-as_adjacency_matrix(g2,attr="w")
dA<-rowSums(A)
a<-0.5 #alpha in network propagation. Higher a means broader propagation.
Wprime<-t(A/dA) #colSums will be 1


##########################################
#network propagation with the observed SNV's
##########################################
Y0<-frequency
Y<-Y0
for (i in 1:20) Y<-a*Wprime %*% Y + (1-a)*Y0
R0<-Y[,1]


##########################################
#generate nmc random samples of the null model
##########################################
nmc<-10000
set.seed(2357)
Ymc<-matrix(0,nrow=nv,ncol=nmc) #input matrix
for (mc in 1:nmc) {
   if ((mc %% 100)==0) cat(mc,"\n")
   Y<-matrix(0,nrow=nv,ncol=ns) #input matrix
   for (j in 1:ns) {
      i<-sample(nv,m[j])
      Y[i,j]<-1
   }
   Ymc[,mc]<-rowSums(Y)
}
Rmc<-matrix(0,nrow=nv,ncol=nmc)
for (mc in 1:nmc) {
   if ((mc %% 100)==0) cat(mc,"\n")
   Y0<-Ymc[,mc]
   Y<-Y0
   for (i in 1:20) Y<-a*Wprime %*% Y + (1-a)*Y0 #works
   Rmc[,mc]<-Y[,1]
}
 
##########################################
#calculate z from log10-transformed values
##########################################
lR0<-log10(R0)
lRmc<-log10(Rmc)
z<-(lR0-rowMeans(lRmc))/apply(lRmc,1,sd)

##########################################
#use B. Efron's empirical Bayes local fdr method
##########################################
library(locfdr)

efron<-locfdr(z,df=29)
lfdr<-efron$fdr
#because genes with significantly fewer observed mutations than expected are not interesting to us:
lfdr[z<0]<-1 

#save result for GSEA
res<-data.frame(vsym,z,lfdr,frequency)

