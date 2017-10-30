args=commandArgs(trailingOnly=TRUE)

require("hicrep")
require("reshape2")
#=============================

testing=function(){
	f1='/oak/stanford/groups/akundaje/oursu/code/forencode/3DChromatin_ReplicateQC/examples/output_test/data/edges/HIC001/HIC001.chr21.gz'
	f2='/oak/stanford/groups/akundaje/oursu/code/forencode/3DChromatin_ReplicateQC/examples/output_test/data/edges/HIC002/HIC002.chr21.gz'
	out='test'
	c1=1
	c2=2
	c3=3
	maxdist=5000000
	resol=40000
	nodefile='/oak/stanford/groups/akundaje/oursu/code/forencode/3DChromatin_ReplicateQC/examples/output_test/data/nodes/nodes.chr21.gz'
	h=5
	m1name='m1'
	m2name='m2'
}

f1=args[1]
f2=args[2]
out=args[3]

c1=1
c2=2
c3=3

maxdist=as.numeric(args[4])
resol=as.numeric(args[5])
nodefile=args[6]
h=as.numeric(args[7])
m1name=args[8]
m2name=args[9]

#=============================

options(scipen=999)

# read in data sets
data1=read.table(f1)
data2=read.table(f2)
d1=data1[,c(c1,c2,c3)]
d2=data2[,c(c1,c2,c3)]

# read in nodes
nodedata=read.table(nodefile)
nodedata=data.frame(nodedata,nodename=as.numeric(as.character(nodedata[,4])))
rownames(nodedata)=nodedata[,'nodename']
nodes=as.character(nodedata[,'nodename'])

# construct matrices
m1=array(0,dim=c(length(nodes),length(nodes)))
rownames(m1)=colnames(m1)=nodes
m2=array(0,dim=c(length(nodes),length(nodes)))
rownames(m2)=colnames(m2)=nodes

d1[,1]=as.character(d1[,1])
d1[,2]=as.character(d1[,2])
d2[,1]=as.character(d2[,1])
d2[,2]=as.character(d2[,2])
d1cast=acast(d1, V1~V2, value.var="V3",fill=0,fun.aggregate=sum)
d2cast=acast(d2, V1~V2, value.var="V3",fill=0,fun.aggregate=sum)
m1[rownames(d1cast),colnames(d1cast)]=m1[rownames(d1cast),colnames(d1cast)]+d1cast
m1[colnames(d1cast),rownames(d1cast)]=m1[colnames(d1cast),rownames(d1cast)]+t(d1cast)

m2[rownames(d2cast),colnames(d2cast)]=m2[rownames(d2cast),colnames(d2cast)]+d2cast
m2[colnames(d2cast),rownames(d2cast)]=m2[colnames(d2cast),rownames(d2cast)]+t(d2cast)

m1_big=data.frame(chr='chromo',n1=as.numeric(as.character(nodedata[,2])),n2=as.numeric(as.character(nodedata[,3])),m1)
m2_big=data.frame(chr='chromo',n1=as.numeric(as.character(nodedata[,2])),n2=as.numeric(as.character(nodedata[,3])),m2)

colnames(m1_big)=gsub('X','',colnames(m1_big))
colnames(m2_big)=gsub('X','',colnames(m2_big))

m1_big[which(is.na(m1_big))]=0
m2_big[which(is.na(m2_big))]=0

#old code
#========
# prepare matrices
#Pre_HiC <- prep(m1_big, m2_big, resol, h, maxdist)
# compute score
#SCC.out = get.scc(Pre_HiC, resol, maxdist)

#new code
#========
#sessionInfo()
SCC.out = get.scc(m1_big[,-c(1:3)], m2_big[,-c(1:3)], resol, h, 0, maxdist)
#print('here')

# write score
scores=data.frame(M1=m1name,M2=m2name,score=SCC.out[['scc']],sd=SCC.out[['std']])
write.table(scores,file=out,quote=FALSE,row.names=FALSE,col.names=FALSE,sep='\t')

#also, plot the correlations and the weights
#png(paste(out,'.png',sep=''))
#par(mfrow=c(2,1))
#plot(0.000001*resol*c(1:length(SCC.out[['corr']])),SCC.out[['corr']],xlab='Genomic distance (Mb)',
#ylab='Pearson correlation',ylim=c(-1,1))
#plot(0.000001*resol*c(1:length(SCC.out[['wei']])),SCC.out[['wei']],xlab='Genomic distance (Mb)',
#ylab='Weights',ylim=c(0,max(SCC.out[['wei']])))
#dev.off()

