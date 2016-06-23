setwd("~/Documents/2015/nlnet/workspace/R Files/R Files/co")

##
## for ute
##
data.src.ute<-read.table("data_ute.txt")
gene.name<-read.table("gene.name.txt")
gene.ute<-cbind(gene.name,data.src.ute)
result.ute<-gene.ute
table<-apply(data.src.ute,1,table)
conn<-1-0.66
input<-vector()
for(i in 1 : length(table)){
  if( table[[i]][1] > conn * ncol(gene.ute)){
    input<-c(input,i)
  }
}
result.ute<-result.ute[-input,]
input.ute<-log10(1+result.ute[,2:ncol(result.ute)])
library(nlnet)
times<-1
for(n in 1:times){
  output.ute<-nlnet(input.ute,min.module.size = 100,plot.method="community")
  output.ute.comm<-output.ute$community
  comm.ute.max<-max(output.ute.comm)
for(i in 1: comm.ute.max) {
  tmp<-vector()
  for(j in 1: length(output.ute.comm) ){
    if(output.ute.comm[j] == i){
      tmp<-c(tmp,as.character(result.ute[j,1]))
    } 
  }
  output.txt.name<-paste(as.character(n),"output_ute_name")
  output.txt.name<-paste(output.txt.name,as.character(i))
  output.txt.name<-paste(output.txt.name,".txt")
  write.table(tmp,output.txt.name)
 }
}
write.table(as.vector(result.ute[,1]),"allname_ute.txt")
##
## now we are using gostats to analysis
##
library(org.Mm.eg.db)
id.ute.all<-unlist(mget(as.vector(result.ute[,1]),org.Mm.egSYMBOL2EG, ifnotfound=NA))
id.ute.all<-id.ute.all[-which(is.na(id.ute.all))]
##please change the n and i range manually
n = 2
for(i in 1: 6) {
  output.txt.name<-paste(as.character(n),"output_ute_name")
  output.txt.name<-paste(output.txt.name,as.character(i))
  output.txt.name<-paste(output.txt.name,".txt")
  name.ute<-read.table(output.txt.name)
  name.ute<-as.vector(name.ute[[1]])
  id.ute<-unlist(mget(name.ute, org.Mm.egSYMBOL2EG, ifnotfound=NA))
  id.ute.re<- -which(is.na(id.ute))
  if(length(id.ute.re) != 0){
    id.ute<-id.ute[-which(is.na(id.ute))]
  }
  sel.entrez<-as.character(id.ute)
  all.entrez<-as.character(id.ute.all)
  params <- new("GOHyperGParams", geneIds=sel.entrez, universeGeneIds=all.entrez, ontology="BP", pvalueCutoff=0.01,conditional=T, testDirection="over", annotation="org.Mm.eg.db")
  over.pres<-hyperGTest(params)
  ov<-summary(over.pres)
  output.txt.name<-paste(as.character(n),"ov_ute")
  output.txt.name<-paste(output.txt.name,as.character(i))
  output.txt.name<-paste(output.txt.name,".txt")
  write.table(ov,output.txt.name)
  print(i)
}
##
## semantic distance
##
library(GOSemSim)
output.ute.graph<-output.ute[2]
output.ute.graph<-output.ute.graph[[1]]
dis.ute.two<-c(as.character(result.ute[7,1]),as.character(result.ute[2425,1]))
dis.ute.two<-c(as.character(result.ute[1000,1]),as.character(result.ute[3609,1]))
id.ute.dis<-unlist(mget(dis.ute.two, org.Mm.egSYMBOL2EG, ifnotfound=NA))
id.ute.dis.out<-geneSim(as.character(id.ute.dis[1]),as.character(id.ute.dis[2]), organism='mouse', ont="BP", measure="Resnik")
##
## transform igraph to networks manually
##
library(network)
output.ute.graph<-output.ute$graph
matrix<-matrix(ncol=3724,nrow=3724)
matrix[1:3724,]<-as.vector(output.ute.graph[1:3724])
network<-as.network.matrix(matrix)







##
##
##
data.src.co<-read.table("data_co.txt")
#gene.name<-read.table("gene.name.txt")
gene.co<-cbind(gene.name,data.src.co)
result.co<-gene.co
table<-apply(data.src.co,1,table)
conn<-1-0.66
input<-vector()
for(i in 1 : length(table)){
  if( table[[i]][1] > conn * ncol(gene.co)){
    input<-c(input,i)
  }
}
result.co<-result.co[-input,]
input.co<-log10(1+result.co[,2:ncol(result.co)])
library(nlnet)
times<-1
output.co<-nlnet(input.co,min.module.size = 100,plot.method = "none")
output.co.comm<-output.co$community
comm.co.max<-max(output.co.comm)
for(i in 1: comm.co.max) {
  tmp<-vector()
  for(j in 1: length(output.co.comm) ){
    if(output.co.comm[j] == i){
      tmp<-c(tmp,as.character(result.co[j,1]))
    } 
  }
  output.txt.name<-paste(as.character(1),"output_co_name")
  output.txt.name<-paste(output.txt.name,as.character(i))
  output.txt.name<-paste(output.txt.name,".txt")
  write.table(tmp,output.txt.name)
}
write.table(as.vector(result.co[,1]),"allname_co.txt")

##
## transform igraph to networks manually
##
library(network)
output.co.graph<-output.co$graph
matrix<-matrix(ncol=2912,nrow=2912)
matrix[1:2912,]<-as.vector(output.co.graph[1:2912])
network<-as.network.matrix(matrix)
##
## now we are using gostats to analysis
##
library(org.Mm.eg.db)
id.co.all<-unlist(mget(as.vector(result.co[,1]),org.Mm.egSYMBOL2EG, ifnotfound=NA))
id.co.all<-id.co.all[-which(is.na(id.co.all))]
##please change the n and i range manually
n<-2
  output.txt.name<-paste(as.character(1),"output_co_name")
  output.txt.name<-paste(output.txt.name,as.character(n))
  output.txt.name<-paste(output.txt.name,".txt")
  name.co<-read.table(output.txt.name)
  name.co<-as.vector(name.co[[1]])
  id.co<-unlist(mget(name.co, org.Mm.egSYMBOL2EG, ifnotfound=NA))
  id.co.re<- -which(is.na(id.co))
  if(length(id.co.re) != 0){
    id.co<-id.co[-which(is.na(id.co))]
  }
  sel.entrez<-as.character(id.co)
  all.entrez<-as.character(id.co.all)
  params <- new("GOHyperGParams", geneIds=sel.entrez, universeGeneIds=all.entrez, ontology="BP", pvalueCutoff=0.01,conditional=T, testDirection="over", annotation="org.Mm.eg.db")
  over.pres<-hyperGTest(params)
  ov<-summary(over.pres)
  output.txt.name<-paste(as.character(n),"ov_co")
  output.txt.name<-paste(output.txt.name,as.character(n))
  output.txt.name<-paste(output.txt.name,".txt")
  write.table(ov,output.txt.name)



#
#For specific ute gene
#
gene.ute.spe.name<-read.table("ute-spe-name.txt")
gene.ute.spe.name<-as.character(gene.ute.spe.name[[1]])
id.ute<-unlist(mget(gene.ute.spe.name, org.Mm.egSYMBOL2EG, ifnotfound=NA))
id.ute.re<- -which(is.na(id.ute))
if(length(id.ute.re) != 0){
  id.ute<-id.ute[-which(is.na(id.ute))]
}
sel.entrez<-as.character(id.ute)
all.entrez<-as.character(id.ute.all)
params <- new("GOHyperGParams", geneIds=sel.entrez, universeGeneIds=all.entrez, ontology="BP", pvalueCutoff=0.01,conditional=T, testDirection="over", annotation="org.Mm.eg.db")
over.pres<-hyperGTest(params)
ov<-summary(over.pres)
output.txt.name<-paste(as.character(n),"ov_ute_spec")
output.txt.name<-paste(output.txt.name,".txt")
write.table(ov,output.txt.name)
#
#Analysis Go similarity between gene and gene-spc
#
ov_id_transfer<-function(spc){
  spc_id<-spc[1]
  spc_id<-spc_id[[1]]
  spc_id<-substr(spc_id,4,11)
  spc_id<-as.numeric(spc_id)
  return(spc_id)
}
ov_id_similarity<-function(ov_nor,ov_spc){
  result<-vector()
  for(i in i : length(ov_spc)){
     if(length(which(ov_nor == ov_spc[i])) != 0){
       result<-append(ov_spc[i])
     }
  }
  return(result)
}

ov_ute_spc<-read.table("1 ov_ute_spec .txt")
ov_ute_spc_id<-ov_id_transfer(ov_ute_spc)


ov_ute<-read.table("1 ov_ute 1 .txt")
ov_ute_id<-ov_id_transfer(ov_ute)
ov_id_similarity(ov_ute_id,ov_ute_spc_id)

ov_ute<-read.table("1 ov_ute 2 .txt")
ov_ute_id<-ov_id_transfer(ov_ute)
ov_id_similarity(ov_ute_id,ov_ute_spc_id)

ov_ute<-read.table("1 ov_ute 3 .txt")
ov_ute_id<-ov_id_transfer(ov_ute)
ov_id_similarity(ov_ute_id,ov_ute_spc_id)

ov_ute<-read.table("1 ov_ute 4 .txt")
ov_ute_id<-ov_id_transfer(ov_ute)
ov_id_similarity(ov_ute_id,ov_ute_spc_id)

ov_ute<-read.table("1 ov_ute 5 .txt")
ov_ute_id<-ov_id_transfer(ov_ute)
ov_id_similarity(ov_ute_id,ov_ute_spc_id)

ov_co_spc<-read.table("1 ov_co_spec .txt")
ov_co_spc_id<-ov_id_transfer(ov_co_spc)

ov_co<-read.table("1 ov_co 1 .txt")
ov_co_id<-ov_id_transfer(ov_co)
ov_id_similarity(ov_co_id,ov_co_spc_id)

#
#For specific co gene
#
gene.co.spe.name<-read.table("co-spe-name.txt")
gene.co.spe.name<-as.character(gene.co.spe.name[[1]])
id.co<-unlist(mget(gene.co.spe.name, org.Mm.egSYMBOL2EG, ifnotfound=NA))
id.co.re<- -which(is.na(id.co))
if(length(id.co.re) != 0){
  id.co<-id.co[-which(is.na(id.co))]
}
sel.entrez<-as.character(id.co)
all.entrez<-as.character(id.co.all)
params <- new("GOHyperGParams", geneIds=sel.entrez, universeGeneIds=all.entrez, ontology="BP", pvalueCutoff=0.01,conditional=T, testDirection="over", annotation="org.Mm.eg.db")
over.pres<-hyperGTest(params)
ov<-summary(over.pres)
output.txt.name<-paste(as.character(n),"ov_co_spec")
output.txt.name<-paste(output.txt.name,".txt")
write.table(ov,output.txt.name)



##
## To see if gene level similarity
##
gene.name.tmp<-read.table("1 output_co_name 1 .txt")
gene.name.tmp<-gene.name.tmp[[1]]
gene.name.tmp<-as.matrix(gene.name.tmp)
gene.co.spe.name
num<-0
tmp<-0
for(i in 1 : length(gene.name.tmp)){
  for(j in 1:length(gene.co.spe.name)){
    if(gene.co.spe.name[j] == gene.name.tmp[i]){
       num<-num+1
       print(i)
       print(j)
       next;
    }
  }
}
num

g<-output.co$graph
plotcord <- data.frame(layout.fruchterman.reingold(g))
colnames(plotcord) = c("X1","X2")

#get edges, which are pairs of node IDs
edgelist <- get.edgelist(g)

#convert to a four column edge data frame with source and destination coordinates
edges <- data.frame(plotcord[edgelist[,1],], plotcord[edgelist[,2],])
colnames(edges) <- c("X1","Y1","X2","Y2")

ggplot() + geom_segment(aes(x=X1, y=Y1, xend = X2, yend = Y2), data=edges, size = 0.5, colour="grey") + geom_point(aes(X1, X2), data=plotcord)

