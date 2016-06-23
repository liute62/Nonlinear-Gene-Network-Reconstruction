setwd("/Users/liuhaodong1/simu")
source("clusGap.R")
source("DCOL_Kmeans.R")

############# the function to conduct simulation for one scenario
# epsilon: noise level
# n.depend: data generation dependence structure. can be 0, 1, 2
# n.samples: the number of columns of the matrix
# n.grps: the number of hidden clusters
# aver.grp.size: averge number of genes in a cluster
# n.fun.types: number of function types to use
# total.n: number of simulations
# min.size.ratio: to be passed to the dynamic tree cutting

one.simu<-function(cutoff.quantile,module.size, epsilon, n.depend, n.samples=100, n.grps=10, aver.grp.size=100, n.fun.types=4, total.n=20, min.size.ratio)
{
  library(nlnet)
  #print(cutoff.quantile)
  r<-rep(0,4)
  for(N in 1:total.n)
  {
    print(N)
    n.genes<- (1+n.grps) * aver.grp.size
    dataset=data.gen(n.genes, n.samples, n.grps, aver.grp.size, n.fun.types, epsilon, n.depend)
    if (0 %in% unique(dataset$grps))
    {
      ngrp = length(unique(dataset$grps))-1
    }
    else
    {
      ngrp = length(unique(dataset$grps))
    }
    
    print("1.start---mutal info")
    time.1.1<-proc.time()
    rec<-aracne.a(knnmi.all(dataset$data))
    rel<-matrix(ncol = ncol(rec), nrow = nrow(rec))
    #print(table(rec == 0))
    for(i in 1:nrow(rec)){
      for(j in 1:ncol(rec)){
        if(rec[i,j] != 0){
          rel[i,j]<-1
        }else{
          rel[i,j]<-0
        }
      }
    }
    
    gene.graph<-graph.adjacency(rel, mode="undirected", weighted=NULL)
    commu<-multilevel.community(gene.graph, weights=NA)
    print(mem_change(commu))
    time.1.2<-proc.time()
    print(time.1.2-time.1.1)
    print("1.end---mutal info")
    print(mem_change(aracne.a(knnmi.all(dataset$data))))
    
    mutal.info<-commu$membership
    #print(dataset)
    
    print("2.start---data.cluster.2")
    time.2.1<-proc.time()
    results<-data.cluster.2(dataset$data, ngrp, 100)$cluster
    print(mem_change(results))
    time.2.2<-proc.time()
    print(time.2.2-time.2.1)
    print("2.end---data.cluster.2")
    #print(results)
    
    print("3.start---nlnet.mul")
    time.3.1<-proc.time() 
    nlnet.mul<-nlnet(min.module.size=module.size,dataset$data,gene.community.method="multilevel",plot.method="false",conn.proportion=cutoff.quantile)$community
    print(mem_change(nlnet.mul))
    time.3.2<-proc.time()
    print(time.3.2-time.3.1)
    print("3.end---nlnet.mul")
    
    #print(nlnet.mul)
    #print(table(nlnet.mul))
    print("4.start---nlnet.label")
    time.4.1<-proc.time() 
    nlnet.label<-nlnet(min.module.size=module.size,dataset$data,gene.community.method="label.propagation",plot.method="false",conn.proportion=cutoff.quantile)$community
    print(mem_change(nlnet.label))
    time.4.2<-proc.time()
    print(time.4.2-time.4.1)
    print("4.end---nlnet.label")
    
    processable<-"haha"
    print("5.start---GDHC-cuttreeDynamic")
    time.5.1<-proc.time() 
    hc= nlhc(dataset$data)
    suppressWarnings(hgrps<-cutreeDynamic(hc, minClusterSize= aver.grp.size*min.size.ratio))
    print(mem_change(hgrps))
    processable<-try(hgrps2<-cutree(hc, h=-2))
    time.5.2<-proc.time()
    print(mem_change(hgrps2))
    print(time.5.2-time.5.1)
    print("5.end---GDHC-cuttreeDynamic")
    #print(results)
    #table(results)
    #print(hgrps)
    #table(hrgps)
    #print(hgrps2)
    #print(adjustedRandIndex(dataset$grps + 1, hgrps))
    #print(adjustedRandIndex(dataset$grps + 1, hgrps2))
    
    if(substr(processable,1,5)=="Error")
    {
      r<-rbind(r,c(
        adjustedRandIndex(dataset$grps +1 , results+1),
        adjustedRandIndex(dataset$grps + 1, hgrps),
        NA,
        adjustedRandIndex(dataset$grps +1 , nlnet.mul),
        adjustedRandIndex(dataset$grps +1 , nlnet.label),
        adjustedRandIndex(dataset$grps +1 , mutal.info)))
      #adjustedRandIndex(dataset$grps +1 , nlnet.leading)))
    }else{
      r<-rbind(r,c(
        adjustedRandIndex(dataset$grps +1 , results+1),
        adjustedRandIndex(dataset$grps + 1, hgrps),
        adjustedRandIndex(dataset$grps + 1, hgrps2),
        adjustedRandIndex(dataset$grps +1 , nlnet.mul),
        adjustedRandIndex(dataset$grps +1 , nlnet.label),
        adjustedRandIndex(dataset$grps +1 , mutal.info)))
      #adjustedRandIndex(dataset$grps +1 , nlnet.leading)))
    }
    #cat('----r:',r)
  }
  #print("end")
  r<-r[-1,]
  r
}


one.simu.2<-function(cutoff.quantile,module.size, epsilon, n.depend, n.samples=100, n.grps=10, aver.grp.size=100, n.fun.types=4, total.n=20, min.size.ratio)
{
  library(nlnet)
  #print(cutoff.quantile)
  r<-rep(0,4)
  for(N in 1:total.n)
  {
    print(N)
    n.genes<- (1+n.grps) * aver.grp.size
    dataset=data.gen(n.genes, n.samples, n.grps, aver.grp.size, n.fun.types, epsilon, n.depend)
    if (0 %in% unique(dataset$grps))
    {
      ngrp = length(unique(dataset$grps))-1
    }
    else
    {
      ngrp = length(unique(dataset$grps))
    }
    rec<-aracne.a(knnmi.all(dataset$data))
    rel<-matrix(ncol = ncol(rec), nrow = nrow(rec))
    print(table(rec == 0))
    for(i in 1:nrow(rec)){
      for(j in 1:ncol(rec)){
        if(rec[i,j] != 0){
          rel[i,j]<-1
        }else{
          rel[i,j]<-0
        }
      }
    }
    gene.graph<-graph.adjacency(rel, mode="undirected", weighted=NULL)
    commu<-multilevel.community(gene.graph, weights=NA)
    mutal.info<-commu$membership
    processable<-"haha"
    if(substr(processable,1,5)=="Error")
    {
      r<-c(r,adjustedRandIndex(dataset$grps +1 , mutal.info))
      #adjustedRandIndex(dataset$grps +1 , nlnet.leading)))
    }else{
      r<-c(r,adjustedRandIndex(dataset$grps +1 , mutal.info))
      #adjustedRandIndex(dataset$grps +1 , nlnet.leading)))
    }
  }
  r
}

############# actual simulation
library("parmigene")
library(igraph)
all.epsilon<-c(0.2,0.4,0.6,0.8)
all.fdr.quantile<-c(0.035,0.030,0.027,0.015)
n.grps<-c(10, 20)
all.grp.size<-c(100)
#min.size.ratio<-c(0.5, 0.7, 0.9, 1)
min.size.ratio<-c(0.5)
all.n.depend<-c(0,1,2)
all.module.size<-c(10)
#all.epsilon<-c(1.6)
#all.n.depend<-0
#all.grp.size<-100
#all.module.size<-c(10)
#min.size.ratio<-1
#n.grps<-10
#combos<-expand.grid(all.epsilon, all.n.depend,all.grp.size,min.size.ratio, n.grps,all.module.size)
combos<-expand.grid(all.epsilon, all.n.depend,all.grp.size,min.size.ratio, n.grps,all.module.size)
#combos
#foreach(m=nrow(combos):1) %dopar%
#cutoffs<-c(0.1,0.15,0.2,0.25,0.3)
for(m in 1:nrow(combos))
{
  #rec<-one.simu(epsilon=epsilon, n.depend=n.depend, aver.grp.size=grp.size, min.size.ratio=min.size.ratio, n.grps=n.grps, n.fun.types=4, total.n=10)
  cutoff.tmp<-0
  if(combos[m,1] == all.epsilon[1]){
    cutoff.tmp<-all.fdr.quantile[1]
  }else if(combos[m,1] == all.epsilon[2]){
    cutoff.tmp<-all.fdr.quantile[2]
  }else if(combos[m,1] == all.epsilon[3]){
    cutoff.tmp<-all.fdr.quantile[3]
  }else if(combos[m,1] == all.epsilon[4]){
    cutoff.tmp<-all.fdr.quantile[4]
  }
  rec<-one.simu(cutoff.quantile=cutoff.tmp,module.size=combos[m,6],epsilon=combos[m,1], n.depend=combos[m,2], aver.grp.size=combos[m,3], min.size.ratio=combos[m,4], n.grps=combos[m,5], n.fun.types=4, total.n=15)
  #rec<-one.simu.2(cutoff.quantile=cutoff.tmp,module.size=combos[m,6],epsilon=combos[m,1], n.depend=combos[m,2], aver.grp.size=combos[m,3], min.size.ratio=combos[m,4], n.grps=combos[m,5], n.fun.types=4, total.n=15)
  #save(rec, file=paste("simu epsilon",epsilon, "n.depend", n.depend, "grp size",grp.size, "min size ratio",min.size.ratio," n grps",n.grps,".bin"))
  #save(rec, file=paste("simu epsilon", combos[m,1], "n.depend", combos[m,2],"grp size",combos[m,3], "min size ratio", combos[m,4], "n grps", combos[m,5], " mutual.bin"))
}
#}

rel<-matrix(ncol = ncol(rec), nrow = nrow(rec))
for(i in 1:nrow(rec)){
  for(j in 1:ncol(rec)){
    if(rec[i,j]<0.2){
      rel[i,j]<-1
    }else{
      rel[i,j]<-0
    }
  }
}
gene.graph<-graph.adjacency(rel, mode="undirected", weighted=NULL)
commu<-multilevel.community(gene.graph, weights=NA)
mem<-commu$membership
r.mutal<-c(0,0,0,0,0,0,0)
#####for mutal information
for(grp.size in c(100))
{
  #for(n.depend in c(0,1,2))
  #for(n.grps in c(10, 20))
  for(n.depend in c(0,1,2))
  {
    for(n.grps in c(10,20))
    {
      r<-rep(NA,7)
      #for(epsilon in c(0.2, 0.6, 0.8, 1.2, 1.6))
      for(epsilon in c(0.2,0.4,0.6,0.8))
      {
        #for(min.size.ratio in c(0.5, 0.7, 0.9, 1))
        for(min.size.ratio in c(0.5))
        {
          #this.name<-paste("linear simu epsilon", epsilon, "n.depend", n.depend,"grp size",grp.size, "min size ratio", min.size.ratio, "n grps", n.grps, ".bin")
          this.name.m<-paste("simu epsilon",epsilon, "n.depend", n.depend, "grp size",grp.size, "min size ratio",min.size.ratio,"n grps", n.grps," mutual.bin")
          load(this.name.m)
          #load("simu epsilon 1.2 n.depend 2 grp size 10 min size ratio 1 n grps 20.bin")
          this<-mean(rec[which(rec[1:length(rec)] != 0)])
        }
        r<-rbind(r, c(epsilon,this)) 
      }  
      print(r)
      r<-r[-1,]
      r.mutal<-rbind(r.mutal,r)
    }
  }
}
r.mutal<-r.mutal[-1,]
############# plotting results
par(mfcol=c(2,3))
num<-1
for(grp.size in c(100))
{
  #for(n.depend in c(0,1,2))
  #for(n.grps in c(10, 20))
  for(n.depend in c(0,1,2))
  {
    for(n.grps in c(10,20))
    {
      r<-rep(NA,7)
      #for(epsilon in c(0.2, 0.6, 0.8, 1.2, 1.6))
      for(epsilon in c(0.2,0.4,0.6,0.8))
      {
        #for(min.size.ratio in c(0.5, 0.7, 0.9, 1))
        for(min.size.ratio in c(0.5))
        {
          #this.name<-paste("linear simu epsilon", epsilon, "n.depend", n.depend,"grp size",grp.size, "min size ratio", min.size.ratio, "n grps", n.grps, ".bin")
          this.name<-paste("simu epsilon",epsilon, "n.depend", n.depend, "grp size",grp.size, "min size ratio",min.size.ratio,"n grps", n.grps,".bin")
          #this.name.m<-paste("simu epsilon",epsilon, "n.depend", n.depend, "grp size",grp.size, "min size ratio",min.size.ratio,"n grps", n.grps," mutual.bin")
          load(this.name)
          #load(this.name.m)
          this.name.m
          #load("simu epsilon 1.2 n.depend 2 grp size 10 min size ratio 1 n grps 20.bin")
          this<-apply(rec[,1:5],2,mean)
          size<-nrow(rec)
          sum<-0
          for(i in 1:nrow(rec)){
            #get rid of the noise node.
            if(rec[i,6] < 0.1){
              size<-size-1
            }else{
              sum<-sum+rec[i,6]
            }
          }
          this[6]<-sum/size
        }
        r<-rbind(r, c(epsilon, this,r.mutal[num,2])) 
        num<-num+1
      }  
      r<-r[-1,]
      print(r)
      plot(r[,1], r[,2], xlab="
           ", ylab="adjusted Rand Index", type="l", col="red",ylim=c(0,1),lwd=2, main=paste("# clusters", n.grps,", n.depend", n.depend))
      lines(r[,1], r[,3],col="yellow",lwd=2, lty=1)
      #lines(r[,1], r[,4],col="blue",lwd=2, lty=1)
      lines(r[,1], r[,5],col="blue",lwd=2, lty=1)
      lines(r[,1], r[,6],col="green",lwd=2, lty=1)
      #lines(r[,1], r[,7],col="purple",lwd=2, lty=1)
      lines(r[,1], r[,8],col="purple",lwd=2, lty=1)
      #lines(r[,1], r[,8],col="cyan",lwd=2, lty=1)
      #lines(r[,1], r[,9],col="yellow",lwd=2,lty=1)
    }
  }
}


