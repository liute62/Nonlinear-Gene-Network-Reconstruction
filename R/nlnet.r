# source this code to conduct the non-linear relations network reconstruction.
# array: the data matrix with no missing values
# gene.fdr.cutoff: the local false discovery cutoff in establishing links between genes
# module.min.size: the min number of genes together as a module.
# gene.community.method: it provides three kinds of community detection method: mutilevel, label.propagation and leading.eigenvector.
# normalization: whether to use the normal approximation to for the null hypothesis. 
# If TRUE, normal approximation is used for every feature, AND all covariances are assumed to be zero. If FALSE, generates permutation based null distribution - mean vector and a variance-covariance matrix.
# return the community 
################################################################################

nLRNC<-function(array, gene.fdr.cutoff=0.05,module.min.size=40,gene.community.method="multilevel",normalization="standardize")
{
  
  normrow<-function(array)
  {
    m<-apply(array,1,mean,na.rm=T)
    s<-apply(array,1,sd,na.rm=T)
    array<-(array-m)/s
    return(array)
  }
  
  gene.specific.null<-function(array, B=500)
  {
    null.mat<-matrix(0, nrow=nrow(array), ncol=B)
    l<-ncol(array)
    d.array<-array[,1:(l-1)]
    for(i in 1:B)
    {
      this.order<-sample(l, l, replace=FALSE)
      for(j in 1:(l-1)) d.array[,j]<-abs(array[,this.order[j+1]]-array[,this.order[j]])
      null.mat[,i]<-apply(d.array, 1, sum)
    }
    r<-cbind(apply(null.mat, 1, mean), apply(null.mat, 1, sd))
    return(r)
  }
  
  scol.matrix.order<-function(array,x) # x is the vector, a is the matrix, find ordered distance of rows.of.a|x
  {
    if(is.null(nrow(array)) | nrow(array) == 1)
    {
      array<-as.vector(array)
      array<-array[order(x)]
      d<-array[2:length(array)]-array[1:(length(array)-1)]
      dd<-sum(abs(d),na.rm=T)
    }else{
      array<-a[,order(x)]
      d<-array[,2:ncol(array)]-a[,1:(ncol(array)-1)]
      dd<-apply(abs(d),1,sum,na.rm=T)
    }
    return(dd)
  }
  
  scol.matrix<-function(a, direction=2)  # when direction is 1, scol.matrix[i,j] = SCOL(a[i,], a[j,]), j|i
  {  	
    
    rdmat<-matrix(0, ncol=nrow(a), nrow=nrow(a))
    for(j in 1:nrow(a))
    {
      rdmat[j,]<-scol.matrix.order(a, a[j,])
    }
    
    if(direction == 2)
    {
      rdmat.diff<-rdmat-t(rdmat)
      sel<-which(rdmat.diff > 0)
      rdmat[sel]<-t(rdmat)[sel]
    }
    return(rdmat)
  }
  gene.specific.p<-function(null.distr, new.d)
  {
    for(i in 1:length(new.d))
    {
      new.d[i]<-pnorm(new.d[i], mean=null.distr[i,1], sd=null.distr[i,2], lower.tail=TRUE)
    }
    return(new.d)
  }
  gene.specific.q<-function(new.d)
  {
    for(i in 1:length(new.d))
    {
      new.d[i]<-qnorm(new.d[i], mean=0, sd=1, lower.tail=TRUE)
    }
    return(new.d)
  }
  
  library(fdrtool)
  library(igraph)
  n<-nrow(input)
  m<-ncol(input)
  array<-normrow(input)
  orig.array <-array
  null.distr<-gene.specific.null(array)
  sim.mat<-scol.matrix(array,direction=1)  ## similarity matrix by SCOL, asymmetric, column given row
  d.mat<-sim.mat
  d.tmp.mat<-d.mat
  for(i in 1:nrow(sim.mat)) d.tmp.mat[i,]<-gene.specific.p(null.distr, sim.mat[i,])
  for(i in 1:nrow(sim.mat)) d.mat[i,]<-gene.specific.q(d.tmp.mat[i,])
  diag(d.mat)<-0
  gene.rel.mat<-matrix(0,nrow=n,ncol=n)#the matrix store the relationship between two genes, 1 means having relationship while 0 means no relationship
  for(i in 1:n) {
    sim.vec<-d.mat[i,]
    suppressWarnings(t.locfdr<-fdrtool(sim.vec, statistic="normal", plot=1,color.figure=TRUE,verbose=TRUE,cutoff.method="pct0",pct0=0.75))
    t.row.lfdr<-as.vector(t.locfdr$lfdr)
    for(j in 1:n){
      if(t.row.lfdr[j] < gene.fdr.cutoff){
        gene.rel.mat[i,j]<-1 
      } 
    }
  }
  gene.graph<-graph.adjacency(gene.rel.mat, mode="undirected", weighted=NULL)
 if(gene.community.method=="multilevel"){
    commu<-multilevel.community(gene.graph, weights=NA)
  }else if(gene.community.method=="label.propagation"){
    commu<-label.propagation.community(gene.graph)
  }else if(gene.community.method=="leading.eigenvector.community"){
    commu<-leading.eigenvector.community(gene.graph)
  }
  return(commu);
}

