\name{nlnet}
\alias{nlnet}
\title{
  Non-Linear Network reconstruction
}
\description{
  Non-Linear Network reconstruction method
}
\usage{
nlnet(input, gene.fdr.cutoff=0.05, gene.fdr.plot=FALSE, module.min.size=40, gene.community.method="multilevel", use.normal.approx=FALSE, normalization="standardize", plot.method="communitygraph")
}
\arguments{
  \item{input}{the data matrix with no missing values.
}
  \item{gene.fdr.cutoff}{
   the local false discovery cutoff in establishing links between genes.
}
 \item{gene.fdr.plot}{
  whether plot a figure with estimated densities, distribution functions, and (local) false discovery rates.
}
  \item{module.min.size}{
  the min number of genes together as a module.
}
  \item{gene.community.method}{the method for community detection.
}
  \item{use.normal.approx}{whether to use the normal approximation for the null hypothesis. 
}
  \item{normalization}{the normalization method for the array. 
}
  \item{plot.method}{
  the method for graph and community ploting
}
}
\details{
gene.community.method: It provides three kinds of community detection method:
  "mutilevel", "label.propagation" and "leading.eigenvector".
  
use.normal.approx: If TRUE, normal approximation is used for every feature, AND all covariances are assumed to be zero. If FALSE, generates permutation based null distribution - mean vector and a variance-covariance matrix.

normalization: There are three choices:
"standardize" means removing the mean of each row and make the standard deviation one; 
"normal_score" means normal score transformation; 
"none" means do nothing. In that case we still assume some normalization has been done by the user such that each row has approximately mean 0 and sd 1. 

plot.method: It provides three kinds of ploting method:
  "none" means ploting no graph,
  "communitygraph" means ploting community with graph, 
  "graph" means ploting graph, 
  "membership" means ploting membership of the community 
}
\value{
  it returns a community membership of the graph.
 \item{membership }{Numeric vector, one value for each vertex, the membership vector of the community structure.}
 \item{memberships}{The memberships of the community structure}
 \item{modularity}{Numeric scalar or vector, the modularity value of the community structure. It can also be NULL, if the modularity of the (best) split is not available.}
 \item{vcount}{The number of the item in graph}
 \item{algorithm}{The algorithm name for community detection}
}
\author{
  Haodong Liu <liute62@sina.cn>
}


\seealso{
  \code{\link{data.gen}}
}
\examples{
 
 ## generating the data matrix & hiden clusters as a sample
 input<-data.gen()
 ##Now input includes data matrix and hiden clusters, so get the matrix as input.
 input<-input$data
 
 ##change the ploting method
 nlnet(input,plot.method="graph")
 ## see fdrtool image
 nlnet(input,gene.fdr.plot=TRUE)
 
 ## use different community detection method
 nlnet(input,gene.community.method="label.propagation")
 
 ## change the localfdr cutoff or adjust the modularity size
 nlnet(input,gene.fdr.cutoff=0.10,module.min.size=10)
 
 ## get the result and see it values
 comm<-nlnet(input)
 mem<-comm$membership ##membership of community
 modul<-comm$modularity ##modularity of community
}
