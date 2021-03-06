\name{data.gen}
\alias{data.gen}
\title{
  Simulated Data Generation 
}

\description{
 Generating gene matrix as a example of input. 
}

\usage{ 
data.gen(n.genes=100, n.samples=10, n.grps=10, aver.grp.size=10, n.fun.types=3, epsilon=1, n.depend=0)
}

\arguments{
\item{n.genes}{the number of rows of the matrix.
}
 \item{n.samples}{the number of columns of the matrix.
}
 \item{n.grps}{the number of hidden clusters.
}
  \item{aver.grp.size}{averge number of genes in a cluster.
}
  \item{n.fun.types}{number of function types to use.
}
  \item{epsilon}{noise level.
}
  \item{n.depend}{data generation dependence structure. can be 0, 1, 2.
} 
}
\details{
	The data generation scheme is described in detail in IEEE ACM Trans. Comput. Biol. Bioinform. 10(4):1080-85.
}
\value{
 return the data including gene and clustering.
 \item{data}{the gene matrix}
 \item{grps}{the predicted clustering}
}
\author{
Tianwei Yu<tyu8@emory.edu>
}

\examples{
##generating a gene matrix with 100 gene and 100 samples per gene.
data.gen(100,100)
##gene matrix with more than 10 hiden clusters.
data.gen(n.grps=20)
##get the gene matrix from the source of data.
output<-data.gen()
matrix<-output$data
##get the hiden clusters from the source of data.
grps<-output$grp
}

