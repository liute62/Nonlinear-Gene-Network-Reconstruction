\name{nLRNC-package}
\alias{nLRNC-package}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Non-Linear gene Relationship Network re-Construction}
\description{To find the non-linear relation between genes and then detect the community}
\usage{
nLRNC(array, gene.fdr.cutoff = 0.05, module.min.size = 40, gene.community.method = "multilevel", normalization = "standardize")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{array}{
%%     ~~Describe \code{array} here~~
}
  \item{gene.fdr.cutoff}{
%%     ~~Describe \code{gene.fdr.cutoff} here~~
}
  \item{module.min.size}{
%%     ~~Describe \code{module.min.size} here~~
}
  \item{gene.community.method}{
%%     ~~Describe \code{gene.community.method} here~~
}
  \item{normalization}{
%%     ~~Describe \code{normalization} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (array, gene.fdr.cutoff = 0.05, module.min.size = 40, 
    gene.community.method = "multilevel", normalization = "standardize") 
{
    normrow <- function(array) {
        m <- apply(array, 1, mean, na.rm = T)
        s <- apply(array, 1, sd, na.rm = T)
        array <- (array - m)/s
        return(array)
    }
    gene.specific.null <- function(array, B = 500) {
        null.mat <- matrix(0, nrow = nrow(array), ncol = B)
        l <- ncol(array)
        d.array <- array[, 1:(l - 1)]
        for (i in 1:B) {
            this.order <- sample(l, l, replace = FALSE)
            for (j in 1:(l - 1)) d.array[, j] <- abs(array[, 
                this.order[j + 1]] - array[, this.order[j]])
            null.mat[, i] <- apply(d.array, 1, sum)
        }
        r <- cbind(apply(null.mat, 1, mean), apply(null.mat, 
            1, sd))
        return(r)
    }
    scol.matrix.order <- function(array, x) {
        if (is.null(nrow(array)) | nrow(array) == 1) {
            array <- as.vector(array)
            array <- array[order(x)]
            d <- array[2:length(array)] - array[1:(length(array) - 
                1)]
            dd <- sum(abs(d), na.rm = T)
        }
        else {
            array <- a[, order(x)]
            d <- array[, 2:ncol(array)] - a[, 1:(ncol(array) - 
                1)]
            dd <- apply(abs(d), 1, sum, na.rm = T)
        }
        return(dd)
    }
    scol.matrix <- function(a, direction = 2) {
        rdmat <- matrix(0, ncol = nrow(a), nrow = nrow(a))
        for (j in 1:nrow(a)) {
            rdmat[j, ] <- scol.matrix.order(a, a[j, ])
        }
        if (direction == 2) {
            rdmat.diff <- rdmat - t(rdmat)
            sel <- which(rdmat.diff > 0)
            rdmat[sel] <- t(rdmat)[sel]
        }
        return(rdmat)
    }
    gene.specific.p <- function(null.distr, new.d) {
        for (i in 1:length(new.d)) {
            new.d[i] <- pnorm(new.d[i], mean = null.distr[i, 
                1], sd = null.distr[i, 2], lower.tail = TRUE)
        }
        return(new.d)
    }
    gene.specific.q <- function(new.d) {
        for (i in 1:length(new.d)) {
            new.d[i] <- qnorm(new.d[i], mean = 0, sd = 1, lower.tail = TRUE)
        }
        return(new.d)
    }
    library(fdrtool)
    library(igraph)
    n <- nrow(input)
    m <- ncol(input)
    array <- normrow(input)
    orig.array <- array
    null.distr <- gene.specific.null(array)
    sim.mat <- scol.matrix(array, direction = 1)
    d.mat <- sim.mat
    d.tmp.mat <- d.mat
    for (i in 1:nrow(sim.mat)) d.tmp.mat[i, ] <- gene.specific.p(null.distr, 
        sim.mat[i, ])
    for (i in 1:nrow(sim.mat)) d.mat[i, ] <- gene.specific.q(d.tmp.mat[i, 
        ])
    diag(d.mat) <- 0
    gene.rel.mat <- matrix(0, nrow = n, ncol = n)
    for (i in 1:n) {
        sim.vec <- d.mat[i, ]
        suppressWarnings(t.locfdr <- fdrtool(sim.vec, statistic = "normal", 
            plot = 1, color.figure = TRUE, verbose = TRUE, cutoff.method = "pct0", 
            pct0 = 0.75))
        t.row.lfdr <- as.vector(t.locfdr$lfdr)
        for (j in 1:n) {
            if (t.row.lfdr[j] < gene.fdr.cutoff) {
                gene.rel.mat[i, j] <- 1
            }
        }
    }
    gene.graph <- graph.adjacency(gene.rel.mat, mode = "undirected", 
        weighted = NULL)
    if (gene.community.method == "multilevel") {
        commu <- multilevel.community(gene.graph, weights = NA)
    }
    else if (gene.community.method == "label.propagation") {
        commu <- label.propagation.community(gene.graph)
    }
    else if (gene.community.method == "leading.eigenvector.community") {
        commu <- leading.eigenvector.community(gene.graph)
    }
    return(commu)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
