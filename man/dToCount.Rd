\name{dToCount.Rd}
\alias{dToCount}
\alias{countToDmixed}
\alias{ff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Functions to switch between count and array representation of an array
}
\description{
dToCount rearranges an array into count vector format. countToDmixed rearranges a count vector representation into an array.
}
\usage{
dToCount(d, nlevels=NULL, startfrom1=FALSE)
countToDmixed(nlevels, count)
ff(...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{
an array (matrix with runs as rows and factors as columns)
}
  \item{nlevels}{
vector of integers (>=2); numbers of factor levels
}
  \item{startfrom1}{
logical; if TRUE, values are from 1 to \code{nlevels}; the default is that values are from 0 to \code{nlevel -1}
}
  \item{count}{
vector of \code{prod(nlevels)} nonnegative integers; \cr
the number of times the respective run of the full factorial in lexicographic order (as produced by \code{ff(nlevels)}) occurs in the array
}
  \item{...}{
vector of integers (>=2) or the integers themselves, separated by commata
}
}
\details{
\code{dToCount} transforms an array into count representation. If all array columns contain all potential factor levels, nlevels does not need to be specified. Otherwise, \code{nlevels} is needed.

\code{countToDmixed} transforms the count representation of an array (counts refer to the rows of \code{ff(nlevels)}) into an array
}
\value{
\code{dToCount} produces a vector of length \code{prod(nlevels)},\cr
\code{countToDmixed} produces a matrix with \code{sum(count)} rows and \code{length(nlevels)} columns,\cr
\code{countToDmixed} produces a matrix with \code{prod(nlevels)} rows and \code{length(nlevels)} columns.\cr
}
\references{
Groemping, U. (2020). DoE.MIParray: an R package for algorithmic creation of orthogonal arrays. \emph{Journal of Open Research Software}, \bold{8: 24}. DOI: https://doi.org/10.5334/jors.286
}

\note{
The size of the full factorial design (produced with \code{ff}) is a limiting factor for using the functionality of this package.
}

\author{
Ulrike Groemping
}


\examples{
d <- ff(c(2,2,4))[1:6,]   ## first six rows of the full factorial only
d
## the count vector must have 2*2*4=16 elements,
## the first six must be 1, the last ten must be zero
dToCount(d)               ## does not produce the desired result,
                          ## because the first column of d
                          ## does not contain both levels
(d_as_count <- dToCount(d, nlevels=c(2,2,4)))
                          ## corresponds to the above way of creating d
dToCount(d, nlevels=c(2,2,5))  ## would correspond to a different reality,
                               ## where the third factor has in fact 5 levels,
                               ## of which only four are in the array
countToDmixed(c(2,2,4), d_as_count)
                          ## creates d from the count representation
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ array }% use one of  RShowDoc("KEYWORDS")
\keyword{ design }% __ONLY ONE__ keyword per line
