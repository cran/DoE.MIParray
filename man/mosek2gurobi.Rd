\name{mosek2gurobi}
\alias{mosek2gurobi}
\alias{gurobi2mosek}
\alias{qco}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Functions to recast quadratically constrained MIP in different format, and class qco
}
\description{
The functions recast a Mosek model into Gurobi format and vice versa,
for use with objects of class qco from package DoE.MIParray. The class is also 
documented here.
}
\usage{
mosek2gurobi(qco, ...)
gurobi2mosek(qco, ...)
}
\arguments{
  \item{qco}{
    a mixed integer optimization problem of class \code{qco}, as generated
    from package \pkg{DoE.MIParray}
}
  \item{\dots}{
   not used so far
}
}
\details{
   The functions treat the special \code{qco} objects created by package 
   \pkg{DoE.MIParray}: these are minimization problems with linear equality 
   constraints and possibly conic quadratic constraints, as suitable for the 
   problems treated in \pkg{DoE.MIParray}.
   
   Class \code{qco} objects on their own only occur as interim results of the 
   optimization functions \code{mosek_MIParray}, 
   \code{mosek_MIPcontinue}, 
   \code{gurobi_MIParray} or \code{gurobi_MIPcontinue}. Where it might be useful, 
   the class \code{link[DoE.base]{oa}} output objects of the optimization functions 
   contain an attribute \code{MIPinfo} of class \code{qco}. For reducing the size of 
   an object that is not going to be used for further improvement, 
   the following command can be run for extracting the 
   the useful information content from the \code{qco} object and replacing the large
   \code{MIPinfo} attribute with this much smaller object: \cr
   \code{attr(obj, "MIPinfo") <- attr(obj, "MIPinfo")$info}\cr
   Make sure to only run this command if \code{MIPinfo} attribute is indeed of class \code{qco} 
   and further optimization is not intended.
}
\value{
an object of S3 class \code{qco} (see Details section)
}
\author{
  Ulrike Groemping
  }

\seealso{
See also as \code{\link{mosek_MIParray}}, \code{\link{gurobi_MIParray}}.
}
