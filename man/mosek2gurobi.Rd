\name{mosek2gurobi}
\alias{mosek2gurobi}
\alias{gurobi2mosek}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Functions to recast quadratically constrained MIP in different format
}
\description{
The functions recast a Mosek model into Gurobi format and vice versa,
for use with objects of class qco from package DoE.MIParray.
}
\usage{
mosek2gurobi(qco, ...)
gurobi2mosek(qco, ...)
}
\arguments{
  \item{qco}{
    a mixed integer optimization problem of class qco, as generated
    from package \pkg{DoE.MIParray}
}
  \item{\dots}{
   not used so far
}
}
\details{
   The functions only treat the special cases treated in package DoE.MIParray,
   (minimization problem, all linear constraints are equality constraints,
   cones are created in the form done in that package).
}
\value{
an object of S3 class \code{qco}
}
\author{
  Ulrike Groemping
  }

\seealso{
See also as \code{\link{mosek_MIParray}}, \code{\link{gurobi_MIParray}}.
}
