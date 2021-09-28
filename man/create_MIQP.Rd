\name{create_MIQP.Rd}
\alias{write_MPSILP}
\alias{create_MIQP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Internal functions that support the export of (mixed) integer quadratic or linear problems in MPS format
}
\description{
create_MIQP creates an MPS format character object,
write_MPSILP saves a single integer linear problem in MPS format (used in the user-visible write_MSPILPlist).
The functions work, if neither Mosek nor Gurobi is available.
}
\usage{
create_MIQP(nruns, nlevels, resolution=3, distinct=TRUE,
              start=NULL, forced=NULL, forMosek=FALSE)
write_MPSILP(file, qco=NULL, qcoinfo=NULL,
              cvec=NULL, bvec=NULL, Amat=NULL, bsense=NULL, boundmat=NULL,
              name = "EnforceResolution", commentline="* no comment")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nruns}{
positive integer; number of runs
}
  \item{nlevels}{
vector of integers (>=2); numbers of factor levels
}
  \item{resolution}{
positive integer; the minimum resolution requested
}
  \item{distinct}{
logical; if TRUE (default), restricts counting vector to 0/1 entries, which means that the resulting array is requested to have distinct rows; otherwise, duplicate rows are permitted, i.e. the counting vector can have arbitrary non-negative integers. Designs with distinct runs are usually better; in addition, binary variables are easier to handle by the optimization algorithm. Nevertheless, there are occasions where a better array is found faster with option \code{distinct=FALSE}, even if it has distinct rows.
}
  \item{start}{
a starting value for the algorithm: can be a array matrix with entries 1 to number of levels for each column, or a counting vector for the full factorial in lexicographic order; if specified, \code{start} must specify an array with the appropriate number of rows and columns, the requested resolution and, if \code{distinct = TRUE}, also contain distinct rows (matrix) or 0/1 elements only. \code{start} cannot be combined with \code{search.orders=TRUE}.
}
  \item{forced}{
for \code{resolution > 1} only;\cr
runs to force into the solution design; can be given as an array matrix with the appropriate number of columns and less than \code{nruns} rows or a counting vector for the full factorial in lexicographic order with sum smaller than \code{nruns}; if \code{distinct=TRUE}, \code{forced} must have distinct rows (matrix) or 0/1 elements only.
}
  \item{forMosek}{
logical}
  \item{file}{
character string that provides the path to the save file
}
  \item{qco}{
an object of class \code{qco} that can be used as input for function \code{Rmosek::mosek}
}
  \item{qcoinfo}{
contains information on the \code{qco} object (keeping that information in \code{qco} itself hinders execution of \code{Rmosek::mosek})
}
  \item{cvec}{
coefficients of linear objective function
}
  \item{bvec}{
RHS of the constraints
}
  \item{Amat}{
matrix of coefficients for the LHS of the constraints
}
  \item{bsense}{
directions of constraints (character vector with elements \code{"E"} for equal, \code{"G"} for greater, \code{"L"} for less or \code{"N"} for not restricted.
}
  \item{boundmat}{
matrix of lower (first row) and upper (second row) bounds for the variables
}
  \item{name}{
name to be written into the name field of the mps file
}
  \item{commentline}{
comment to be written directly underneath the name (comments have to start with a *)
}
  \item{prefix}{
file name prefix to be supplemented with a number and the .mps suffix later. In addition, a further file with a different suffix may be created (see details).
}
  \item{qcolist}{
list of ILP objects in Mosek notation (as e.g. produced by function \code{create_ILP_list})
}
  \item{toc}{
logical (default TRUE) that indicates whether a table of content should be produced (txt file)
  }
}
\details{
The functions do \emph{not} do any problem solving, they serve the sole purpose of exporting problems so that they can be addressed by solvers other than Mosek or Gurobi. The target solver must be able to process the MPS format.

\code{write_MPSILP} writes a single problem to an MPS file. The function is used by function \code{\link{write_MPSILPlist}}.

\code{create_MIQP} creates a quadratic problem in MPS format, using a general quadratic representation (section \code{Q}) for external solvers (default) or Mosek format with conic quadratic constraints (\code{forMosek=TRUE}); the function is used by function \code{\link{write_MPSMIQP}}.

The start array for \code{create_MIQP} can be made available from a previous linear optimization (though this is not enforced).
If an admissible start array is provided, \code{write_MPSMIQP} initially prints the GWLP of that start array.
Otherwise, it warns of inadmissibility. The start value cannot be stored in the MPS file. Instead, for a non-NULL start array,
a separate file (suffix \code{.start}) is created, and users have to work out how they can make their solver use that start solution.
Note that the availability of a start array can improve the ability of a solver to find an optimum solution. However, this is not always the case,
there are also instances for which a better solution is found without providing a start solution.

Note that it can take a long time to write the problems, if the problem has many variables (the number of variables is \code{prod(nlevels)}).
}
\value{
Function \code{create_MIPQ} creates a list object of class \code{qco}, see also section "Details".

Function \code{write_MPSILP} does not return a value.
}
\references{
Groemping, U. and Fontana, R. (2019). An Algorithm for Generating Good Mixed Level Factorial Designs. \emph{Computational Statistics & Data Analysis} \bold{137}, 101-114. doi:10.1016/j.csda.2019.01.020.

Mosek ApS (2017a). MOSEK version w.x.y.z documentation. Accessible at: \url{https://www.mosek.com/documentation/}. This package has been developed using version 8.1.0.23 (accessed August 29 2017).

Mosek ApS (2017b). MOSEK Rmosek Package 8.1.y.z. \url{https://docs.mosek.com/8.1/rmosek/index.html}.
\emph{!!! In normal R speak, this is the documentation of the Rmosek package version 8.0.69 (or whatever comes next), when applied on top of the Mosek version 8.1.y.z (this package has been devoloped with Mosek version 8.1.0.23 and will likely not work for Mosek versions before 8.1). !!! (accessed August 29 2017)}

}

\note{
The functions are not meant for situations, for which a full factorial design would be huge. Even though the functions do not solve anything, MPS files will be very large and writing them will be quite slow for such cases.
}

\author{
Ulrike Groemping
}

\seealso{
See also \code{\link{write_MPSILPlist}} and \code{\link{write_MPSMIQP}}.
}
\examples{
\dontrun{
## can also be run with gurobi_MIParray instead of mosek_MIParray
## there are of course better ways to obtain good arrays for these parameters
## (e.g. function FrF2 from package FrF2)
problemlist <- create_ILPlist(16, c(2,2,4,8), resolution=2, kmax=4)
feld
names(attributes(feld))
attr(feld, "MIPinfo")$info

## using a start value
start <- DoE.base::L16.2.8.8.1[,1:5]
feld <- mosek_MIParray(16, rep(2,5), resolution=4, start=start)

## counting vector representation of the start value could also be used
DoE.MIParray:::dToCount(start-1)
   ## "-1", because the function requires values starting with 0
   ## 32 elements for the full factorial in lexicographic order, 16 ones for the runs

## extending an existing array
force <- matrix(as.numeric(as.matrix(DoE.base::undesign(DoE.base::oa.design(L8.2.7)))), nrow=8)
feld <- mosek_MIParray(16, rep(2,7), resolution=3, kmax=4, forced=force)
attr(feld, "MIPinfo")$info
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ array }% use one of  RShowDoc("KEYWORDS")
\keyword{ design }% __ONLY ONE__ keyword per line
\keyword{ internal }% __ONLY ONE__ keyword per line
