\name{write_MPSILPlist.Rd}
\alias{write_MPSILPlist}
\alias{write_MPSMIQP}
\alias{create_ILPlist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Functions to create and write lists of (mixed) integer quadratic or linear problems related to orthogonal arrays
}
\description{
create_ILPlist creates a list of problems in Mosek formatting.
write_MPSILPlist saves a list of integer linear problems as separate MPS files that are accompanied by a table of content txt file. write_MPSMIQP creates and writes a single quadratic mixed integer problem in MPS format.
All functions work, even if neither Mosek nor Gurobi is available.
}
\usage{
create_ILPlist(nruns, nlevels, resolution=3, distinct=TRUE,
              search.orders=TRUE, start=NULL, forced=NULL, orders=NULL)
write_MPSILPlist(prefix, qcolist, toc=TRUE)
write_MPSMIQP(prefix, nruns, nlevels, resolution=3, distinct=TRUE, start=NULL,
              forced=NULL, name="ImproveAR", commentline="* quadratic problem")
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
  \item{search.orders}{
logical (default TRUE); determines whether a list of arrays for different level orderings is produced or a single array is output only
}
  \item{start}{
a starting value for the algorithm: can be a matrix with entries 1 to number of levels for each column, or a counting vector for the full factorial in lexicographic order; if specified, \code{start} must specify an array with the appropriate number of rows and columns, the requested resolution and, if \code{distinct = TRUE}, also contain distinct rows (matrix) or 0/1 elements only. \code{start} cannot be combined with \code{search.orders=TRUE}.
}
  \item{forced}{
for \code{resolution > 1} only;\cr
runs to force into the solution design; can be given as an array matrix with the appropriate number of columns and less than \code{nruns} rows or a counting vector for the full factorial in lexicographic order with sum smaller than \code{nruns}; if \code{distinct=TRUE}, \code{forced} must have distinct rows (matrix) or 0/1 elements only.
}
  \item{orders}{
a list of level orderings to be considered; if \code{orders} is not specified but \code{search.orders=TRUE}, all distinct level orderings
}
  \item{prefix}{
file name prefix to be supplemented with a number and the .mps suffix later. In addition, a further file with a different suffix may be created (see details).
}
  \item{qcolist}{
list of ILP objects in Mosek notation (as e.g. produced by function \code{create_ILPlist})
}
  \item{toc}{
logical (default TRUE) that indicates whether a table of content should be produced (txt file)
  }
  \item{name}{
name to be written into the name field of the mps file
}
  \item{commentline}{
comment to be written directly underneath the name (comments have to start with a *)
}
}
\details{
The functions do \emph{not} do any problem solving, they serve the sole purpose of exporting problems so that they can be addressed by solvers other than Mosek or Gurobi. The target solver must be able to process the MPS format.

\code{create_ILPlist} creates a list with one or more integer linear problems for creating designs of at least the specified resolution. The problems are in Mosek formatting.

\code{write_MPSILPlist} writes a list produced by \code{create_ILPlist} to one or several MPS files, one file per list element. If not suppressed, a toc file is also created (suffix \code{_toc.txt}), which contains the number of runs, the target resolution and the elements of \code{nlevels} for each MPS file.

\code{write_MPSMIQP} creates and writes a problem for improving the number of shortest words, in a form that corresponds to general MPS format, without extensions by Mosek or Gurobi. The problem has a quadratic objective like in Groemping and Fontana (2019) Eq.(2Q) (instead of Eq.(2L), which is solved by Gurobi and Mosek).

The example section demonstrates on a small example how these functions can be used.

For \code{write_MPSMIQP}, a start array can be provided from a previous linear optimization (this is not enforced).
If an admissible start array is provided, \code{write_MPSMIQP} initially prints the GWLP of that start array.
Otherwise, it warns of inadmissibility. The start value cannot be stored in the MPS file. Instead, for a non-NULL start array,
a separate file (suffix \code{.start}) is created, and users have to work out how they can make their solver use that start solution (which is stored in the form of a counting vector, see example section).
Note that the availability of a start array can improve the ability of a solver to find an optimum solution. However, this is not always the case,
there are also instances for which a better solution is found without providing a start solution.

For \code{write_MPSMIQP}, a lower bound for the objective value
can substantially improve the run time, if the solution achieves that lower bound. The lower bound is not provided in the MPS file, and its detail depends on the optimizer's way to implement quadratic problems (see example section for how to obtain it). Users have to work out how to provide that bound to their optimizer.

Note that it can take a long time to write the mps files, if the problem has many variables (the number of variables is \code{prod(nlevels)}).
}
\value{
Function \code{create_ILPlist} creates a list of one or more optimization problems,
which can be used in a call to function \code{write_MPSILPlist}.

Function \code{write_MPSILPlist} returns a table of contents matrix for the written files. If \code{toc=TRUE}, that value is also written to a separate file with suffix \code{toc}.

Function \code{write_MPSMIQP} returns the start vector (in counting vector representation), if \code{start} is not NULL; that vector is also written to a separate file with suffix \code{start}.
}
\references{
Groemping, U. and Fontana, R. (2019). An Algorithm for Generating Good Mixed Level Factorial Designs. \emph{Computational Statistics & Data Analysis} \bold{137}, 101-114. doi:10.1016/j.csda.2019.01.020.

Mosek ApS (2017a). MOSEK version w.x.y.z documentation. Accessible at: \url{https://www.mosek.com/documentation/}. This package has been developed using version 8.1.0.23 (accessed August 29 2017).

Mosek ApS (2017b). MOSEK Rmosek Package 8.1.y.z. \url{https://docs.mosek.com/8.1/rmosek/index.html}.
\emph{!!! In normal R speak, this is the documentation of the Rmosek package version 8.0.69 (or whatever comes next), when applied on top of the Mosek version 8.1.y.z (this package has been devoloped with Mosek version 8.1.0.23 and will likely not work for Mosek versions before 8.1). !!! (accessed August 29 2017)}

}

\note{
It can take an optimizer a long time to confirm optimality after finding an optimum. For the quadratic problem, it is therefore very beneficial to provide an optimal value to the algorithm, where possible (see example section).

The functions are not meant for situations, for which a full factorial design would be huge. Even though the functions do not solve anything, MPS files will be very large and writing them will be quite slow for such cases.
}

\author{
Ulrike Groemping
}

\seealso{
See also \code{\link{create_MIQP}} and \code{\link{write_MPSILP}} for internal functions that support these exported functions, and \code{\link{dToCount}} and \code{\link{countToDmixed}} for switching back and forth between an array and its counting vector representation.
}
\examples{
###################################################################
## an array and its counting vector

## arrays (starting the coding with 1)
## and their counting vectors can be used interchangeably
myarr <- cbind(c(1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2),
               c(1,1,1,1,2,2,2,2,1,1,1,1,2,2,2,2),
               c(1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4),
               c(1,5,3,7,2,6,4,8,8,4,6,2,7,3,5,1))

## we want to see it w.r.t. a 2,2,4,8 level full factorial
## determine the counting vector representation of this array
##     nlevels is needed,
##     because the third column of myarr has only 2 levels
(myarr_cv <- dToCount(myarr, nlevels=c(2,2,4,8), startfrom1=TRUE))

###################################################################
## demo: counting vector represents the array runs
###################################################################
## full factorial in lexicographic order
fullfac <- ff(2,2,4,8) +  1   ### ff levels start with 0
##
## pick the selected runs from fullfac
selfac <- fullfac[which(myarr_cv==1),]
##
## order both variants in the same way and compare them
## (in this case, they are equal without reordering)
ord1 <- DoE.base::ord(selfac)  ## order them
ord2 <- DoE.base::ord(myarr)   ## order them
selfac[ord1,] == myarr[ord2,]
#######################################################

#######################################################
## We go for an array in 16 runs with four factors in
## 2,2,4,8 levels.

## Is a strength 2 oa feasible?
##
oa_feasible(16, c(2,2,4,8), 2)  ## FALSE
##
## consequence: use resolution 2 (=strength 1),
##              minimize number of words of length 2

problemlist <- create_ILPlist(16, nlevels = c(2,2,4,8), resolution = 2)
length(problemlist) ## 12 distinct search orders
names(problemlist[[3]])
problemlist[[3]][-2]  ## ILP is too long for printing
problemlist1 <- create_ILPlist(16, nlevels = c(2,2,4,8),
                  resolution = 2, search.orders = FALSE)
                  ## only the pre-specified search order
problemlist2 <- create_ILPlist(16, nlevels = c(2,2,4,8),
                  resolution = 2, orders = list(c(2,2,4,8),
                                                c(8,2,4,2)))
                  ## the two specified search orders
\dontrun{
write_MPSILPlist(prefix="miniprob", problemlist)
## writes miniprob01.mps, ..., miniprob12.mps and miniprob_toc.txt
write_MPSILPlist(prefix="miniprob", problemlist1, toc=FALSE)
## writes miniprob1.mps
}

## The MPS files can be read by various optimizers.
## The ILP problems aim for a feasible solution.
## Start values are possible, but usually not useful.
## The best solution (lowest target value) can be imported into R.

## the solution  a counting vector
## its format depends on the optimizer
## import it into R and calculated array from it
importedsol <- myarr_cv # for demo only
solarray <- countToDmixed(myarr_cv, nlevels=c(2,2,4,8))
##
## it is crucial to use the order of the levels
## that corresponds to the problem that the solver solved

GWLP(solarray)

#######################################################
## providing a lower bound for the number of
## length 2 words in a strength 1 (resolution 2) array
#######################################################
##
lowerbound_AR(nruns = 16, nlevels = c(2,2,4,8), R = 2)  # 1
##
## In this example, we have immediately hit on a solution
## with optimum A2-value (see GWLP)

#######################################################
## using a quadratic problem for optimizing A2
##
\dontrun{
write_MPSMIQP("quadprob", 16, c(2,2,4,8), resolution=2)
## writes quadprob.mps
}

## Run time for solving the quadratic problem exported by write_MPSMIQP
## may substantially (!) benefit from providing the lower bound of the
## objective function, if that bound is attained.
##
## The lower bound for the minimum of the quadratic problem
##     created by write_MPSMIQP
##   is the lower bound for the word length, multiplied with n^2,
##     here 16 ^ 2 * 1 = 256,
##   or half that value,
##     depending on how the optimizer handles quadratic objectives.
#######################################################

## Depending on the optimizer, it is useful or even crucial to provide a
## starting value to write_MPSMIQP. This starting value can be obtained
## as the solution to a linear problem (that was exported using functions
## create_ILPlist and write_MPSILPlist).

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ array }% use one of  RShowDoc("KEYWORDS")
\keyword{ design }% __ONLY ONE__ keyword per line
