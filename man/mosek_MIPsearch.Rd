\name{mosek_MIPsearch}
\alias{mosek_MIPsearch}
\alias{gurobi_MIPsearch}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Functions to Search for optimum MIP Based Array Using Gurobi or Mosek
}
\description{
The functions search through different orderings of the nlevels vector 
with the goal to create an array with minimum resolution and optimized shortest 
word length. They create the orders and call gurobi_MIParray or mosek_MIParray 
for each order.
}
\usage{
mosek_MIPsearch(nruns, nlevels, resolution = 3, maxtime = 60, 
   stopearly=TRUE, listout=FALSE, orders=NULL, 
   distinct = TRUE, detailed = 0, forced=NULL, find.only=TRUE,
   nthread=2, mosek.opts = list(verbose = 1, soldetail = 1),
   mosek.params = list(dparam = list(LOWER_OBJ_CUT = 0.5, MIO_TOL_ABS_GAP = 0.2,
      INTPNT_CO_TOL_PFEAS = 1e-05, INTPNT_CO_TOL_INFEAS = 1e-07),
      iparam = list(PRESOLVE_LINDEP_USE="OFF", LOG_MIO_FREQ=100)))
gurobi_MIPsearch(nruns, nlevels, resolution = 3, maxtime = 60, 
   stopearly=TRUE, listout=FALSE, orders=NULL, 
   distinct = TRUE, detailed = 0, forced=NULL, find.only=TRUE, 
   nthread = 2, heurist=0.5, MIQCPMethod=0, MIPFocus=1,
   gurobi.params = list(BestObjStop = 0.5, OutputFlag=0))
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
  \item{maxtime}{
the maximum run time in seconds per Gurobi or Mosek optimization request (the overall run time may become (much) larger); in case of conflict between \code{maxtime} and an explicit timing request in
\code{gurobi.params$TimeLimit} or \code{mosek.params$dparam$$MIO_MAX_TIME}, the stricter request prevails; the default values differ between Gurobi (60) and Mosek (Inf), because Mosek runs can be easily escaped, while Gurobi runs cannot.
}
\item{stopearly}{
  logical; if TRUE, the search stops if the shortest word length hits the lower bound; set to FALSE if you want longer word lengths to be optimized among several choices with the same shortest word length
}
\item{listout}{
  logical; if TRUE, all experimental plans are stored, instead of only the best one; if \code{stopearly=TRUE}, \code{listout=TRUE} does not make sense
}
\item{orders}{
  NULL (in which case distinct level orders are automatically determined) or a list of level orders to be searched
}
  \item{distinct}{
logical; if TRUE (default), restricts counting vector to 0/1 entries, which means that the resulting array is requested to have distinct rows; otherwise, duplicate rows are permitted, i.e. the counting vector can have arbitrary non-negative integers. Designs with distinct runs are usually better; in addition, binary variables are easier to handle by the optimization algorithm. Nevertheless, there are occasions where a better array is found faster with option \code{distinct=FALSE}, even if it has distinct rows.
}
  \item{detailed}{
integer (default 0); determines the output detail: positive values imply inclusion of a problem and solution history (attribute \code{history}), values of at least 3 add the lists of optimization matrices (Us and Hs, attribute \code{matrices}).
}
  \item{forced}{
for \code{resolution > 1} only;\cr
runs to force into the solution design; can be given as an array matrix with the appropriate number of columns and less than \code{nruns} rows or a counting vector for the full factorial in lexicographic order with sum smaller than \code{nruns}; if \code{distinct=TRUE}, \code{forced} must have distinct rows (matrix) or 0/1 elements only.
}
  \item{find.only}{
logical; if TRUE (default), the function only attempts to find an array of the requested resolution, without optimizing word lenghts; otherwise, a design of the requested resolution is found, which is subsequently improved in terms of its word lengths up to words of length \code{kmax}.
}
  \item{nthread}{
the number of threads (=cores) to use;
there are also the Mosek parameter \code{NUM_THREADS} and the Gurobi parameter \code{Threads}; in case of conflict, the smaller request prevails.
For using Gurobi's or Mosek's default (which is in most cases the use of all available cores), choose \code{nthread=0}.\cr
CAUTION: \code{nthread} should not be chosen larger than the available number of cores. Gurobi warns that performance will deteriorate, but was observed to perform OK. For Mosek, performance will strongly deteriorate, and for extreme choices the R session might even crash (even for small problems)!
}
\item{mosek.opts}{
list of Mosek options; these have to be looked up in Mosek documentation
}
\item{mosek.params}{
list of mosek parameters, which can have the list-valued elements \code{dparam}, \code{iparam} and/or \code{sparam}; their use has to be looked up in the RMosek documentation.
The arguments \code{maxtime} and \code{nthread} correspond to the \code{dparam$MIO_MAX_TIME} and
\code{iparam$NUM_THREADS} specifications. Conflicts are resolved as stated in their documentation.

The element \code{dparam$LOWER_OBJ_CUT} can be used to incorporate a best
bound found in an earlier successless optimization attempt; per default, it is set to 0.5, since
the target function can take on integer values only and cannot be negative.

If a valid starting value is not accepted by Mosek, it may be worthwhile to increase \code{dparam$INTPNT_CO_TOL_PFEAS}.

Users of Mosek versions 9 and higher may want to play with \code{iparam$MIO_SEED}, which was introduced as a new parameter with Mosek version 9 (default: 42); 
different seeds modify the path taken through the search space for a given level ordering; thus, varying seeds can also be the route to choose where searching over level 
orderings is not feasible. 

Note that a user specified \code{mosek.params} should always contain the specifications shown under Usage. Exceptions: \code{LOWER_OBJ_CUT} is always specified to be at least 0.5, i.e. this option can be safely omitted without loosing anything, and intentional changes can of course be made.
}
  \item{heurist}{
the proportion heuristics time used by Gurobi in quadratic objective optimization (default 0.5; Gurobi default is 0.05);
there is also the Gurobi parameter Heuristics; in case of conflict, the larger request prevails;
the setting for heurist is deactivated for the initial linear problem which is always run with the Gurobi default. It can be worthwhile playing with this option for improving the run time for certain settings; for example, with \code{nruns=48} and \code{nlevels=c(2,2,3,4,4)}, \code{heurist=0.05} performs better than the default 0.5.
}
  \item{MIQCPMethod}{
the method used by Gurobi for quadratically constrained optimization (default 0; other possibilities -1 (Gurobi decides) or 1);
there is also the Gurobi parameter MIQCPMethod; in case of conflict, the method is set to "0"; this choice is made because it proved beneficial in many cases explored (although there also were a few cases which fared better with Gurobi's default).
}
  \item{MIPFocus}{
the strategy used by Gurobi for quadratically constrained optimization (default 1: focus on finding good feasible solutions fast; other possibilities: 0 (Gurobi decides/compromise), 2 or 3 (focus on increasing the lower bound fast));
there is also the Gurobi parameter MIPFocus; in case of conflict, MIPFocus is set to "0";
the setting for MIPFocus is deactivated for the initial linear problem which is always run with the Gurobi default.
}
  \item{gurobi.params}{
list of gurobi parameters; these have to be looked up in Gurobi documentation;
the arguments \code{maxtime}, \code{heurist}, \code{MIQCPMethod} and \code{MIPFocus} refer to the Gurobi parameters "TimeLimit", "Heuristics", "MIQCPMethod" and "MIPFocus", respectively. See their documentation for what happens in case of conflict.

The Gurobi parameter \code{BestObjStop} can be used to incorporate a best
bound found in an earlier successless optimization attempt; per default, it is set to 0.5, since
the objective function can take on integer values only and cannot be negative.
}
}
\details{
The search functions have been implemented, because the algorithm's behavior may 
strongly depend on the order of factors in case of mixed level arrays. 
In many examples, Mosek quickly improved the objective function which 
then stayed constant for a long time; thus, it may be promising to run 
\code{mosek_MIPsearch} with \code{maxtime=60} (or even less). See also 
Groemping and Fontana (2019) for examples of successful applications of the search 
functionality.

Even though Gurobi was less successful as a search tool in the examples that 
were examined so far, it may be helpful for other examples.

The options suppress printed output from the optimizers themselves.

Mosek Version 9 has gained a seed argument (\code{iparam$MIO_SEED}, 
which implements the Mosek parameter \code{MSK_IPAR_MIO_SEED}). 
Playing with seeds in \code{mosek_MIParray} may be an alternative to using the search approach, 
because it may lead to different paths through the search space for a fixed ordering of the \code{nlevels} vector. 
So far, I have only very little experience with using seeds; user reports are very welcome.
}
\value{
an array of class \code{\link[DoE.base]{oa}} with the attributes added 
by \code{mosek_MIParray} or \code{gurobi_MIParray}, resp. \cr
In addition, the attribute \code{optorder} contains the vector of level orders that yielded the best design; 
if \code{listout=TRUE}, also the attributes \code{orders} and \code{allplans}.

Objects with the attribute \code{allplans} are quite large. If the attribute is 
no longer needed, it can be removed from an object named \code{obj} 
(replace with the name of your object) by the command \cr
\code{attr(obj, "allplans") <- NULL}
}

\author{
Ulrike Groemping
}

\references{
Groemping, U. and Fontana R. (2019). An Algorithm for Generating Good Mixed Level Factorial Designs. \emph{Computational Statistics & Data Analysis} \bold{137}, 101-114.

Groemping, U. (2020). DoE.MIParray: an R package for algorithmic creation of orthogonal arrays. \emph{Journal of Open Research Software}, \bold{8: 24}. DOI: https://doi.org/10.5334/jors.286
}
\seealso{
See also \code{\link{mosek_MIParray}} and \code{\link{gurobi_MIParray}},
\code{\link[DoE.base]{oa_feasible}} from package \pkg{DoE.base} for checking 
feasibility of requested array strength (resolution - 1) 
for combinations of \code{nruns} and \code{nlevels}, and 
\code{\link[DoE.base]{lowerbound_AR}} from package \pkg{DoE.base} 
for a lower bound for the length R words in a resolution R array.
}
\examples{
\dontrun{
## can also be run with gurobi_MIParray instead of mosek_MIParray
## there are of course better ways to obtain good arrays for these parameters
## (e.g. function FrF2 from package FrF2)
oa_feasible(18, c(2,3,3,3,3), 2)  ## strength 2 array feasible
lowerbound_AR(18, c(2,3,3,3,3), 3)  ## lower bound for A3
## of course not necessary here, the design is found fast
feld <- mosek_MIPsearch(18, c(2,3,3,3,3), stopearly=FALSE, listout=TRUE, maxtime=30)
## even stopearly=TRUE would not stop, because the lower bound 2 is not achievable
feld
names(attributes(feld))
attr(feld, "optorder")
## even for this simple case, running optimization until confirmed optimality 
## would be very slow
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ array }% use one of  RShowDoc("KEYWORDS")
\keyword{ design }% __ONLY ONE__ keyword per line
