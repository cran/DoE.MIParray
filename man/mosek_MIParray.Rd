\name{mosek_MIParray}
\alias{mosek_MIParray}
\alias{gurobi_MIParray}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Functions to Create a MIP Based Array Using Gurobi or Mosek
}
\description{
The functions create an array with specified minimum resolution and optionally also optimized word length pattern based on mixed integer programming with the commercial software Gurobi (free academic license available) or Mosek (free academic license available).
Creation is done from scratch, or using a user-specified starting value, or extending an existing array. Important: Installation of Gurobi and/of Mosek as well as the corresponding R packages is necessary. The R package gurobi comes with the software, the current version of the R package Rmosek has to be obtained from vendor's website (CRAN version is outdated!).
}
\usage{
mosek_MIParray(nruns, nlevels, resolution = 3, kmax = max(resolution, 2),
   distinct = TRUE, detailed = 0, start=NULL, forced=NULL, find.only = FALSE,
   maxtime = Inf, nthread=2, mosek.opts = list(verbose = 10, soldetail = 1),
   mosek.params = list(dparam = list(LOWER_OBJ_CUT = 0.5, MIO_TOL_ABS_GAP = 0.2,
      INTPNT_CO_TOL_PFEAS = 1e-05, INTPNT_CO_TOL_INFEAS = 1e-07),
      iparam = list(PRESOLVE_LINDEP_USE="OFF", LOG_MIO_FREQ=100)))
gurobi_MIParray(nruns, nlevels, resolution = 3, kmax = max(resolution, 2),
  distinct = TRUE, detailed = 0, start=NULL, forced=NULL, find.only = FALSE,
  maxtime = 60, nthread = 2, heurist=0.5, MIQCPMethod=0, MIPFocus=1,
  gurobi.params = list(BestObjStop = 0.5, LogFile=""))
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
  \item{kmax}{
integer, \code{kmax >= resolution} and \code{kmax >= 2} are required; the largest number of words to be optimized (default: \code{kmax = resolution}). (If \code{find.only=TRUE}, optimization of numbers of words is suppressed.)
}
  \item{distinct}{
logical; if TRUE (default), restricts counting vector to 0/1 entries, which means that the resulting array is requested to have distinct rows; otherwise, duplicate rows are permitted, i.e. the counting vector can have arbitrary non-negative integers. Designs with distinct runs are usually better; in addition, binary variables are easier to handle by the optimization algorithm. Nevertheless, there are occasions where a better array is found faster with option \code{distinct=FALSE}, even if it has distinct rows.
}
  \item{detailed}{
integer (default 0); determines the output detail: positive values imply inclusion of a problem and solution history (attribute \code{history}), values of at least 3 add the lists of optimization matrices (Us and Hs, attribute \code{matrices}).
}
  \item{start}{
for \code{resolution > 1} only;\cr
a starting value for the algorithm: can be a array matrix with entries 1 to number of levels for each column, or a counting vector for the full factorial in lexicographic order; if specified, \code{start} must specify an array with the appropriate number of rows and columns, the requested resolution and, if \code{distinct = TRUE}, also contain distinct rows (matrix) or 0/1 elements only.
}
  \item{forced}{
for \code{resolution > 1} only;\cr
runs to force into the solution design; can be given as an array matrix with the appropriate number of columns and less than \code{nruns} rows or a counting vector for the full factorial in lexicographic order with sum smaller than \code{nruns}; if \code{distinct=TRUE}, \code{forced} must have distinct rows (matrix) or 0/1 elements only.
}
  \item{find.only}{
logical; if FALSE (default), a design of the requested resolution is found, which is subsequently improved in terms of its word lengths up to words of length \code{kmax}; otherwise, the function only attempts to find an array of the requested resolution, without optimizing word lenghts.
}
  \item{maxtime}{
the maximum run time in seconds per Gurobi or Mosek optimization request (the overall run time may become (much) larger); in case of conflict between \code{maxtime} and an explicit timing request in
\code{gurobi.params$TimeLimit} or \code{mosek.params$dparam$$MIO_MAX_TIME}, the stricter request prevails; the default values differ between Gurobi (60) and Mosek (Inf), because Mosek runs can be easily escaped, while Gurobi runs cannot.
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
different seeds modify the path taken through the search space. Varying the seed may be an alternative to searching over different level orderings with function \code{\link{mosek_MIPsearch}}.

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
The functions initially solve a linear optimization problem for obtaining a design with the requested resolution;
if a start value is provided, this step is skipped; if \code{find.only=TRUE}, the functions stop after this step.
If \code{find.only=FALSE}, the number of shortest words is optimized, followed by the numbers of words of length up to \code{kmax}.
The argument \code{forced} allows to specify an existing array that is to be extended (e.g. to double or triple size;
extension by a small number of runs will usually not be possible) in an optimized way.

For all but very small problems, it is likely advisable to choose \code{kmax} equal to the requested resolution (the default), and to proceed to longer words only after it has been made sure that the shortest word length has been optimized (as far as possible with reasonable effort). Further improvements of the same can be attempted by applying \code{gurobi_MIPcontinue} or \code{mosek_MIPcontinue} to the result object returned by the function. Note that it is possible to switch from using Mosek to using Gurobi or vice versa.
An example for an optimization sequence can be found in the package overview at \code{\link{DoE.MIParray}}.

In case of long run times, escaping from the gurobi run will most likely be unsuccessful and might even leave R in an unstable state; thus, one should think carefully about the affordable run times. On the contrary, it should usually be doable to escape a Mosek run; the remaining code of the function will still be executed and will return the final state.

Besides the run time, the number of threads is a very important resource control parameter. The default assumes that the user wants to use two of the computer's multiple cores. For using Gurobi's or Mosek's default (which is in most cases the use of all available cores), choose \code{nthread=0}.

The default Gurobi parameters have been chosen after systematic experimentation with a limited set of scenarii und Gurobi version 7.5.1.  The choice of \code{MIQCPMethod=0} was instrumental in many cases; changing it to -1 may occasionally be tried. The MIPFocus parameter also appears beneficial in many cases; changing it to 0 (leave choice to Gurobi) can be an option; it was slightly better than choice 1 for mixed level designs with relatively small run sizes, while choice \dQuote{1} was substantially better for the other cases. The heuristics proportion has been chosen as 0.5, because this choice seemed the best compromise for the situations considered. Note, however, that these parameters deteriorate performance for very simple cases, e.g. the test cases of Fontana (2017). For such cases, using \code{MIQCPMethod=-1}, \code{MIPFocus=0} and \code{heurist=0.05} will be preferrable; the defaults were chosen in this way, since doubling or even tripling very short run times was decided to be less detrimental than making more difficult problems completely intractable.

For Gurobi, several optimization parameters are switched off for the initial linear optimization step:
the parameters \code{Heuristics} and \code{MIPFocus} are reset to their defaults (0.05 and 0).

Gurobi always stores the file "gurobi.log" in the working directory; even if storage of the log is
suppressed with the default option \code{LogFile=""} or directed to another location by specifying a path, the default file "gurobi.log" is created and filled with a small amount of content. Thus, make sure to use a different file name when intentionally storing some log.

For Mosek, storing log output can be accomplished by directing the printed output to a suitable storing location. Note that the setting \code{iparam$LOG_MIO_FREQ = 100} reduces the frequency of printing a log line for branch-and-cut optimization by the factor 10 versus the default. Parallelization in Mosek is not well-protected against interference from screen activity (for example). Thus, one should switch off logging to screen or otherwise, when working with many (all available) threads in parallel (\code{LOG=0} instead of \code{LOG_MIO_FREQ = 100} in the list \code{iparam}).
}
\value{
an array of class \code{\link[DoE.base]{oa}}, possibly with the following attributes:
\code{MIPinfo}, which is either an object of class \code{qco} or a simple list with information (which would be the \code{info} element of the object of class \code{qco} in case the last optimization was not successful), \code{history} as a list of problem and solution lists, and \code{matrices} as a list of matrix lists. Presence or absence of \code{history} and \code{matrices} is controlled by option \code{detailed}, while \code{MIPInfo} is present if the optimization can be potentially improved by improving the last step (stop because of time limit and not because of optimal value) or by improving the number of longer words.
}
\section{Installation}{
Gurobi and Mosek need to be separately installed; please follow vendors' instructions; see also \code{\link{mosek_MIPsearch}} for more comments.
}
\references{
Fontana, R. (2017). Generalized Minimum Aberration mixed-level orthogonal arrays: a general approach based on sequential integer quadratically constrained quadratic programming. \emph{Communications in Statistics – Theory Methods} \bold{46}, 4275-4284.

Groemping, U. and Fontana R. (2019). An Algorithm for Generating Good Mixed Level Factorial Designs. \emph{Computational Statistics & Data Analysis} \bold{137}, 101-114.

Groemping, U. (2020). DoE.MIParray: an R package for algorithmic creation of orthogonal arrays. \emph{Journal of Open Research Software}, \bold{8: 24}. DOI: https://doi.org/10.5334/jors.286

Gurobi Optimization Inc. (2017). Gurobi Optimizer Reference Manual. \url{https://www.gurobi.com:443/documentation/}.

Mosek ApS (2017a). MOSEK version w.x.y.z documentation. Accessible at: \url{https://www.mosek.com/documentation/}. This package has been developed using version 8.1.0.23 (accessed August 29 2017).

Mosek ApS (2017b). MOSEK Rmosek Package 8.1.y.z. \url{https://docs.mosek.com/8.1/rmosek/index.html}.
\emph{!!! In normal R speak, this is the documentation of the Rmosek package version 8.0.69 (or whatever comes next), when applied on top of the Mosek version 8.1.y.z (this package has been devoloped with Mosek version 8.1.0.23 and will likely not work for Mosek versions before 8.1). !!! (accessed August 29 2017)}

}

\note{
The functions are not meant for situations, for which a full factorial design would be huge; the mixed integer problem to be solved has at least \code{prod(nlevels)} binary or general integer variables and will likely be untractable, if this number is too large. (For extending an existing designs, since some variables are fixed, the limit moves out a bit.)

Please be aware that escaping a Gurobi run will be not unlikely to leave the computer in an unstable situation.  If function \code{gurobi_MIParray} is successfully interrupted by the <ESC> key or <Ctrl>-<C>, it will usually be necessary to restart R in order to free all CPU usage.

If a Mosek run is interrupted by the <ESC> key, it is advised to execute the command \code{Rmosek::mosek_clean()} afterwards; this may help prevent problems from unclean closes of mosek runs.

}

\author{
Ulrike Groemping
}

\seealso{
See also \code{\link{mosek_MIPsearch}} and \code{\link{gurobi_MIPsearch}} for searching over \code{nlevels} orderings,
\code{\link{mosek_MIPcontinue}} and \code{\link{gurobi_MIPcontinue}} for continuing an uncompleted optimization,
\code{\link[DoE.base]{show.oas}} from package \pkg{DoE.base} for catalogued orthogonal arrays,
and \code{\link[DoE.base]{oa_feasible}} from package \pkg{DoE.base} for checking feasibility of requested array strength (resolution - 1) for combinations of \code{nruns} and \code{nlevels}.
}
\examples{
\dontrun{
## can also be run with gurobi_MIParray instead of mosek_MIParray
## there are of course better ways to obtain good arrays for these parameters
## (e.g. function FrF2 from package FrF2)
feld <- mosek_MIParray(16, rep(2,7), resolution=3, kmax=4)
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
