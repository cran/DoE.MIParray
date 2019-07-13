\name{DoE.MIParray-package}
\alias{DoE.MIParray}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Package to Create a MIP Based Array
}
\description{
The package creates an array with specified minimum resolution and optimized word length pattern, using mixed integer programming (MIP). It requires the presence of at least one of the commercial softwares Gurobi (free academic license available, R package gurobi included with the software) or Mosek (free academic license available, outdated CRAN package Rmosek must be updated with current version provided by vendor). Arrays can be created from scratch, or using a user-specified starting array, or extending an existing array.
}
\details{
The package implements generalized minimum aberration (GMA) according to Xu and Wu (2001), which is a way to minimize the confounding potential in a factorial design, as measured by the generalized word length pattern (GWLP). Reduction of short words has priority over reduction of long words, because it pertains to confounding of lower order interactions, which is assumed to be more severe than confounding of higher order interactions.

There is always one word of length zero. An array is said to have strength t (resolution R=t+1), if it has no words of lengths 1,...,t, but has words of length R=t+1; in that case, t-factor interactions are not confounded with the overall mean, (t-1)-factor interactions are not confounded with main effects, and so forth. Groemping and Xu (2014) provided an interpretation of the number of shortest words (i.e. words of length R) in terms of coefficients of determination of linear models with main effects model matrix columns on the LHS and full factorial t-factor models on the right-hand side; for example, in an array of strength 2 with three 3-level factors, the number of words of length 3 is the sum of the two R^2 values obtained from regressing the two main effects model matrix columns of one of the factors on a full model in the other two factors, provided main effect model matrix columns are coded orthogonally (to each other and the overall mean). GMA considers a design better than another one, if it has larger strength or resolution. In case of ties, a design is better if it has fewer shortest words; further ties are resolved by comparing words of lengths t+2, t+3, ...

Any array found by a function of this package will have the requested resolution (if not possible, an error will be thrown), and optimization of the number of shortest words will be attempted. If \code{kmax} is chosen larger than the resolution, optimization of longer words will also be attempted. Mixed integer optimization is very resource-intensive and often fails to provide a confirmed optimum (see also below). Choosing \code{kmax} larger than the resolution is therefore advisable for very small problems only; in most cases, one should attempt an optimization of the number of shortest words only. Only after this has been achieved, possibly with several sequential attempts, a subsequent attempt to improve the number of longer words should be undertaken.

Functions \code{\link{gurobi_MIParray}} and \code{\link{mosek_MIParray}} create an array from scratch (\code{start=NULL} and \code{forced=NULL}), by trying to improve a starting array specified with the argument \code{start}, or by restricting a portion of the array to correspond to a pre-existing array with the argument \code{forced}. If no starting array is given, the functions initially obtain one by solving a linear optimization problem for obtaining a design with the requested resolution. Subsequently, the number of shortest words of the starting array is optimized, followed by the numbers of words of length up to \code{kmax}.

Where functions \code{\link{gurobi_MIParray}} and \code{\link{mosek_MIParray}} do not easily find an optimal array, it may be worthwhile to consider using functions \code{\link{gurobi_MIPsearch}} and \code{\link{mosek_MIPsearch}}, which can search over different orderings of the factor levels; these have been provided because a brief search for a fortunate level orderings may be successful where a very long search for an unfortunate level ordering fails.

The strategy used in the package is a modification of Fontana (2017). Modifications include enforcing the requested resolution via a linear optimization step (much faster than sequentially optimizing all those word lengths until they are zeroes), and using a result on coding invariance from Groemping (2018) for a more parsimonious formulation of constraints. Mixed integer programming can use a lot of time and resources; in particular, the confirmation of the optimality of a solution that has been found can take prohibitively long, even if the optimal solution itself has been found fast (which is also not necessarily so). Functions \code{\link{gurobi_MIPcontinue}} and \code{\link{mosek_MIPcontinue}} can be used to continue optimization for larger problems,
where optimization was previously aborted, e.g. due to a time limit (however, a lot of effort is lost and has to be repeated when continuing a previous attempt). Note that it is possible to continue an optimization effort started with Gurobi using Mosek and vice versa, because the functions \code{gurobi_MIPcontinue} and \code{mosek_MIPcontinue} can convert problems using the internal functions \code{\link{mosek2gurobi}} and \code{\link{gurobi2mosek}}. 


}

\references{
Fontana, R. (2017). Generalized Minimum Aberration mixed-level orthogonal arrays: a general approach based on sequential integer quadratically constrained quadratic programming. \emph{Communications in Statistics – Theory Methods} \bold{46}, 4275-4284.

Groemping, U. and Xu, H. (2014). Generalized resolution for orthogonal arrays. \emph{The Annals of Statistics} \bold{42}, 918-939.

Groemping, U. (2018). Coding Invariance in Factorial Linear Models and a New Tool for Assessing Combinatorial Equivalence of Factorial Designs. \emph{Journal of Statistical Planning and Inference} \bold{193}, 1-14.

Groemping, U. and Fontana R. (2019). An Algorithm for Generating Good Mixed Level Factorial Designs. \emph{Computational Statistics & Data Analysis} \bold{137}, 101-114.

Gurobi Optimization Inc. (2018). Gurobi Optimizer Reference Manual. \url{http://www.gurobi.com/documentation/}.

Mosek ApS (2017a). MOSEK version w.x.y.z documentation. Accessible at: \url{https://www.mosek.com/documentation/}. This package has been developed using version 8.1.0.23 (accessed August 29 2017).

Mosek ApS (2017b). MOSEK Rmosek Package 8.1.y.z. \url{http://docs.mosek.com/8.1/rmosek/index.html}.
\emph{!!! In normal R speak, this is the documentation of the Rmosek package version 8.0.69, when applied on top of the Mosek version 8.1.y.z (this package has been devoloped with Mosek version 8.1.0.23 and will likely not work for Mosek versions before 8.1). !!! (accessed August 29 2017)}

Xu, H. and Wu, C.F.J. (2001). Generalized minimum aberration for asymmetrical fractional factorial designs. \emph{Annals of Statistics} \bold{29}, 549-560.
}

\note{
The package is not meant for situations, for which a full factorial design would be huge; the mixed integer problem to be solved has at least \code{prod(nlevels)} binary or general integer variables and will likely be untractable, if this number is too large. (For extending an existing designs, since some variables are fixed, the limit moves out a bit.)
}
\section{Warning}{
Escaping from a Gurobi run will most likely be unsuccessful, might leave R in an unstable state, and usually fails to release the entire CPU usage; thus, one should think carefully about the affordable run times. \cr
Mosek can usually be escaped using the <ESC> key, and one can even hope to get a valid output. 
However, after such escapes, it is also advisable to use \pkg{RMosek}'s internal clean function
(\code{Rmosek:::mosek_clean()}), since computer instabilities after repeated escapes from Mosek have been observed
(the package functions execute the \code{mosek_clean} command after conducting Mosek runs).
}
\section{Installation}{
Gurobi and Mosek need to be separately installed; please follow vendors' instructions; it is necessary to obtain a license; for academic use, free academic licenses are available in both cases.

Both Gurobi and Mosek provide R packages (\pkg{gurobi} and \pkg{Rmosek}) for accessing the software, and they also provide instructions on their installation.

For Gurobi, the R package \pkg{gurobi} is provided with the software installation files and can simply be installed like usual for a non-Web package.

For Mosek, there is a CRAN package \pkg{Rmosek}, which is a stump only and can be used for installing the suitable package \pkg{Rmosek} from the Mosek website. 
Its installation requires compilation from source.  Several prerequisites are needed, basically the same ones needed for compiling packages. 
Make sure to set up the R environment for this purpose (for Windows, see \url{https://cran.r-project.org/bin/windows/Rtools/}; you may have to set some paths yourself; 
I am not familiar with the proper process for other platforms).  If this is accomplished and package \pkg{Matrix} is at the latest level (as recommended in the \pkg{Rmosek} 
online documentation), install the CRAN package \pkg{Rmosek}, whose purpose it is to support installation of the appropriate \pkg{Rmosek} package for your platform and Mosek version. 
Once this CRAN package is available, run function \code{mosek_attachbuilder}, specifying the appropriate path as 
pointed out in the function's documentation, and subsequently run the custom-made function \code{install.rmosek}.  After this activity (if all went well), the CRAN package \pkg{Rmosek} 
will have been replaced by a working version of package \pkg{Rmosek} whose version number depends on your version of Mosek.  In case of problems, 
running the \code{install.packages} command provided in the \pkg{Rmosek} online documentation for your version of Mosek (Mosek ApS 2017b) may also be worth a try; 
if all else fails, you may have to contact \pkg{Rmosek} support.

Don't be confused by Mosek ApS's somewhat strange communication regarding the package \pkg{Rmosek}: the CRAN version (stump for the purpose of supporting installation of the working version) 
has its own separate numbering that currently starts with \dQuote{1}. The working \pkg{Rmosek} packages have version numbers whose first digit is kept in sync with the Mosek major version number; 
apart from that, version numbers of \pkg{Rmosek} versions and Mosek versions are not kept in sync; for example, the current version (July 13 2019) of package Rmosek for Mosek 
version 8 is 8.0.69, while the Mosek version is 8.1.0.81 (the revision versions of Mosek change quite frequently).  Whenever Mosek itself is re-installed (at least with a change of minor version 
like 8.0 to 8.1), the \pkg{Rmosek} sources must be re-compiled, even if their version has not changed; 
this is presumably why Mosek ApS speak of a version number for the binary that corresponds to the Mosek version number, while the sources keep their version number 
(somewhat confusing for the R community).  From within R, version numbers can be queried by \code{packageVersion("Rmosek")} and \code{Rmosek::mosek_version()}, respectively.}
\author{
Ulrike Groemping
}

\seealso{
See also as \code{\link{gurobi_MIParray}}, \code{\link{mosek_MIParray}}, \code{\link{gurobi_MIPsearch}}, \code{\link{mosek_MIPsearch}}, \code{\link{gurobi_MIPcontinue}} and \code{\link{mosek_MIPcontinue}}. Furthermore, the following links from package \pkg{DoE.base} are of interest:
\code{\link[DoE.base]{show.oas}} queries catalogued orthogonal arrays and \code{\link[DoE.base]{oa_feasible}} carries out feasibility checks for runs, levels and strength requirements.
}
\examples{
\dontrun{
## ideal sequence of optimization problems
## shown here for Mosek,
## for Gurobi analogous, if necessary increasing maxtime to e.g. 600 or 3600 or ...

## very small problem
plan <- mosek_MIParray(16, rep(2,6), resolution=4, kmax=6)

## an example approach for a larger problem
## optimize shortest word length
plan3 <- mosek_MIParray(24, c(2,4,3,2,2,2,2), resolution=3, maxtime=20)
## feasible solution was found, no confirmed optimum, 7/3 words of length 3
## try to optimize further or confirm optimality (improve=TRUE does this),
##                      give it 10 minutes
plan3b <- mosek_MIPcontinue(plan3, improve=TRUE, maxtime=600)
##     no improvement has been found, and the gap is still very large
##     (the time limit makes the result non-deterministic, of course,
##      because it depends on the computer's power and availability of its resources)

## For large problems, it cannot be expected that a *confirmed* optimum is found.
## Of course, one can put more effort into the optimization, e.g. by running overnight.
## It is also advisable to compare the outcome to other ways for obtaining a good array,
##    e.g. function oa.design from package DoE.base with optimized column allocation.
require(DoE.base)
show.oas(nruns=24, nlevels=c(2,4,3,2,2,2,2), show=Inf)
GWLP(plan_oad <- oa.design(nruns=24, nlevels=c(2,4,3,2,2,2,2), col="min34"))
## here, plan3b has a better A3 than plan_oad

## one might also try to confirm optimality by switching to the other optimizer
plan3c <- gurobi_MIPcontinue(plan3b, improve=TRUE, maxtime=600, MIPFocus=3)
   ## focus on improved bound with option MIPFocus
   ## still same value with very large gap after running this
   ## thus, now assume this as best practically feasible value

## one might now try to improve words of length 4 (improve=FALSE turns to the next word length)
plan4 <- mosek_MIPcontinue(plan3b, improve=FALSE, maxtime=600)
   ## this does not yield any improvement
   ## working on longer words is not considered worthwhile
   ## thus, plan3 or plan3b are used for pragmatic reasons,
   ## without confirmed optimality
}
}
\keyword{ array }% use one of  RShowDoc("KEYWORDS")
\keyword{ design }% use one of  RShowDoc("KEYWORDS")
