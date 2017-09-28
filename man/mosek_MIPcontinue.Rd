\name{mosek_MIPcontinue}
\alias{mosek_MIPcontinue}
\alias{gurobi_MIPcontinue}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Functions to Continue Optimization from Stored State
}
\description{
These functions continue optimization for a MIP-based array from a stored state.
}
\usage{
mosek_MIPcontinue(qco, improve = TRUE, maxtime = Inf, nthread = 2,
      mosek.opts = list(verbose = 10, soldetail = 1),
      mosek.params = list(dparam = list(LOWER_OBJ_CUT = 0.5,
      MIO_TOL_ABS_GAP = 0.2, INTPNT_CO_TOL_PFEAS = 1e-05,
      INTPNT_CO_TOL_INFEAS = 1e-07),
      iparam = list(PRESOLVE_LINDEP_USE="OFF", LOG_MIO_FREQ=100)))
gurobi_MIPcontinue(qco, improve = TRUE, maxtime = 60, nthread = 2,
      heurist = 0.05, MIQCPMethod = 0, MIPFocus = 0,
      gurobi.params =list(BestObjStop = 0.5, LogFile=""))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{qco}{
   object of class \code{qco}, created by a function in the package; or object of class \code{oa} that has a \code{qco} object as its \code{MIPinfo} attribute
}
  \item{improve}{
   logical; if TRUE (default), try to improve the already obtained solution for word length \code{qco$info$last.k}, otherwise optimize the next word length
   }
  \item{maxtime}{
  time in seconds for the optimization call; defaults differ for Mosek (Inf) and Gurobi (60), because a Mosek run can be easily escaped (<ESC>-key), contrary to a Gurobi run
}
  \item{nthread}{
  number of cores to use (0=all cores)\cr
CAUTION: \code{nthread} should not exceed the available number of cores. Gurobi warns that performance might deteriorate. For Mosek, performance WILL strongly deteriorate, and for extreme choices the R session might even crash (even for small problems)!
}
  \item{heurist}{
  for \code{gurobi_MIPcontinue} only: the percentage of time (number between 0 and 1) spent on heuristics
}
  \item{MIQCPMethod}{
  for \code{gurobi_MIPcontinue} only: the choice of optimization method; the default "0" has been observed to be better than Gurobi's default for most cases (Gurobi version 7.5.1); "-1" leaves the choice to Gurobi, "1" chooses the other method
}
  \item{MIPFocus}{
  for \code{gurobi_MIPcontinue} only: the choice of strategy; the default "0" leaves this choice to Gurobi; for finding better feasible solutions, "1" is recommended; for improving the speed of increasing the lower bound for eventually proving optimality, "3" can be tried
}
  \item{mosek.opts}{
  mosek options
}
  \item{mosek.params}{
  Mosek parameters
}
  \item{gurobi.params}{
  Gurobi parameters
}
}
\details{
Note that it is possible to continue optimization with Gurobi, if it was started with Mosek, and vice versa. The tool will transform the problem into the respective other format.

Usage of options is analogous to functions \code{\link{mosek_MIParray}} and \code{\link{gurobi_MIParray}}, respectively, where these are described in more detail.

For some applications, usability of \code{mosek_MIPcontinue} is currently hampered by the fact that Mosek's presolve routines identify additional integer variables and fail to recognise user-specified starting values for these that are not exactly integer-valued. According to Mosek ApS, this is scheduled to be remedied with Mosek Version 9.
}
\value{
an array of class "oa", if not optimized to GMA with info for further continuation (see documentation of \code{\link{mosek_MIParray}} or \code{\link{gurobi_MIParray}})
}

\author{
Ulrike Groemping
}

\seealso{
See also \code{\link{DoE.MIParray}} for examples of the role of the \code{MIPcontinue} functions, \code{\link{mosek_MIParray}} and \code{\link{gurobi_MIParray}} for more detail on the optimization arguments.
}
