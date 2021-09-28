\name{functionsFromDoE.base}
\alias{oa_feasible}
\alias{lowerbound_AR}
\alias{length2}
\alias{length3}
\alias{length4}
\alias{length5}
\alias{contr.XuWu}
\alias{GWLP}
\alias{ICFTs}
\alias{SCFTs}

\title{
Functions from package DoE.base
}
\description{
These functions from DoE.base are exported from DoE.MIParray, 
because they are especially important for its use.
}
\usage{
oa_feasible(nruns, nlevels, strength = 2, verbose = TRUE, returnbound = FALSE)
lowerbound_AR(nruns, nlevels, R, crit = "total")
length2(design, with.blocks = FALSE, J = FALSE)
length3(design, with.blocks = FALSE, J = FALSE, rela = FALSE)
length4(design, with.blocks = FALSE, separate = FALSE, J = FALSE, rela = FALSE)
length5(design, with.blocks = FALSE, J = FALSE, rela = FALSE)
contr.XuWu(n, contrasts=TRUE)
GWLP(design, ...)
SCFTs(design, digits = 3, all = TRUE, resk.only = TRUE, kmin = NULL, kmax = ncol(design),
   regcheck = FALSE, arft = TRUE, cancors = FALSE, with.blocks = FALSE)
ICFTs(design, digits = 3, resk.only = TRUE, kmin = NULL, kmax = ncol(design), 
   detail = FALSE, with.blocks = FALSE, conc = TRUE)
}
\arguments{
  \item{nruns}{see \pkg{\link{DoE.base}}}
  \item{nlevels}{see \pkg{\link{DoE.base}}}
  \item{strength}{see \pkg{\link{DoE.base}}}
  \item{verbose}{see \pkg{\link{DoE.base}}}
  \item{returnbound}{see \pkg{\link{DoE.base}}}
  \item{R}{see \pkg{\link{DoE.base}}}
  \item{crit}{see \pkg{\link{DoE.base}}}
  \item{design}{see \pkg{\link{DoE.base}}}
  \item{with.blocks}{see \pkg{\link{DoE.base}}}
  \item{J}{see \pkg{\link{DoE.base}}}
  \item{rela}{see \pkg{\link{DoE.base}}}
  \item{n}{see \pkg{\link{DoE.base}}}
  \item{contrasts}{see \pkg{\link{DoE.base}}}
  \item{separate}{see \pkg{\link{DoE.base}}}
  \item{digits}{see \pkg{\link{DoE.base}}}
  \item{all}{see \pkg{\link{DoE.base}}}
  \item{resk.only}{see \pkg{\link{DoE.base}}}
  \item{kmin}{see \pkg{\link{DoE.base}}}
  \item{kmax}{see \pkg{\link{DoE.base}}}
  \item{regcheck}{see \pkg{\link{DoE.base}}}
  \item{arft}{see \pkg{\link{DoE.base}}}
  \item{cancors}{see \pkg{\link{DoE.base}}}
  \item{detail}{see \pkg{\link{DoE.base}}}
  \item{conc}{see \pkg{\link{DoE.base}}}
  \item{...}{see \pkg{\link{DoE.base}}}
}
\details{
for documentation of the functions, see the links under "See also"
}
\value{
for documentation of the functions, see the links under "See also"
}
\references{
for documentation of the functions, see the links under "See also"
}
\author{
Ulrike Groemping
}

\seealso{
See also \code{\link[DoE.base]{oa_feasible}}, \code{\link[DoE.base]{lowerbound_AR}}, 
\code{\link[DoE.base]{length2}}, \code{\link[DoE.base]{length3}}, \code{\link[DoE.base]{length4}}, 
\code{\link[DoE.base]{length5}}, \code{\link[DoE.base]{GWLP}}, \code{\link[DoE.base]{SCFTs}}, 
\code{\link[DoE.base]{ICFTs}}.
}
\examples{
oa_feasible(24, c(2,3,4,6),2)
lowerbound_AR(24, c(2,3,4,6),2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ design }% use one of  RShowDoc("KEYWORDS")
\keyword{ array }% __ONLY ONE__ keyword per line
