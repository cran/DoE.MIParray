create_ILPlist <- function(nruns, nlevels, resolution=3,
              distinct=TRUE, search.orders=TRUE, start=NULL, forced=NULL, orders=NULL){
  ## problem is formulated for minimization
  ## this is only relevant for strength=0 (resolution=1)
       ## think about separating that case
  ## nlevels must be in ascending order, if orders is used
  kmax <- max(resolution, 2)
  qco1 <- list(sense="min")
  nlev <- nlevels
  ## INTPNT_CO_TOL_PFEAS helps, if mosek declares previous start value as infeasible
     ## default 10^-8 is too strict, in the 48 2.4.3.4.2 example even 10^-6 is too strict

  ## the function creates the optimization problem and exports it in MPS format

  ## check inputs

  if (!is.numeric(nruns)) stop("nruns must be an integer number")
  if (!is.numeric(nlev)) stop("nlev must have integer elements")
  nfac <- length(nlev)
  if (nfac < 2) stop("nlev must have at least two elements")
  if (!is.numeric(resolution)) stop("resolution must be an integer numcber")
  if (!is.numeric(kmax)) stop("kmax must be an integer number")

  if (!length(nruns)==1) stop("nruns must be scalar")
  if (!length(resolution)==1) stop("resolution must be scalar")
  strength <- resolution - 1
  if (!is.logical(distinct)) stop("distinct must be logical")
  if (!is.logical(search.orders)) stop("search.orders must be logical")

  if (distinct && nruns>prod(nlev)) stop("too many runs for design with distinct runs")

  if (!is.null(start)){
    if (search.orders) stop("a start array cannot be provided if search.orders is TRUE")
    if (strength==0) stop("a start array requires resolution > 1")
    if (!is.numeric(start)) stop("start must be numeric")
    if (!all(start%%1==0)) stop("start must have integer elements")
    if (!all(start>=0)) stop("start must have non-negative elements")
    if (is.matrix(start)){
      if (!all(dim(start)==c(nruns, nfac))) stop("matrix start has wrong dimensions")
      if (!all(levels.no(start)==nlev)) stop("start and nlevels do not match")
      for (i in 1:length(nlev)) if (length(setdiff(start[,i],1:nlev[i]))>0)
        stop("invalid entries in column ", i, " of matrix start")
      start <- dToCount(start-1)
    }
    if (!length(start)==prod(nlev)) stop("vector start has wrong length")
    if (!sum(start)==nruns) stop("vector start is incompatible with nruns")
    if (!round(DoE.base::GWLP(countToDmixed(nlev, start), kmax=strength), 8)[strength+1]==0)
      stop("resolution of start array too low")
    if (distinct && !all(start %in% c(0,1)))
      stop("start array does not comply with option distinct")
    ## fixed text in May 2021, "not" was missing
  }

  if (!is.null(forced)){
    if (strength==0) stop("forcing runs requires resolution > 1")
    if (!is.numeric(forced)) stop("forced must be numeric")
    if (!all(forced%%1==0)) stop("forced must have integer elements")
    if (!all(forced>=0)) stop("forced must have non-negative elements")
    if (is.matrix(forced)){
      if (!ncol(forced)==nfac) stop("matrix forced has wrong number of columns")
      if (nrow(forced)>=nruns) stop("matrix forced has too many rows")
      if (!all(levels.no(forced)==nlev)) stop("forced and nlevels do not match")
      for (i in 1:length(nlev)) if (length(setdiff(forced[,i],1:nlev[i]))>0)
        stop("invalid entries in column ", i, " of matrix forced")
      forced <- dToCount(forced-1)
    }
    if (!length(forced)==prod(nlev)) stop("vector forced has wrong length")
    if (!sum(forced)<nruns) stop("vector forced fixes all runs", "use start for a start array")
    if (distinct && !all(forced %in% c(0,1)))
      stop("forced array does comply with option distinct")
    if (!is.null(start)){
      if (!all(forced <= start)) stop("start array does not contain all forced runs.")
    }
  }

if (search.orders){
    stopifnot(is.null(orders) || is.list(orders))
    hilf <- sort(nlev)
    if (!is.null(orders))
      stopifnot(all(sapply(lapply(orders, sort), function(obj) all(obj ==
                                                                   hilf))))
    if (is.null(orders))
      orders <- unique(combinat::permn(nlevels))
  } else
    orders <- list(nlevels)

  forced.orig <- forced  ## needed because of potential search orders

    ## preliminary exclusion of definitely infeasible cases
    N <- prod(nlev)
    if (strength > 0)
      if (!DoE.base::oa_feasible(nruns, nlev, strength)) stop("requested array does not exist")

  problist <- vector(mode="list")

  ## orders is now a list
  ## and is looped over, even if no search is conducted
  ## loop over orders
  for (l in 1:length(orders)){
    nlev <- orders[[l]]
    if (!is.null(forced)){
      forced <- rearrange_forced(forced.orig, nlevels, nlev)
    }

  ## preparation of matrices needed in formulating the optimization problem
  df1 <- nlev-1
  D <- ff(nlev)
  df <- 1
  for (j in 1:nfac) df <- c(df,sum(apply(matrix(df1[nchoosek(nfac,j)],nrow=j),2,prod)))
  dfweg <- cumsum(df)
  Dfac <- as.data.frame(D)
  for (j in 1:ncol(Dfac)){
    Dfac[[j]] <- factor(Dfac[[j]])
    contrasts(Dfac[[j]]) <- DoE.base::contr.XuWu(nlev[j])
  }

  ## model matrix in normalized orthogonal coding
  mm <- model.matrix(formula(paste("~.^",kmax,sep="")), Dfac)

  ## prepare H and U matrices (Hs[[i]]==t(Us[[i]])%*%Us[[i]])
  ## Us contains transposed model matrices for main effects, 2fis, 3fis etc.
  ## Hs contains their inner products
  Us <- vector(mode="list")
  Hs <- vector(mode="list")
  ## for lower triangular representation of Hs entries
  #i <- unlist(lapply(1:N, function(obj) obj:N))
  #j <- unlist(lapply(N:1, function(obj) rep(N+1-obj,obj)))
  for (kk in 1:kmax){
    Us[[kk]] <- t(mm[,(dfweg[kk]+1):dfweg[kk+1]]) ## crossprod is H
    #Hs[[kk]] <- round(tcrossprod(Us[[kk]]),8)
  }


  ## linear constraint: sum equal to nruns
  qco1$A <- Matrix::Matrix(rep(1,N),nrow=1)      ## constrain sum
  qco1$bc <- matrix(nruns, 2,1)                  ## equal to nruns

  ## implement forced, if applicable
  if (!is.null(forced)){
    ## additional rows of A
    hilf <- which(forced>0)
    mathilf <- matrix(0, length(hilf), ncol(qco1$A))
    mathilf[,hilf] <- diag(length(hilf))  ## identity
    if (distinct) qco1 <- mosek_modelAddLinear(qco1, mathilf, forced[hilf])
    else qco1 <- mosek_modelAddLinear(qco1, mathilf, forced[hilf], sense=">=")
  }

  ## variable related dimension
  if (distinct)
    qco1$bx <- rbind(rep(0,N), rep(1,N))
  else qco1$bx <- rbind(rep(0,N), rep(Inf,N))         ## all nonnegative
  qco1$c <- rep(0,N)                             ## objective is constant
  qco1$intsub <- 1:N  ## give indices            ## all integer

  ## info elements
  info <- list(nruns = nruns,
              nlev = nlev,
              reso = resolution)
  class(qco1) <- "qco"

  ## enforced resolution by linear equality constraints
  if (strength==0){
    ## not even balance enforced, cone variables needed for A1 minimization
    qco1 <- mosek_modelAddLinear(qco1, Us[[1]])
    qco1 <- mosek_modelAddConeQobj(qco1, df[2])
  }
  else
    for (ii in 1:strength)
      qco1 <- mosek_modelAddLinear(qco1, Us[[ii]])
  qco1$info <- NULL
  problist[[l]] <- list(nlevels=nlev, ILP=qco1, info=info)
  }
    return(problist)
}
