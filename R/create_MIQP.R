create_MIQP <- function(nruns, nlevels, resolution=3, distinct=TRUE, start=NULL, forced=NULL, forMosek=FALSE){
    ## start is expected to have N elements, i.e. to be a count vector for an initial design
    ## this should already have the resolution enforced
    ## start should be provided !!!
    ## problem is formulated for minimization
    ## this internal function creates an MPS problem with quadratic objective,
    ## for general MPS usage
    ## (Gurobi and Mosek would work on a linearized version with conic quadratic constraints,
    ## which is not part of the general MPS syntax)
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
    if (!is.numeric(resolution)) stop("resolution must be an integer number")
    if (!length(nruns)==1) stop("nruns must be scalar")
    if (!length(resolution)==1) stop("resolution must be scalar")
    strength <- resolution - 1
    kmax <- resolution
    if (kmax < 2) stop("This function only works for resolution at least 2.")
    if (!is.logical(distinct)) stop("distinct must be logical")
    if (distinct && nruns>prod(nlev)) stop("too many runs for design with distinct runs")

    if (is.null(start)) message("A start value for an array with suitable resolution might be helpful.")
    else{
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
        stop("start array does comply with option distinct")
         }
    if (!is.null(forced)){
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
      if (!sum(forced) < nruns) stop("vector forced fixes all runs", "use start for a start array")
      if (distinct && !all(forced %in% c(0,1)))
        stop("forced array does comply with option distinct")
      if (!is.null(start))
      if (!all(forced <= start)) stop("start array does not contain all forced runs.")
    }

    ## the stricter time request wins, LOWER_OBJ_CUT is always set
    ##      for A_resolution, the bound can and will be made
    ##      using hidden function lowerbounds from DoE.base
    ##      this happens in the loop

    ## preliminary exclusion of definitely infeasible cases
    N <- prod(nlev)
    if (strength > 0)
      if (!DoE.base::oa_feasible(nruns, nlev, strength)) stop("requested array does not exist")

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
    Hs[[kk]] <- round(tcrossprod(mm[,(dfweg[kk]+1):dfweg[kk+1]]),8)
    Us[[kk]] <- t(mm[,(dfweg[kk]+1):dfweg[kk+1]]) ## crossprod is H
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
  qco1$info$nruns <- nruns
  qco1$info$nlev <- nlev
  qco1$info$reso <- resolution
  qco1$info$last.k <- max(1,strength)
  qco1$info$optimizer <- "mosek"
  qco1$info$conecur <- FALSE
  qco1$info$start <- start  ## for the use in write_MPSMIQP
  class(qco1) <- "qco"

  for (ii in 1:strength)
       qco1 <- mosek_modelAddLinear(qco1, Us[[ii]])

    ## start array specified
  if (!is.null(start)){
    opt <- start
    sl <- vector(mode="list")
    sl$sol$int$xx <- opt
    ## should be 0
    ## improvement of word counts for lengths > R not implemented here
    vcur <- round(opt%*%Hs[[strength]]%*%opt)
    qco1$info$timelinear <- 0  ## no time in linear optimization
  }
  else {
    vcur <- Inf
    opt <- rep(NA, N)
  }

  vhistory <- list(objective=vcur, counts=cbind(round(opt[1:N])))
  nvar <- ncol(qco1$A)

  if (!is.null(start)){
    ## admissible?
    if (!all(round(qco1$A[-1,] %*% start, 8)==0))
        warning("start array is not admissible")
    else{
        ## pathological cases that yield wrong GWLP cannot occur
        ## because of the previous admissibility check
        cat("=== GWLP of start array ===\n")
        print(round(DoE.base::GWLP(countToDmixed(nlev,start)),3))
        cat("=======================================\n")
    }
    }

  ## initial step completed
  ## now optimize shortest words

  kk <- resolution
  if (forMosek){
    ## now extend by current quadratic objective
    ## df(kk+1) + 2 additional variables
    qco1 <- mosek_modelAddLinear(qco1, Us[[kk]])
    qco1 <- mosek_modelAddConeQobj(qco1, df[kk+1])
    qco1 <- mosek_modelAddStart(qco1, sl$sol, add=Hs[[kk]])
       ## extend dimensions of sol object for new cone,
         ## including the current value of the (to be minimized) objective
    #nvar <- as.integer(ncol(qco1$A))   #update nvar
  }
  else{
    ## this is for general MPS that does not know conic constraints
   qco1$Q <- Hs[[kk]]
  }
      bound <- sum(lowerbounds(nruns, nlev, kk)) + 0.5
      qco1$dparam$LOWER_OBJ_CUT <- max(qco1$dparam$LOWER_OBJ_CUT, bound)
   return(qco1)
}

