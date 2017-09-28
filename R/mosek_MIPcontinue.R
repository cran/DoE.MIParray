mosek_MIPcontinue <- function(qco, improve=TRUE,
              maxtime=Inf, nthread = 2, mosek.opts=list(verbose=10, soldetail=1),
              mosek.params=list(dparam=list(LOWER_OBJ_CUT=0.5,
                                           MIO_TOL_ABS_GAP=0.2,
                                           INTPNT_CO_TOL_PFEAS=0.00001,
                                           INTPNT_CO_TOL_INFEAS=0.0000001),
                                iparam=list(PRESOLVE_LINDEP_USE="OFF", LOG_MIO_FREQ=100)
                                            )){
  aufruf <- sys.call()

  h <- function(w)
    if( any( grepl( "printing of extremely long output is truncated", w) ) )
      invokeRestart( "muffleWarning" )
  # withCallingHandlers( f(5), warning = h )
  ## around mosek calls will suppress annoying warnings

  ## check inputs
  if (!("qco" %in% class(qco) || "oa" %in% class(qco))) stop("qco must be of class qco or oa")
  if (!"qco" %in% class(qco)) {
    detailed <- 0
    if (!is.null(attr(qco, "history"))){
      detailed <- 1
      history <- attr(qco, "history")
      if (!is.null(attr(qco, "matrices"))){
        matrices <- attr(qco, "matrices")
        detailed <- 3
      }
    }
    qco <- attr(qco, "MIPinfo")
  }
  if (is.null(qco)) stop("invalid object qco")
  if (!is.logical(improve)) stop("improve must be logical")
  nruns <- qco$info$nruns
  nlev <- qco$info$nlev
  nfac <- length(nlev)
  sl <- qco$sl
  last.k <- qco$info$last.k
  if (last.k==nfac && !improve){
    warning("improve was set to TRUE, because last.k=nfac")
    improve <- TRUE
  }

  ## refuse to improve, if already optimal
  if (improve && rev(qco$info$stati)[1] %in% c("INTEGER_OPTIMAL", "OPTIMAL"))
    stop("improvement unnecessary, ", "already optimal")
  ## append or replace optimizer info for last optimization step
  ## gurobi2mosek always replaces the last entry
  ## this needs to be changed for improve=FALSE
  if (rev(qco$info$optimizer)[1]=="gurobi"){
    if (improve)
      qco <- gurobi2mosek(qco)
    else{
      hilf <- qco$info$optimizer
      if (length(hilf)==1)
        hilf <- rep(hilf, length(qco$info$stati))
      hilf <- c(hilf, "mosek")
      qco <- gurobi2mosek(qco)
      qco$info$optimizer <- hilf
    }
  }

  if (!is.null(maxtime)){
    if (!is.numeric(maxtime)) stop("maxtime must be numeric")
    if (!length(maxtime)==1) stop("maxtime must be scalar")
    if (!maxtime>0) stop("maxtime must be positive")
  }
  if (!is.numeric(nthread)) stop("nthread must be numeric")
  if (!length(nthread)==1) stop("nthread must be scalar")
  if (!nthread>=0) stop("nthread must be non-negative")

  if (!is.list(mosek.opts)) stop("mosek.opts must be a named list")
  if (!is.list(mosek.params)) stop("mosek.params must be a named list")
  if (!all(names(mosek.params) %in% c("iparam","dparam","sparam")))
      stop("invalid entries in list mosek.params")

   ## the stricter time request wins
  ## the stricter time request wins, LOWER_OBJ_CUT is always set
  if (maxtime<Inf){
    if (!is.null(mosek.params$dparam)){
      mosek.params$dparam$MIO_MAX_TIME <-
        min(maxtime, mosek.params$dparam$MIO_MAX_TIME)
      mosek.params$dparam$LOWER_OBJ_CUT <- max(0.5, mosek.params$dparam$LOWER_OBJ_CUT)
    }
    else mosek.params$dparam <- list(LOWER_OBJ_CUT=0.5, MIO_MAX_TIME=maxtime)
  }
  else{
    if (!is.null(mosek.params$dparam))
      mosek.params$dparam$LOWER_OBJ_CUT <- max(0.5, mosek.params$dparam$LOWER_OBJ_CUT)
    else mosek.params$dparam <- list(LOWER_OBJ_CUT=0.5)
  }

  ## the smaller thread request wins
  if (!is.null(mosek.params$iparam))
    mosek.params$iparam$NUM_THREADS <-
      min(nthread, mosek.params$iparam$NUM_THREADS)
    else mosek.params$iparam <- list(NUM_THREADS=nthread)


  opts <- mosek.opts
    qco$iparam <- mosek.params$iparam
    qco$dparam <- mosek.params$dparam
    qco$sparam <- mosek.params$sparam

    N <- prod(nlev)
    kmax <- last.k + 1

    if (improve){
      ## can start immediately to try further improvement
      ## increase LOWER_OBJ_CUT to current largest bound, if larger
      ##    (add small margin, because everything should be integer)
      if (qco$info$objbound > qco$dparam$LOWER_OBJ_CUT)
        qco$dparam$LOWER_OBJ_CUT <- qco$info$objbound + 0.3

      hilf <- qco$info
      qco$info <- NULL
      withCallingHandlers(
        {sl <- Rmosek::mosek(qco, opts=opts)},
      warning = h )
      qco$info <- hilf
      opt <- sl$sol$int$xx  ## A_last.k improved
    }
    else{
      ## optimize the next variant
      if (qco$info$conecur){
        nc <- ncol(qco$cones)
        poscopt <- qco$cones[,nc][[2]][1]
        vcur <- qco$sol$int$xx[poscopt]
        ## remove conic constraint where linear is sufficient
        ##    (vcur is integral, thus liberal tolerance OK)
        ## else make sure A_last.k cannot deteriorate
        if (round(vcur,2)==0) qco <- mosek_modelLastQuadconToLinear(qco)
        else qco$bx[2, poscopt] <- vcur + 0.3
      }

     ## preparation of matrices needed in formulating the optimization problem
     df1 <- nlev-1
     D <- ff(nlev)
     df <- 1
     for (j in 1:nfac) df <- c(df,sum(apply(matrix(df1[DoE.base:::nchoosek(nfac,j)],nrow=j),2,prod)))
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
        Hs[[kmax]] <- round(tcrossprod(mm[,(dfweg[kmax]+1):dfweg[kmax+1]]),8)
        Us[[kmax]] <- t(mm[,(dfweg[kmax]+1):dfweg[kmax+1]]) ## crossprod is H
      ## we start from the existing model qco and add a cone constraint for the
        ## new variables (Us[[kmax]]%*%x)
      qco <- mosek_modelAddLinear(qco, Us[[kmax]])
      qco <- mosek_modelAddConeQobj(qco, df[kmax + 1])
     # if (rem>0) sl$sol$int$xx <- sl$sol$int$xx[1:(length(sl$sol$int$xx)-rem)]
      qco <- mosek_modelAddStart(qco, qco$sol, add=Hs[[kmax]])

      ## run optimization
      hilf <- qco$info
      qco$info <- NULL
      withCallingHandlers(
        {sl <- Rmosek::mosek(qco, opts=opts)},
      warning = h )
      qco$info <- hilf
      opt <- sl$sol$int$xx     ## A_kmax minimized
      vcur <- round(opt[1:N])%*%Hs[[kmax]]%*%round(opt[1:N])
       ## next optimum
  }
    if (!improve) cat(paste("=== GWLP after optimizing A",kmax," ===\n",sep=""))
    else cat(paste("=== GWLP after further improving A",last.k," ===\n",sep=""))
    print(round(DoE.base::GWLP(countToDmixed(nlev,round(opt[1:N]))),3))
    cat("=================================\n")

  feld <- countToDmixed(nlev, round(opt[1:N])) + 1
  class(feld) <- c("oa", "matrix")
  attr(feld, "origin") <- list(package="DoE.MIParry", call=aufruf)
  nsteps <- length(qco$info$stati)
  if (improve) qco$info$stati[nsteps] <- sl$sol$int$solsta
  else{
    qco$info$stati <- c(qco$info$stati, sl$sol$int$solsta)
    names(qco$info$stati)[nsteps+1] <- paste("A",kmax,sep="")
  }
  opt[qco$intsub] <- round(opt[qco$intsub])
      ## necessary for the solution to be accepted as starting value
  qco$sol <- list(int=list(xx=opt))   ## add latest solution
  if (!improve) qco$info$last.k <- kmax
  qco$info$objval <- sl$sol$int$pobjval  ## current objective value
  if (improve)
     qco$info$objbound <- max(sl$sol$int$pobjbound, qco$info$objbound)
  else qco$info$objbound <- sl$sol$int$pobjbound
            ## for improve, sharpest objective bound seen so far

  if (qco$info$last.k < nfac ||
      !all(qco$info$stati %in% c("INTEGER_OPTIMAL","OPTIMAL","USER_OBJ_LIMIT")))
       attr(feld, "MIPinfo") <- qco
  ## other detail
  if (detailed && improve) {
    history$sols[[length(history$sols)]] <- sl
    ## only solution changed
    ## is it beneficial to also change problem in case of changed optimizer?
    attr(feld, "history") <- history
    if (detailed==3){
      ## unchanged matrices
      attr(feld, "matrices") <- matrices
    }
  }
  if (detailed && !improve) {
    history$probs <- c(history$probs, qco)
    history$sols <- c(history$sols, sl)
    attr(feld, "history") <- history
    if (detailed==3){
      matrices$Us <- c(matrices$Us, Us[[kmax]])
      matrices$Hs <- c(matrices$Hs, Hs[[kmax]])
      attr(feld, "matrices") <- matrices
    }
  }
  Rmosek::mosek_clean()  ## prevent problems that lead to machine instabilities
  feld
}

