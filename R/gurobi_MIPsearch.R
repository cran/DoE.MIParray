gurobi_MIPsearch <- function(nruns, nlevels, resolution=3, maxtime=60,
              stopearly=TRUE, listout=FALSE, orders=NULL,
              distinct=TRUE, detailed=0,
              forced=NULL, find.only=TRUE, nthread=2,
              heurist=0.5, MIQCPMethod=0, MIPFocus=1,
              gurobi.params = list(BestObjStop = 0.5, OutputFlag=0)){
  aufruf <- sys.call()
  nlev <- nlevels

  ## the function ensures the requested resolution and throws an error, if that is not possible
  ## it is in principle possible to interrupt the process and keep the result, but unstable situations may occur

  ## check inputs
  stopifnot(is.logical(listout))
  stopifnot(is.logical(stopearly))
  stopifnot(is.null(orders) || is.list(orders))
  if (!is.null(orders)) stopifnot(all(sapply(lapply(orders, sort), function(obj) all(obj==nlev))))
  if (!is.numeric(detailed)) stop("detailed must be numeric")
  detailed <- max(0, detailed)
  detailed <- min(3, detailed)
  detailed <- floor(detailed)  ## 0,1,2,3 possible, 1.999 -> 1 (e.g.)
  if (!is.numeric(nruns)) stop("nruns must be an integer number")
  if (!is.numeric(nlev)) stop("nlev must have integer elements")
  nfac <- length(nlev)
  if (nfac < 2) stop("nlev must have at least two elements")
  if (!is.numeric(resolution)) stop("resolution must be an integer number")
  if (!is.null(maxtime)){
    if (!is.numeric(maxtime)) stop("maxtime must be numeric")
    if (!length(maxtime)==1) stop("maxtime must be scalar")
    if (!maxtime>0) stop("maxtime must be positive")
  }
  if (!is.numeric(nthread)) stop("nthread must be numeric")
  if (!length(nthread)==1) stop("nthread must be scalar")
  if (!nthread>=0) stop("nthread must be non-negative")
  if (!is.numeric(heurist)) stop("heurist must be numeric")
  if (!length(heurist)==1) stop("heurist must be scalar")
  if (!(heurist>=0 && heurist<=1)) stop("heurist must be between 0 and 1")
  if (!is.numeric(MIQCPMethod)) stop("MIQPCMethod must be a valid integer number")

  if (!length(nruns)==1) stop("nruns must be scalar")
  if (!length(resolution)==1) stop("resolution must be scalar")
  if (resolution < 2) stop("resolution must be at least 2")
  if (resolution > nfac) stop("resolution must not be larger than the number of factors")
  strength <- resolution - 1
  if (!is.logical(distinct)) stop("distinct must be logical")
    if (distinct && nruns>prod(nlev)) stop("too many runs for design with distinct runs")

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
    if (!sum(forced)<nruns) stop("vector forced fixes all runs")
    if (distinct && !all(forced %in% c(0,1)))
      stop("forced array does comply with option distinct")
  }

     ## preliminary exclusion of definitely infeasible cases
    N <- prod(nlev)
    if (strength > 0)
      if (!DoE.base::oa_feasible(nruns, nlev, strength)) stop("requested array does not exist")

    ## global optimum of n^2*A_R
    globopt <- round(lowerbound_AR(nruns, nlevels, resolution)*nruns^2)

  if (!is.list(gurobi.params)) stop("gurobi.params must be a named list")
  params <- gurobi.params
  if (is.null(params)) params <- vector(mode="list")
  if (is.null(params$BestObjStop)) params$BestObjStop <- 0.5
     ## BestObjStop 0.5 enforced (user can specify larger value, if desired)
     ##      for A_resolution, the bound can and will be made
     ##      using internal function lowerbounds from DoE.base
     ##      this happens in the loop
  if (params$BestObjStop < 0) {
       params$BestObjStop <- 0.5      ## objective values smaller than zero are nonsense
       warning("negative gurobi parameter BestObjStop was replaced by 0.5")
  }
  ## the strictest restriction of maximum time wins
  params$TimeLimit <- min(params$TimeLimit, maxtime)
  if (params$TimeLimit == Inf) params$TimeLimit <- NULL

  ## the largest permissible choice of Heuristics wins
  params$Heuristics <- max(params$Heuristics, heurist)
  if (params$Heuristics > 1) params$Heuristics <- 1
  if (params$Heuristics < 0) params$Heuristics <- 0

  ## the smaller thread request wins
  params$Threads <- min(nthread, params$Threads)

  if (!MIQCPMethod %in% c(-1,0,1)){
    MIQCPMethod <- 0
    warning("invalid specification of MIQCPMethod ignored")
  }
  if (!is.null(params$MIQCPMethod)){
  if (!MIQCPMethod==params$MIQCPMethod){
     warning("conflicting specifications of MIQCPMethod ignored")
     params$MIQCPMethod <- 0  ## default
  }
  }
  else params$MIQCPMethod <- MIQCPMethod

  if (!MIPFocus %in% 0:3){
    MIPFocus <- 0
    warning("invalid specification of MIPFocus ignored")
  }
  if (!is.null(params$MIPFocus)){
    if (!MIPFocus==params$MIPFocus){
      warning("conflicting specifications of MIPFocus ignored")
      params$MIPFocus <- 0  ## default
    }
  }
  else params$MIPFocus <- MIPFocus

  if (is.null(orders)) ll <- unique(combinat::permn(nlevels))
  else ll <- orders

  opt <- Inf
  planopt <- NULL
  iGWLPopt <- NULL
  optorder <- NULL
  message(paste(length(ll), "orders to be examined"))
  i <- 0
  if (listout) liste <- vector(mode="list", length=length(ll))
  for (l in ll){
    i <- i+1
    cat(      "\n***********************************\n")
    cat(paste("*********** order", i, "***********\n"))
    cat(      "***********************************\n\n")
    if (is.null(forced))
      plan <- try(gurobi_MIParray(nruns, l, resolution, maxtime=maxtime,
              distinct=distinct, detailed=detailed,
              find.only=find.only, nthread=nthread,
              gurobi.params=gurobi.params))
    else
      plan <- try(gurobi_MIParray(nruns, l, resolution, maxtime=maxtime,
              distinct=distinct, detailed=detailed,
              forced=rearrange_forced(forced, levels.orig=nlevels, levels.new=l),
              find.only=find.only, nthread=nthread,
              gurobi.params=gurobi.params))

    if (listout) liste[[i]] <- plan
    cur <- Inf
    if (!"try-error" %in% class(plan)){
      suppressWarnings({curGWLP <- round(GWLP(plan)*nruns^2)}) ## integer
      cur <- curGWLP[resolution + 1]

    if (cur < opt){
      opt <- cur
      optorder <- l
      planopt <- plan
      iGWLPopt <- curGWLP
      if (stopearly && !find.only){
        MIPinfo <- attr(plan, "MIPinfo")
        if (!"qco" %in% class(MIPinfo)){
        if (MIPinfo$stati[2]=="OPTIMAL"){
        message("optimum found")
        break
        }
        }
      }
      if (stopearly && globopt==opt){
        message("optimum found")
        break
        }
    }
    else if (cur==opt && opt < Inf && nfac > resolution){
      # resolve ties by GWLP
      for (gg in (resolution+2):(nfac+1)){
         if (curGWLP[gg] > iGWLPopt[gg])
         break
         else if (curGWLP[gg] < iGWLPopt[gg]){
              opt <- cur
              optorder <- l
              planopt <- plan
              iGWLPopt <- curGWLP
              break
         }
      }
    }
  }
  }
  aus <- planopt
  attr(aus, "optorder") <- optorder
  if (listout) {
    attr(aus, "orders") <- ll
    attr(aus, "allplans") <- liste
  }
  aus
}
