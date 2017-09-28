mosek2gurobi <- function(qco, ...){
  if (!"qco" %in% class(qco)) stop("qco must be an object of class qco")
  if (!rev(qco$info$optimizer)[1]=="mosek") stop("wrong optimizer entry for mosek2gurobi")
  aus <- qco
  ## bounds, variable types and starting values for x
  aus$lb <- aus$bx[1,]
  aus$ub <- aus$bx[2,]
  conefixes <- which(aus$lb==0.5 & aus$ub==0.5)
  aus$lb[conefixes] <- 1  ## accomodate different cone coding
  aus$ub[conefixes] <- 1  ## accomodate different cone coding
  aus$bx <- NULL
  if (!is.null(aus$sol)){
    aus$start <- aus$sol$int$xx
    aus$start[conefixes] <- 1
    aus$sol <- NULL
  }
  aus$obj <- aus$c
  aus$c <- NULL
  aus$vtype <- rep("C", length(aus$obj))
  aus$vtype[aus$intsub] <- "I"
  bin <- aus$intsub[aus$lb[aus$intsub]==0 & aus$ub[aus$intsub]==1]
  aus$vtype[bin] <- "B"
  aus$intsub <- NULL

  ## linear constraints
  aus$sense <- rep("=", ncol(aus$bc))
  aus$sense[aus$bc[1,]==-Inf] <- "<="
  aus$sense[aus$bc[2,]==Inf] <- ">="
  aus$rhs <- aus$bc[1,]
  aus$rhs[which(aus$sense=="<=")] <- aus$bc[2,which(aus$sense=="<=")]
  ## all constraints are equality constraints for the problems of interest
  ## unless forced is used with distinct=FALSE ("<=" should not happen)
  aus$bc <- NULL

  ## conic constraints
  ### in Gurobi: a list with each element a list of Qc and rhs
  ### in Mosek: a matrix of lists with each column for one cone
  ###  the two positions with the "-1" element of Gurobi must be the first two
  ###      positions of the Mosek vector, with the lower index the first
  if (!is.null(aus$cones)){
    ## exploits the special form of cone created for the MIParray problem!
    aus$quadcon <- vector(mode="list")
    n <- length(aus$obj)
    for (ii in 1:ncol(aus$cones))
      aus$quadcon[[ii]] <- list(Qc=Matrix::spMatrix(n,n,
                                           i=aus$cones[,ii][[2]][-2],
                                           j=aus$cones[,ii][[2]][-1],
                                           x=c(-1,rep(1, length(aus$cones[,ii][[2]])-2))),
                       rhs=0)
    vectors <- lapply(aus$quadcon, function(curc)
         c(curc$Qc@i[1],curc$Qc@j[1],curc$Qc@i[-1]) + 1)   ## +1, because indexed with 0
    ## sanity check
       if (!length(setdiff(conefixes, sapply(vectors, function(obj) obj[2])))==0)
         stop("too many or wrong positions considered cone fixed ?")
    aus$cones <- NULL
  }

  if (length(aus$info$optimizer)==1) aus$info$optimizer <- rep(aus$info$optimizer, length(aus$info$stati))
  aus$info$optimizer[length(aus$info$optimizer)] <- "gurobi"
  aus
  }
