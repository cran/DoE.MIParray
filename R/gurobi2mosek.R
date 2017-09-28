gurobi2mosek <- function(qco, ...){
  if (!"qco" %in% class(qco)) stop("qco must be an object of class qco")
  if (!rev(qco$info$optimizer)[1]=="gurobi") stop("wrong optimizer entry for gurobi2mosek")
  aus <- qco
  ## bounds, variable types and starting values for x
  aus$bx <- rbind(aus$lb, aus$ub)
  conefixes <- which(aus$lb==1 & aus$ub==1)
  aus$bx[,conefixes] <- 0.5  ## accomodate different cone coding
  aus$lb <- NULL
  aus$ub <- NULL
  if (!is.null(aus$start)){
    aus$sol$int$xx <- aus$start
    aus$sol$int$xx[conefixes] <- 0.5
    aus$start <- NULL
  }
  aus$c <- aus$obj
  aus$obj <- NULL
  aus$intsub <- which(aus$vtype %in% c("B","I"))
  aus$bx[2,aus$vtype=="B"] <- 1  ## binary
  aus$vtype <- NULL

  ## linear constraints
  aus$bc <- rbind(aus$rhs,aus$rhs)
  if (!all(aus$sense=="=")){
    aus$bc[1, aus$sense=="<="] <- -Inf
    aus$bc[2, aus$sense==">="] <- Inf
  }
  ## all constraints are equality constraints for the problems of interest,
  ## unless force is used with distinct=FALSE (i.e. the case "<=" should not happen)
  aus$rhs <- NULL

  ## conic constraints
  ### in Gurobi: a list with each element a list of Qc and rhs
  ### in Mosek: a matrix of lists with each column for one cone
  ###  the two positions with the "-1" element of Gurobi must be the first two
  ###      positions of the Mosek vector, with the lower index the first
  if (!is.null(aus$quadcon)){
    ## exploits the special form of cone created for the MIParray problem!
    aus$cones <- matrix(list(), nrow=2, ncol=length(aus$quadcon))
    vectors <- lapply(aus$quadcon, function(curc)
         c(curc$Qc@i[1],curc$Qc@j[1],curc$Qc@i[-1]) + 1) ## uses 0-based indices internally
    ## sanity check
       if (!length(setdiff(conefixes, sapply(vectors, function(obj) obj[2])))==0)
         stop("too many or wrong positions considered cone fixed ?")
    for (i in 1:length(aus$quadcon)){
      aus$cones[,i] <- list(type="RQuad", sub=vectors[[i]])
    }
    aus$quadcon <- NULL
  }

  aus$sense <- "min"
  if (length(aus$info$optimizer)==1) aus$info$optimizer <- rep(aus$info$optimizer, length(aus$info$stati))
  aus$info$optimizer[length(aus$info$optimizer)] <- "mosek"
  aus
  }
