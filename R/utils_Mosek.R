## not used any more
mosek_rsave <- function(mod, sol, file="mosek_lastsave.rda"){
  ## model and solution must be compatible (user responsibility)
  #if (is.list(sol)){
  #  if (!length(sol$int$xx)==ncol(mod$A)) stop("model and sol incompatible")
  #}
  #else warning("sol not compatible to mod, only mod is saved")
  saved_mosek_model <- list(mod=mod, sol=sol)
  save(saved_mosek_model, file=file)
}

#mosek_rretrieve <- function(file="mosek_lastsave.rda"){
  ## model and solution must be compatible (user responsibility)
  #if (exists("saved_mosek_model") && !replace)
  #  stop("existing object saved_mosek_model not overwritten (replace=FALSE)")
#  load(file=file)
#  saved_mosek_model
#}

mosek_modelLastQuadconToLinear <- function(qco1){
  ## function to remove the last cone element with all implications
  ## so that the cone constraint becomes a linear constraint to 0

  ## to be used for maintaining an optimum 0 that was found underway

  if (!is.null(qco1$cones)){
    if(!qco1$info$conecur) return(qco1)  ## current constraint is already linear
    ncone <- ncol(qco1$cones)
    if (ncone==0) {
       ## unlikely to occur
       qco1$cones <- NULL
       return(qco1)
    }
    lc <- qco1$cones[,ncone]  ## the cone to be removed
    if (ncone==1)
      qco1$cones <- NULL
    else
      qco1$cones <- qco1$cones[,-ncone, drop=FALSE]

    remcols <- lc[[2]]  ## the columns to be removed
    qco1$A <- qco1$A[,-remcols, drop=FALSE]
    nc <- ncol(qco1$A)
    qco1$A <- as(Matrix::Matrix(qco1$A) ,'CsparseMatrix')
    qco1$bx <- qco1$bx[,-remcols, drop=FALSE]
    qco1$c <- qco1$c[-remcols]
    qco1$intsub <- setdiff(qco1$intsub, remcols)
    if (!is.null(qco1$sol))
      qco1 <- mosek_modelAddStart(qco1, qco1$sol, remove=length(remcols))
    qco1$info$conecur <- FALSE
    #qco1$sol$int$xx[intsub] <- round(qco1$sol$int$xx[intsub])
  }
    qco1
}

mosek_modelAddLinear <- function(qco1, A, b=0, sense="="){
  # input in gurobi notation (Ax sense b)
  ## adds a linear constraint with matrix A applied to first columns
  ## b can be a scalar or vector
  ## sense can be a string or vector of strings with elements from =, <= or >=
  ## existing quadcon elements have to be changed because of the dimension change
  nc <- ncol(qco1$A)
  nr <- nrow(qco1$A)
  if (nc < ncol(A)) stop ("A has too many columns")
  if (length(b)==1) b <- rep(b, nrow(A))
  if (!length(b)==nrow(A))
    stop ("wrong number of elements in b")
  if (!(length(sense)==1 || length(sense)==nrow(A)))
    stop ("wrong number of elements in sense")
  if (!is.numeric(A)) stop("A must be a numeric matrix")
  if (!is.matrix(A)) stop("A must be a matrix")
  if (!is.numeric(b)) stop("b must be numeric")
  if (!(all(sense %in% c("=","<=",">="))))
      stop("wrong elements in sense")
  qco1$A <- rbind(qco1$A, cbind(A, matrix(0,nrow(A), nc-ncol(A))))
  qco1$A <- as(Matrix::Matrix(qco1$A) ,'CsparseMatrix')
  hilf <- rbind(rep(-Inf,nrow(A)), rep(Inf, nrow(A)))
  hilf[1, sense %in% c("=",">=")] <- b[sense %in% c("=",">=")]
  hilf[2, sense %in% c("=","<=")] <- b[sense %in% c("=","<=")]
  qco1$bc <- cbind(qco1$bc, hilf)
  qco1$info$conecur <- FALSE   ## no cone for the added linear constraint
  qco1
}

mosek_modelAddConeQobj <- function(qco1, dimsquares, rhs=0){
  ## qco1 a mosek model object
  ## dimsquares the number of variables involved in the squares part of the
  ##    cone specification
  ## function adds 2+dimsquares additional variables
  ## the last dimsquares of which are identical to the last dimsquares linear constraints
  ## the first variable is the square to be optimized
  ##       (and will be constrained to the optimum later),
  ## the second variable is constrained to 1
  ## any existing start vector is no longer removed
  nc <- ncol(qco1$A)
  nr <- nrow(qco1$A)
  if (dimsquares > nr) stop("invalid dimsquares")
  generated.integers <- numeric(0)
  if (dimsquares > 0)
    generated.integers <- which(apply(qco1$A[(nr-dimsquares+1):nr,,drop=FALSE], 1, function(obj) all(round(obj,8)%%1==0))) + nc + 2
               ## added variables with TRUE will be integers, as will be the first position
  qco1$A <- cbind(qco1$A,
            rbind(matrix(0, nr-dimsquares, 2 + dimsquares),
            cbind(matrix(0,dimsquares, 2), -diag(dimsquares)))
          )
  qco1$A <- as(Matrix::Matrix(qco1$A) ,'CsparseMatrix')
  ncneu <- ncol(qco1$A)
  if (is.null(qco1$cones))
    qco1$cones <- matrix(vector(mode="list"),2,0)

  qco1$cones <- cbind(qco1$cones, list(type="RQUAD",sub=(nc+1):ncneu))
  qco1$info$conecur <- TRUE
  qco1$bx <- cbind(qco1$bx, c(0,Inf), c(0.5,0.5), rbind(rep(-Inf, dimsquares),rep(Inf,dimsquares)))
  qco1$c <- c(0*qco1$c, 1, rep(0, 1+dimsquares))
  qco1$intsub <- c(qco1$intsub, nc+1, generated.integers)  ## add target variable and generated integers

#qco1$sol <- NULL
  qco1
}

mosek_modelAddStart <- function(qco1, sl, add=NULL, remove=0){
  ## add a starting value to a model

  ## the solution to an earlier run of a similar model is available in sl

  ## if qco1 has as many variables than present in sl,
  ## the solution is simply added to qco1
  ##    after rounding the integer elements
  ##    (which is necessary for mosek to accept the solution)
  ## otherwise, it is modified according to add or remove

  ## sl is a solution object
  ## add is the matrix for the current quadratic objective that has been added lately
  ## NxN matrix, psd, integer entries, stored in list Hs
  ## remove is a number of variables to be removed from the solution vector

  ## contains a lot of input checks;
  ## if these consume too much time, reduce them after maturing the code

  if (length(sl$int$xx )==ncol(qco1$A)){
    opt <- sl$int$xx
    opt[qco1$intsub] <- round(opt[qco1$intsub])
           ## integrality necessary for acceptance of starting value
    qco1$sol <- list(int=list(xx=opt))
    return(qco1)
  }
  ## else (otherwise function has terminated already)

  if (is.null(add) && remove==0)
    stop("wrong length sl for neither providing add nor remove")
  if (sum(c(!is.null(add), !remove==0))>1)
    stop("conflicting requests")

  ## different cases
  if (!is.null(add)){
    Q <- add
    if (!is.matrix(Q)) stop("add must be a matrix")
    if (!nrow(Q)==ncol(Q)) stop("add must be a square matrix")
    N <- nQ <- nrow(Q)

    nc <- ncol(qco1$cones)
    coneinds <- qco1$cones[,nc][[2]]  ## the variables involved in last-added cone
    nr <- nrow(qco1$A)                ## number of constraints
    U <- qco1$A[(nr-length(coneinds)+3):nr,1:nQ,drop=FALSE]
    nr <- nrow(U)                     ## number of rows of current U (df)

    optpast <- round(sl$int$xx[1:N])  ## the N original integer (binary) variables

    ## the solution is prepared as initial solution
    opt <- c(sl$int$xx, round(optpast%*%Q%*%optpast), 0.5, as.matrix(U)%*%optpast)
  }
  if (!remove == 0){
    if (!is.numeric(remove)) stop("remove must be numeric")
    if (!remove%%1 == 0) stop("remove must be integer")
    if (!remove > 0) stop("remove must be positive")
    keep <- length(sl$int$xx) - remove
    opt <- sl$int$xx[1:keep]
  }
  opt[qco1$intsub] <- round(opt[qco1$intsub])
  ## integrality necessary for acceptance of starting value    ## value of next target under current starting value
  qco1$sol <- list(int=list(xx=opt))
  qco1
}
