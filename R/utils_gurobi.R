## not used any more
#gurobi_rsave <- function(mod, sol, file="gurobi_lastsave.rda"){
  ## model and solution must be compatible (user responsibility)
#  if (is.numeric(sol)) if (!length(sol)==ncol(mod$A)) stop("model and sol incompatible")
#  if (is.list(sol))
#    if (!length(sol$x)==ncol(mod$A)) stop("model and sol incompatible") else sol <- sol$x
#    mod$start <- sol
#    saved_gurobi_model <- mod
#    save(saved_gurobi_model, file=file)
#}

#gurobi_rretrieve <- function(file="gurobi_lastsave.rda"){
  ## the file should contain a gurobi model saved with gurobi_rsave
  #if (exists("saved_gurobi_model") && !replace)
  #  stop("existing object saved_gurobi_model not overwritten (replace=FALSE)")
#  load(file=file)
#  saved_gurobi_model
#}

gurobi_modelLastQuadconToLinear <- function(qco1){
  ## function to remove the last quadcon element with all implications
  ## so that the cone constraint becomes a linear constraint to 0

  ## to be used for maintaining an optimum 0 that was found underway

  if (!is.null(qco1$info$conecur))
    if(!qco1$info$conecur) return(qco1)  ## current constraint is already linear

  if (!is.null(qco1$quadcon)){
    ncone <- length(qco1$quadcon)
    if (ncone==0) {
       ## unlikely to occur
       qco1$quadcon <- NULL
       return(qco1)
    }
    lc <- qco1$quadcon[[ncone]]  ## the cone to be removed
    if (ncone==1) qco1$quadcon <- NULL else qco1$quadcon <- qco1$quadcon[-ncone]
    remcols <- which(colSums(abs(as.matrix(lc$Qc))) > 0)
    qco1$A <- qco1$A[,-remcols]
    nc <- ncol(qco1$A)
    qco1$A <- as(Matrix::Matrix(qco1$A) ,'CsparseMatrix')
    qco1$lb <- qco1$lb[-remcols]
    qco1$ub <- qco1$ub[-remcols]
    qco1$vtype <- qco1$vtype[-remcols]
    qco1$obj <- qco1$obj[-remcols]
    if (!is.null(qco1$start)) qco1$start <- qco1$start[-remcols]
    ## existing quadcon elements have to be modified
    ## because of the modified column number
    if (!is.null(qco1$quadcon))
    for (i in 1:length(qco1$quadcon)) {
      qco1$quadcon[[i]]$Qc@Dim <- c(nc,nc)
    }
    qco1$info$conecur <- FALSE
  }
  qco1
}

gurobi_modelAddLinear <- function(qco1, A, b=0, sense="="){
  ## adds a linear constraint with matrix A applied to first columns
  ## b can be a scalar or vector
  ## sense can be a string or vector of strings
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
  nc <- ncol(qco1$A)
  qco1$A <- as(Matrix::Matrix(qco1$A) ,'CsparseMatrix')
  qco1$rhs <- c(qco1$rhs, b)
  if (length(unique(c(qco1$sense,sense)))==1)
    qco1$sense <- qco1$sense[1]
  else{
    if (length(qco1$sense)==1) qco1$sense <- rep(qco1$sense, nr)
    if (length(sense)==1) sense <- rep(sense, nrow(A))
    qco1$sense <- c(qco1$sense, sense)
  }
  qco1$info$conecur<-FALSE
  qco1
}

gurobi_modelAddConeQobj <- function(qco1, dimsquares, rhs=0){
  ## qco1 a gurobi model object
  ## dimsquares the number of variables involved in the squares part of the
  ##    cone specification
  ## function adds 2+dimsquares additional variables
  ## the last dimsquares of which are identical to the last dimsquares linear constraints
  ## the first variable is the square to be optimized
  ##       (and will be constrained to the optimum later),
  ## the second variable is constrained to 1
  ## any existing start vector is removed
  nc <- ncol(qco1$A)
  nr <- nrow(qco1$A)
  if (dimsquares > nr) stop("invalid dimsquares")
  qco1$A <- cbind(qco1$A,
            rbind(matrix(0, nr-dimsquares, 2 + dimsquares),
            cbind(matrix(0,dimsquares, 2), -diag(dimsquares)))
          )
  qco1$A <- as(Matrix::Matrix(qco1$A) ,'CsparseMatrix')
  ncneu <- ncol(qco1$A)
  if (is.null(qco1$quadcon))
    qco1$quadcon <- vector(mode="list")
  else
    for (i in 1:length(qco1$quadcon)) {
      qco1$quadcon[[i]]$Qc@Dim <- c(ncneu,ncneu)
    }

  qco1$quadcon[[length(qco1$quadcon)+1]] <- list(
    Qc = Matrix::spMatrix(ncneu, ncneu,
                c(nc+1, (nc+3):ncneu),
                (nc+2):ncneu,
                c(-1.0, rep(1.0,dimsquares))),
    rhs=rhs)
  qco1$info$conecur <- TRUE
  qco1$lb <- c(qco1$lb, 0, 1, rep(-Inf, dimsquares))
  qco1$ub <- c(qco1$ub, Inf, 1, rep(Inf, dimsquares))
  qco1$obj <- c(0*qco1$obj, 1, rep(0, 1+dimsquares))
  qco1$vtype <- c(qco1$vtype, c("I",rep("C", 1+dimsquares)))
  qco1$start <- NULL
  qco1
}
