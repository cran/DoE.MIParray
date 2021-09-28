## print.oa in DoE.base has "\\n" instead of "\n",
## which messes up the output
## once this is fixed in DoE.base, these methods should again be identical
print.oa <- function (x, ...)
{
  xnam <- deparse(substitute(x))
  if (!"oa" %in% class(x))
    stop("this print method is for class oa only")
  attrs <- setdiff(names(attributes(x)), c("origin",
                                           "class", "dim", "dimnames"))
  info <- attr(x, "MIPinfo")$info
  for (a in attrs) attr(x, a) <- NULL
  print.default(x, ...)
  if (!is.null(info$stati)) {
    cat("optimization results:\n")
    print(unlist(info$stati))
  }
  if (length(attrs) > 0) {
    cat("\nfurther attribute(s)", "(accessible with attr(",
        xnam, ", attrname)):", fill = TRUE)
    print(attrs)
  }
}


ff <- function(...){
  ein <- list(...)
  if (!is.numeric(unlist(ein))) stop("ff takes only integers as arguments")
  if (length(ein)==1) ein <- unlist(ein)
  hilf <- expand.grid(rev(lapply(ein, function(obj) 0:(obj-1))))
  hilf <- as.matrix(hilf[,ncol(hilf):1])
  colnames(hilf) <- NULL
  hilf
}

#countToD <- function(nlevels, count){
#  ## function to transform a vector of counts into a design matrix
#  dim <- log(length(count), base=nlevels)
#  if (!round(dim%%1,12)==0) stop("wrong length of count")
#  D <- eval(parse(text=paste("ff(",paste(rep(nlevels,dim),collapse=","),")")))
#  d <- matrix(rep(D[1,], count[1]), ncol=dim, byrow=TRUE)
#  csum <- cumsum(count)
#  for (i in 2:length(count)){
#    from <- csum[i-1]+1
#    ci <- count[i]
#    to <- csum[i]
#    if (ci > 0)
#      d <- rbind(d, matrix(rep(D[i,], ci), ncol=dim, byrow=TRUE))
#  }
#  as.matrix(d)
#}

countToDmixed <- function(nlevels, count){
  ## function to transform a vector of counts into a design matrix
  dim <- length(nlevels)
  if (!prod(nlevels)==length(count)) stop("wrong length of count")
  D <- ff(nlevels)
  d <- matrix(rep(D[1,], count[1]), ncol=dim, byrow=TRUE)
  csum <- cumsum(count)
  for (i in 2:length(count)){
    from <- csum[i-1]+1
    ci <- count[i]
    to <- csum[i]
    if (ci > 0)
      d <- rbind(d, matrix(rep(D[i,], ci), ncol=dim, byrow=TRUE))
  }
  as.matrix(d)
}

dToCount <- function(d, nlevels=NULL, startfrom1=FALSE){
  ## function for creating counts of design rows in lexicographic order
  stopifnot(is.matrix(d))
  stopifnot(all(d >= 0))
  if (startfrom1) {
    stopifnot(all(d >= 1))
    d <- d - 1  ## now d starts from 0
  }
  dim <- ncol(d)
  if(is.null(nlevels)) nlevels <- levels.no(d)
  d <- d[DoE.base::ord(as.data.frame(d)),]
  D <- as.matrix(expand.grid(lapply(rev(nlevels), function(obj) 0:(obj-1))))
  D <- D[,dim:1]
  index <- 1:(prod(nlevels))
  aus <- sapply(index, function(obj)
    sum(apply(d, 1, function(obj2) all(obj2==D[obj,]))))
  if (!sum(aus) == nrow(d)) stop("something went wrong, is the coding of d correct?")
  aus
}

find_perms <- function(A, B){
  ## function for finding the adequate resorting of forced columns
  I <- diag(length(A))
  ps <- combinat::permn(1:length(A))
  for (perm in ps)
    if (all(A[perm]==B)) return(perm)
}

rearrange_forced <- function(forced, levels.orig, levels.new){
  ## forced is count vector in level order levels.orig
  ## levels.orig and levels.new are the level orderings that are needed
  ##      perm is calculated from them with function find_perms
  ## result is count vector in permuted level order levels.new
  hilf <- countToDmixed(levels.orig, forced)
  hilf <- hilf[,find_perms(levels.orig, levels.new)]
  dToCount(hilf, levels.new)
}

levels.no <- function(xx)
{
  ## from DoE.base

  # return the number of levels for design xx
  # apply(xx, 2, function(v) {length(table(v))} )
  # changed table to unique, much faster UG 10 May 13
  # added factor distinction, in order to account for
  #     imbalance by missing levels (UG 03 April 17)
  ff <- FALSE
  if (is.data.frame(xx)){
    if (any(ff <- sapply(xx, is.factor)))
      nflevs <- sapply(xx[ff], nlevels)
  }
  aus <- apply(xx, 2, function(v) length(unique(v)))
  if (any(ff)) aus[ff] <- nflevs
  aus
}

lowerbounds <- function(nruns, nlevels, R){
  ## from DoE.base

  ## the function returns a vector of lower bounds
  ## for nruns^2*A_R of the choose(nfac, R) R-factor sets
  ## the entries are integer-valued
  ## if they are all zero,
  ##   the combination of nruns and nlevels
  ##   permits a design of resolution at least R+1
  nfac <- length(nlevels)
  sets <- nchoosek(nfac, R)
  sapply(1:ncol(sets), function(ii){
    pr <- prod(nlevels[sets[,ii]])
    r <- nruns%%pr
    r*(pr-r)})
}

"nchoosek" <-
  function (n, k)
  {
    # from DoE.base

    # function for calculating all subsets of size k from n objects
    # taken from package vsn (provided by Wolfgang Huber under LGPL)
    # slightly modified to also work for n=k=2 by Ulrike Grömping
    if (!is.numeric(n) || !is.numeric(k) || is.na(n) || is.na(k) ||
        length(n) != 1 || length(k) != 1)
      stop("arguments must be non-NA numeric scalars.")
    if (k > n || k < 0)
      stop("Arguments must satisfy 0 <= k <= n.")
    nck = choose(n, k)
    res = matrix(NA, nrow = k, ncol = nck)
    res[, 1] = 1:k
    j = 2
    repeat {
      if (j > nck)
        break
      res[, j] = res[, j - 1]
      i = k
      repeat {
        res[i, j] = res[i, j] + 1
        if (res[i, j] <= n - (k - i))
          break
        i = i - 1
        stopifnot(i >= 1)
      }
      if (i < k)
        res[(i + 1):k, j] = res[i, j] + 1:(k - i)
      j = j + 1
    }
    stopifnot(all(res[, nck] == (n - k + 1):n))
    stopifnot(all(res <= n) && all(res >= 1))
    return(res)
  }
