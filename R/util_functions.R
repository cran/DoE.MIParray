print.oa <- DoE.base:::print.oa

ff <- function(...){
  ein <- list(...)
  if (!is.numeric(unlist(ein))) stop("ff takes only integers as arguments")
  if (length(ein)==1) ein <- unlist(ein)
  hilf <- expand.grid(rev(lapply(ein, function(obj) 0:(obj-1))))
  as.matrix(hilf[,ncol(hilf):1])
}

countToD <- function(nlevels, count){
  ## function to transform a vector of counts into a design matrix
  dim <- log(length(count), base=nlevels)
  if (!round(dim%%1,12)==0) stop("wrong length of count")
  D <- eval(parse(text=paste("ff(",paste(rep(nlevels,dim),collapse=","),")")))
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

dToCount <- function(d){
  ## function for creating counts of design rows in lexicographic order
  dim <- ncol(d)
  nlevels <- DoE.base:::levels.no(d)
  d <- d[DoE.base::ord(as.data.frame(d)),]
  D <- as.matrix(expand.grid(lapply(rev(nlevels), function(obj) 0:(obj-1))))
  D <- D[,dim:1]
  index <- 1:(prod(nlevels))
  sapply(index, function(obj)
    sum(apply(d, 1, function(obj2) all(obj2==D[obj,]))))
}
