write_MPSILPlist <- function (prefix, qcolist, toc=TRUE){
  if (!is.character(prefix)) stop("prefix must be a character string \nthat provides the storage location for the MPS files, \nexcept for the suffix")

    stopifnot(is.list(qcolist))
    stopifnot(all(sapply(qcolist, is.list)))
    stopifnot("qco" %in% class(qcolist[[1]]$ILP))
    if (length(setdiff(c("nlevels", "ILP"), names(qcolist[[1]]) ))>0)
      stop("Entries of qcolist do not contain all necessary elements.")
    nc <- nchar(as.character(length(qcolist)))
    fnnum <- sprintf(paste0("%0", nc, "i"), 1:length(qcolist))

    for (i in 1:length(qcolist)){
      fn <- paste0(prefix, fnnum[i], ".mps")
    message(fn)
      write_MPSILP(file=fn, 
                      qco=qcolist[[i]]$ILP, 
                      qcoinfo=qcolist[[i]]$info)
    }
    if (toc){
      ## write a txt file with an nlevels column for each file
      fn <- paste0(prefix, "_toc.txt")
      message(fn)
      names(qcolist) <- paste0(prefix, fnnum, ".mps")
      if (!is.null(qcolist[[1]]$info)){
      tab <- sapply(qcolist, function(obj) c(obj$info$nruns, obj$info$reso, obj$info$nlev))
      rownames(tab) <- c("nruns", "resotarget", paste0("nlev",1:(nrow(tab)-2)))
      }
      else{
        tab <- sapply(qcolist, function(obj) c(obj$nlevels))
        rownames(tab) <- paste0("nlev",1:(nrow(tab)-2))
      }
      write.table(t(tab), file=fn)
      invisible(t(tab))
    }
}
