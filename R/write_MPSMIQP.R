write_MPSMIQP <- function (prefix, nruns, nlevels, resolution=3, distinct=TRUE, start=NULL, forced=NULL, name = "ImproveAR", commentline="* quadratic problem"){
  if (!is.character(prefix))
    stop("prefix must be a character string that provides the storage location for the MPS file (without the suffix .mps)")
  file <- paste(prefix, "mps", sep=".")
  filestart <- paste(prefix, "start", sep=".")
  qco <- create_MIQP(nruns=nruns, nlevels=nlevels, resolution=resolution,
                     distinct=distinct, start=start, forced=forced)

  if (!is.null(qco)){
    ## cvec, bvec, Amat, Qmat, boundmat and commentline are ignored
    ## qco is a qco object produced by create_MIQP
    stopifnot("qco" %in% class(qco))
    ## a single problem
    cvec <- qco$c
    ## needed for incorporating forced runs
    ## "G" occurs for distinct=FALSE
    bsense <- rep("E", ncol(qco$bc))
    bsense[qco$bc[2,]>qco$bc[1,]] <- "G"
    bvec <- qco$bc[1,]
    Amat <- qco$A
    Qmat <- qco$Q
    boundmat <- qco$bx
  }

 ## adapted from package linprog
  nCon <- length(bvec)
  nVar <- length(cvec)

  stopifnot(is.numeric(cvec), is.numeric(bvec))
  stopifnot(is.numeric(as.matrix(Amat)))
  if (!is.null(boundmat)){
    stopifnot(is.matrix(boundmat))
    if (!nrow(boundmat)==2) stop("boundmat must have two rows")
    stopifnot(ncol(boundmat)==nVar)
  }
  else boundmat <- matrix(c(0,1), 2, nVar, byrow=FALSE)

  if (!is.null(bsense)){
    stopifnot(is.character(bsense))
    stopifnot(length(bsense)==nCon)
    stopifnot(all(bsense %in% c("E","G","L","N")))
  }


  if (is.null(names(bvec))) {
    blab <- paste0("R_", 1:nCon)
  }
  else {
    blab <- names(bvec)
    blab <- substr(gsub(" ", "", blab), 1, 8)
    j <- 2
    if (length(unique(blab)) < nCon){
      for (i in 1:nCon)
        while (i > 1 & blab[i] %in% blab[1:(i - 1)]) {
          blab[i] <- paste(substr(blab[i], 1, 7 - nchar(as.character(j))),
                           "_", as.character(j), sep = "")
          j <- j + 1
        }
    }
  }
  if (is.null(names(cvec))) {
    clab <- rep("", nVar)
    clab <- paste0("C_", 1:nVar)
  }
  else {
    clab <- gsub(" ", "", names(cvec))
    clab <- substr(clab, 1,8)
    j <- 2
    if (length(unique(clab)) < nVar){
      for (i in 1:nVar)
        while (i > 1 & clab[i] %in% clab[1:(i - 1)]) {
          clab[i] <- paste(substr(clab[i], 1, 7 - nchar(as.character(j))),
                           "_", as.character(j), sep = "")
          j <- j + 1
        }
    }
  }
  nc <- nchar(clab) ## for column adjustment
  cv <- as.character(signif(cvec, 10))
  ncv <- nchar(cv)
  nb <- nchar(blab) ## for column adjustment
  bv <- as.character(signif(bvec, 10))
  nbv <- nchar(bv)
  Amatv <- matrix(as.character(signif(Amat, 10)), nCon, nVar)
  Qmatv <- matrix(as.character(signif(Qmat, 10)), nVar, nVar)
  ncA <- nchar(Amatv)
  ncQ <- nchar(Qmatv)
  message("start writing ...")
  cat(paste("NAME          ", name, sep = ""), file=file, fill=TRUE)
  cat(commentline, file=file, append=TRUE, fill=TRUE)
  cat("OBJSENSE", file=file, append=TRUE, fill=TRUE)
  cat("    MIN", file=file, append=TRUE, fill=TRUE)
  message("writing ROWS ...")
  cat("ROWS", file=file, append = TRUE, fill=TRUE)
  cat(" N  obj", file=file, append = TRUE, fill=TRUE)
  for (i in 1:nCon) {
    cat(paste0(" ", bsense[i], "  ", blab[i]), file=file,
          append = TRUE, fill=TRUE)
  }
  message("writing COLUMNS ...")
  cat("COLUMNS", file=file, append = TRUE, fill=TRUE)
  line <- "    COUNTS    'MARKER'                 'INTORG'"
  cat(line, file=file, append = TRUE, fill=TRUE)
  for (i in 1:nVar) {
    line <- paste("    ", clab[i], sep = "")
    line <- paste(line, paste(rep(" ", 10 - nc[i]),
                              collapse = ""), "obj", sep = "")
    ## -3 for nchar("obj")
    line <- paste(line, paste(rep(" ", 22 - 3 - ncv[i]), collapse = ""),
                  cv[i], sep = "")
    cat(line, file=file, append = TRUE, fill=TRUE)

    for (j in 1:nCon) {
      if (Amat[j, i] != 0) {
        line <- paste0("    ", clab[i], sep = "")
        line <- paste0(line, paste(rep(" ", 10 - nc[i]), collapse = ""), blab[j],
                       sep = "")
        line <- paste0(line, paste(rep(" ", 22 -
                                         nb[j] - ncA[j,i]), collapse = ""), Amatv[j, i], sep = "")
        cat(line, file=file, append = TRUE, fill=TRUE)
      }
    }
  }
  line <- "    COUNTS    'MARKER'                 'INTEND'"
  cat(line, file=file, append = TRUE, fill=TRUE)
  ## entire COLUMNS group enclosed in markers

  message("writing RHS ...")
  cat("RHS", file=file, append = TRUE, fill=TRUE)
  for (i in 1:nCon) {
    line <- paste("    RHS       ", blab[i], sep = "")
    line <- paste(line, paste(rep(" ", 22 - nb[i] -
                                    nbv[i]), collapse = ""),
                  bv[i], sep = "")
    cat(line, file=file, append = TRUE, fill=TRUE)
  }
  ### defaults for upper integer bounds (UI) are solver-specific (e.g. 1 for gurobi)
  ### Therefore, bounds are always written.
  ### Infinity bounds are currently replaced by nruns. Is that wise?
  message("writing BOUNDS ...")
  cat("BOUNDS", file=file, append = TRUE, fill=TRUE)
  for (i in 1:nVar){
    cat(paste0(" ", ifelse(boundmat[2,i]<Inf, "BV", "PL"), " BND       ",
                 clab[i]), file=file, append=TRUE, fill=TRUE)
  }
  message("writing QUADOBJ ...")
  ## write non-zero upper diagonal elements of qco$Q
  cat("QUADOBJ", file=file, append = TRUE, fill=TRUE)
  for (i in 1:nVar) {
    for (j in 1:i) {
      if (Qmat[j, i] != 0) {
        line <- paste0("    ", clab[i], sep = "")
        line <- paste0(line, paste(rep(" ", 10 - nc[i]), collapse = ""), clab[j],
                       sep = "")
        line <- paste0(line, paste(rep(" ", 22 -
                                         nc[j] - ncQ[j,i]), collapse = ""), Qmatv[j, i], sep = "")
        cat(line, file=file, append = TRUE, fill=TRUE)
      }
    }
  }
  cat("ENDATA", file=file, append = TRUE, fill=TRUE)
  if (!is.null(start)){
    message("writing start value file ...")
    write.table(qco$info$start, file=filestart, row.names=FALSE, col.names=FALSE)
    invisible(start)
  }
}
