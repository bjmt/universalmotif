check_filter_args <- function(x) {
  for (i in which(is.na(x))) stop("'", names(x[i]),
                                  "' must be NULL or numeric", call. = FALSE)
  for (i in which(!vapply(x, is.null, logical(1)))) {
    if (!is.numeric(x[[i]]) || length(x[[i]]) != 1) {
      stop("'", names(x[i]), "' must be NULL or numeric of length 1",
           call. = FALSE)
    }
  }
}

check_logi_args <- function(x) {
  for (i in which(is.na(x))) stop("'", names(x[i]), "' must be logical",
                                  call. = FALSE)
  for (i in which(!vapply(x, is.logical, logical(1)))) {
    stop("'", names(x[i]), "' must be logical", call. = FALSE)
  }
  for (i in which(!vapply(x, function(x) length(x) == 1, logical(1)))) {
    stop("'", names(x[i]), "' must be a logical of length 1", call. = FALSE)
  }
}
