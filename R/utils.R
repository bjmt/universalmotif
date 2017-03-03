stringify <- function(x) {
  x <- substitute(x)
  if (class(x) == "name") x <- deparse(x)
  return(x)
}
