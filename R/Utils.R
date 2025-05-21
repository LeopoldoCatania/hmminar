MyGrid <- function(iM, iL) {
  
  expand.grid(replicate(iL, 1:iM, simplify = FALSE))
  
}

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}