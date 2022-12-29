#' @export
mylinsysolve <- function(R, r){
  q <- backsolve(R, forwardsolve(t(R), r))
  return(q)
}

