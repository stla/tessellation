#' Objects imported from other packages
#' @description These objects are imported from other packages.
#'   Follow the links to their documentation:
#'   \code{\link[hash:values]{values}},
#'   \code{\link[hash:keys]{keys}}.
#' @importFrom hash values keys
#' @export values keys
#' @name tessellation-imports
#' @aliases values keys
#' @docType import
NULL


uniqueWith <- function(v, f){
  size <- length(v)
  for(i in seq_len(size-1L)){
    j <- i + 1L
    while(j <= size){
      if(f(v[[i]], v[[j]])){
        v <- v[-j]
        size <- size - 1L
      }else{
        j <- j + 1L
      }
    }
  }
  v[1L:size]
}


sameSegments <- function(seg1ids, seg2ids){
  ((seg1ids[[1L]] == seg2ids[[1L]]) && (seg1ids[[2L]] == seg2ids[[2L]])) ||
    ((seg1ids[[1L]] == seg2ids[[2L]]) && (seg1ids[[2L]] == seg2ids[[1L]]))
}
