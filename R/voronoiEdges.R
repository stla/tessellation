#' @title R6 class representing an edge of a Voronoï cell
#'
#' @description An edge is given by two vertices,
#'   named \code{A} and \code{B}.
#'
#' @export
#' @importFrom R6 R6Class
Edge <- R6Class(

  "Edge",

  private = list(
    .A = c(NA_real_, NA_real_, NA_real_),
    .B = c(NA_real_, NA_real_, NA_real_)
  ),

  active = list(
    #' @field A get or set the vertex \code{A}
    A = function(value) {
      if (missing(value)) {
        private[[".A"]]
      } else {
        A <- as.vector(value)
        stopifnot(
          is.numeric(A),
          length(A) == 3L,
          !any(is.na(A))
        )
        private[[".A"]] <- A
      }
    },

    #' @field B get or set the vertex \code{B}
    B = function(value) {
      if (missing(value)) {
        private[[".B"]]
      } else {
        B <- as.vector(value)
        stopifnot(
          is.numeric(B),
          length(B) == 3L,
          !any(is.na(B))
        )
        private[[".B"]] <- B
      }
    }
  ),

  public = list(
    #' @description Create a new \code{Edge} object.
    #' @param A the vertex \code{A}
    #' @param B the vertex \code{B}
    #' @return A new \code{Edge} object.
    #' @examples edge <- Edge$new(c(1, 1, 1), c(1, 2, 3))
    #' edge
    #' edge$A
    #' edge$A <- c(1, 0, 0)
    #' edge
    initialize = function(A, B) {
      A <- as.vector(A)
      stopifnot(
        is.numeric(A),
        length(A) == 3L,
        !any(is.na(A))
      )
      B <- as.vector(B)
      stopifnot(
        is.numeric(B),
        length(B) == 3L,
        !any(is.na(B))
      )
      if(isTRUE(all.equal(A, B))){
        stop(
          "`A` and `B` must be distinct points."
        )
      }
      private[[".A"]] <- A
      private[[".B"]] <- B
    },

    #' @description Show instance of an \code{Edge} object.
    #' @param ... ignored
    #' @examples Edge$new(c(2, 0, 0), c(3, -1, 4))
    print = function(...) {
      cat("Edge:\n")
      cat(" vertex A: ", toString(private[[".A"]]), "\n", sep = "")
      cat(" vertex B: ", toString(private[[".B"]]), "\n", sep = "")
    }
  )
)

#' @title R6 class representing an infinite edge of a Voronoï cell
#'
#' @description An infinite edge is given by a vertex, its origin,
#'   and a vector, its direction.
#'
#' @export
#' @importFrom R6 R6Class
IEdge <- R6Class(

  "IEdge",

  private = list(
    .O = c(NA_real_, NA_real_, NA_real_),
    .direction = c(NA_real_, NA_real_, NA_real_)
  ),

  active = list(
    #' @field O get or set the vertex \code{O}
    O = function(value) {
      if (missing(value)) {
        private[[".O"]]
      } else {
        O <- as.vector(value)
        stopifnot(
          is.numeric(O),
          length(O) == 3L,
          !any(is.na(O))
        )
        private[[".O"]] <- O
      }
    },

    #' @field direction get or set the vector \code{direction}
    direction = function(value) {
      if (missing(value)) {
        private[[".direction"]]
      } else {
        direction <- as.vector(value)
        stopifnot(
          is.numeric(direction),
          length(direction) == 3L,
          !any(is.na(direction))
        )
        private[[".direction"]] <- direction
      }
    }
  ),

  public = list(
    #' @description Create a new \code{IEdge} object.
    #' @param O the vertex \code{O} (origin)
    #' @param direction the vector \code{direction}
    #' @return A new \code{IEdge} object.
    #' @examples iedge <- IEdge$new(c(1, 1, 1), c(1, 2, 3))
    #' edge
    #' edge$O
    #' edge$O <- c(1, 0, 0)
    #' edge
    initialize = function(O, direction) {
      O <- as.vector(O)
      stopifnot(
        is.numeric(O),
        length(O) == 3L,
        !any(is.na(O))
      )
      direction <- as.vector(direction)
      stopifnot(
        is.numeric(direction),
        length(direction) == 3L,
        !any(is.na(direction))
      )
      if(c(crossprod(direction)) == 0){
        stop(
          "`direction` must be a non-null vector."
        )
      }
      private[[".O"]] <- O
      private[[".direction"]] <- direction
    },

    #' @description Show instance of an \code{IEdge} object.
    #' @param ... ignored
    #' @examples IEdge$new(c(2, 0, 0), c(3, -1, 4))
    print = function(...) {
      cat("IEdge:\n")
      cat("  origin O: ", toString(private[[".O"]]), "\n", sep = "")
      cat(" direction: ", toString(private[[".direction"]]), "\n", sep = "")
    }
  )
)
