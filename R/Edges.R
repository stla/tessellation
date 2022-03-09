#' @title R6 class representing an edge in dimension 3.
#'
#' @description An edge is given by two vertices in the 3D space,
#'   named \code{A} and \code{B}. This is for example an edge of a Voronoï cell
#'   of a 3D Delaunay tessellation.
#'
#' @export
#' @importFrom R6 R6Class
#' @importFrom rgl lines3d cylinder3d shade3d
Edge3 <- R6Class(

  "Edge3",

  private = list(
    .A = c(NA_real_, NA_real_, NA_real_),
    .B = c(NA_real_, NA_real_, NA_real_),
    .idA = NA_integer_,
    .idB = NA_integer_
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
    },

    #' @field idA get or set the id of vertex \code{A}
    idA = function(value) {
      if(missing(value)){
        private[[".idA"]]
      } else {
        idA <- as.vector(value)
        stopifnot(
          is.atomic(idA),
          is.numeric(idA),
          length(idA) == 1L
        )
        private[[".idA"]] <- idA
      }
    },

    #' @field idB get or set the id of vertex \code{B}
    idB = function(value){
      if(missing(value)){
        private[[".idB"]]
      } else {
        idB <- as.vector(value)
        stopifnot(
          is.atomic(idB),
          is.numeric(idB),
          length(idB) == 1L
        )
        private[[".idB"]] <- idB
      }
    }

  ),

  public = list(
    #' @description Create a new \code{Edge3} object.
    #' @param A the vertex \code{A}
    #' @param B the vertex \code{B}
    #' @param idA the id of vertex \code{A}, an integer; can be missing
    #' @param idB the id of vertex \code{B}, an integer; can be missing
    #' @return A new \code{Edge3} object.
    #' @examples edge <- Edge3$new(c(1, 1, 1), c(1, 2, 3))
    #' edge
    #' edge$A
    #' edge$A <- c(1, 0, 0)
    #' edge
    initialize = function(A, B, idA, idB) {
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

      if(!missing(idA)){
        idA <- as.vector(idA)
        stopifnot(
          is.atomic(idA),
          is.numeric(idA),
          length(idA) == 1L
        )
        private[[".idA"]] <- as.integer(idA)
      }else{
        private[[".idA"]] <- NA_integer_
      }
      if(!missing(idB)){
        idB <- as.vector(idB)
        stopifnot(
          is.atomic(idB),
          is.numeric(idB),
          length(idB) == 1L
        )
        private[[".idB"]] <- as.integer(idB)
      }else{
        private[[".idB"]] <- NA_integer_
      }

    },

    #' @description Show instance of an \code{Edge3} object.
    #' @param ... ignored
    #' @examples Edge3$new(c(2, 0, 0), c(3, -1, 4))
    print = function(...) {
      idA <- private[[".idA"]]
      idB <- private[[".idB"]]
      if(!is.na(idA)){
        idA <- sprintf(" (%d)", idA)
      }else{
        idA <- ""
      }
      if(!is.na(idB)){
        idB <- sprintf(" (%d)", idB)
      }else{
        idB <- ""
      }
      cat("Edge:\n")
      cat(" vertex A: ", toString(private[[".A"]]), idA, "\n", sep = "")
      cat(" vertex B: ", toString(private[[".B"]]), idB, "\n", sep = "")
    },

    #' @description Plot an \code{Edge3} object.
    #' @param edgeAsTube Boolean, whether to plot the edge as a tube
    #' @param tubeRadius the radius of the tube
    #' @param tubeColor the color of the tube
    #' @examples library(tessellation)
    #' d <- delaunay(centricCuboctahedron())
    #' v <- voronoi(d)
    #' cell13 <- v[[13]] # the point (0, 0, 0), at the center
    #' isBoundedCell(cell13) # TRUE
    #' library(rgl)
    #' open3d(windowRect = c(50, 50, 562, 562))
    #' invisible(lapply(cell13[["cell"]], function(edge) edge$plot()))
    plot = function(edgeAsTube = FALSE, tubeRadius, tubeColor){
      stopifnot(isBoolean(edgeAsTube))
      if(edgeAsTube){
        edge <- cylinder3d(
          rbind(
            private[[".A"]],
            private[[".B"]]
          ),
          radius = tubeRadius,
          sides = 90
        )
        shade3d(edge, color = tubeColor)
      }else{
        lines3d(rbind(
          private[[".A"]],
          private[[".B"]]
        ))
      }
    },

    #' @description Stack the two vertices of the edge (this is for internal
    #' purpose).
    stack = function(){
      rbind(
        private[[".A"]],
        private[[".B"]]
      )
    }
  )
)


#' @title R6 class representing a semi-infinite edge in dimension 3
#'
#' @description A semi-infinite edge is given by a vertex, its origin,
#'   and a vector, its direction. Voronoï diagrams possibly have such edges.
#'
#' @export
#' @importFrom R6 R6Class
IEdge3 <- R6Class(

  "IEdge3",

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
    #' @description Create a new \code{IEdge3} object.
    #' @param O the vertex \code{O} (origin)
    #' @param direction the vector \code{direction}
    #' @return A new \code{IEdge3} object.
    #' @examples iedge <- IEdge3$new(c(1, 1, 1), c(1, 2, 3))
    #' iedge
    #' iedge$O
    #' iedge$O <- c(1, 0, 0)
    #' iedge
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

    #' @description Show instance of an \code{IEdge3} object.
    #' @param ... ignored
    #' @examples IEdge3$new(c(2, 0, 0), c(3, -1, 4))
    print = function(...) {
      cat("IEdge:\n")
      cat("  origin O: ", toString(private[[".O"]]), "\n", sep = "")
      cat(" direction: ", toString(private[[".direction"]]), "\n", sep = "")
    }
  )
)


#' @title R6 class representing an edge in dimension 2.
#'
#' @description An edge is given by two vertices in the 2D space,
#'   named \code{A} and \code{B}. This is for example an edge of a Voronoï cell
#'   of a 2D Delaunay tessellation.
#'
#' @export
#' @importFrom R6 R6Class
#' @importFrom graphics segments
Edge2 <- R6Class(

  "Edge2",

  private = list(
    .A = c(NA_real_, NA_real_),
    .B = c(NA_real_, NA_real_)
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
          length(A) == 2L,
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
          length(B) == 2L,
          !any(is.na(B))
        )
        private[[".B"]] <- B
      }
    }
  ),

  public = list(
    #' @description Create a new \code{Edge2} object.
    #' @param A the vertex \code{A}
    #' @param B the vertex \code{B}
    #' @return A new \code{Edge2} object.
    #' @examples edge <- Edge2$new(c(1, 1), c(2, 3))
    #' edge
    #' edge$A
    #' edge$A <- c(1, 0)
    #' edge
    initialize = function(A, B) {
      A <- as.vector(A)
      stopifnot(
        is.numeric(A),
        length(A) == 2L,
        !any(is.na(A))
      )
      B <- as.vector(B)
      stopifnot(
        is.numeric(B),
        length(B) == 2L,
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

    #' @description Show instance of an \code{Edge2} object.
    #' @param ... ignored
    #' @examples Edge2$new(c(2, 0), c(3, -1))
    print = function(...) {
      cat("Edge:\n")
      cat(" vertex A: ", toString(private[[".A"]]), "\n", sep = "")
      cat(" vertex B: ", toString(private[[".B"]]), "\n", sep = "")
    },

    #' @description Plot an \code{Edge2} object.
    #' @param color the color of the edge
    #' @param ... graphical parameters such as \code{lty} or \code{lwd}
    #' @importFrom graphics segments
    #' @examples library(tessellation)
    #' centricSquare <- rbind(
    #'   c(-1, 1), c(1, 1), c(1, -1), c(-1, -1), c(0, 0)
    #' )
    #' d <- delaunay(centricSquare)
    #' v <- voronoi(d)
    #' cell5 <- v[[5]] # the cell of the point (0, 0), at the center
    #' isBoundedCell(cell5) # TRUE
    #' plot(centricSquare, type = "n")
    #' invisible(lapply(cell5[["cell"]], function(edge) edge$plot()))
    plot = function(color = "black", ...){
      private[[".A"]] -> A
      private[[".B"]] -> B
      segments(A[1L], A[2L], B[1L], B[2L], col = color, ...)
    },

    #' @description Stack the two vertices of the edge (this is for internal
    #' purpose).
    stack = function(){
      rbind(
        private[[".A"]],
        private[[".B"]]
      )
    }
  )
)


#' @title R6 class representing a semi-infinite edge in dimension 2
#'
#' @description A semi-infinite edge is given by a vertex, its origin,
#'   and a vector, its direction. Voronoï diagrams possibly have such edges.
#'
#' @export
#' @importFrom R6 R6Class
IEdge2 <- R6Class(

  "IEdge2",

  private = list(
    .O = c(NA_real_, NA_real_),
    .direction = c(NA_real_, NA_real_)
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
          length(O) == 2L,
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
          length(direction) == 2L,
          !any(is.na(direction))
        )
        private[[".direction"]] <- direction
      }
    }
  ),

  public = list(
    #' @description Create a new \code{IEdge2} object.
    #' @param O the vertex \code{O} (origin)
    #' @param direction the vector \code{direction}
    #' @return A new \code{IEdge2} object.
    #' @examples iedge <- IEdge2$new(c(1, 1), c(2, 3))
    #' iedge
    #' iedge$O
    #' iedge$O <- c(1, 0)
    #' iedge
    initialize = function(O, direction) {
      O <- as.vector(O)
      stopifnot(
        is.numeric(O),
        length(O) == 2L,
        !any(is.na(O))
      )
      direction <- as.vector(direction)
      stopifnot(
        is.numeric(direction),
        length(direction) == 2L,
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

    #' @description Show instance of an \code{IEdge2} object.
    #' @param ... ignored
    #' @examples IEdge2$new(c(2, 0), c(3, -1))
    print = function(...) {
      cat("IEdge:\n")
      cat("  origin O: ", toString(private[[".O"]]), "\n", sep = "")
      cat(" direction: ", toString(private[[".direction"]]), "\n", sep = "")
    }
  )
)

#' @title Edge2 or Edge3
#' @noRd
newEdge <- function(A, B){
  dimension <- length(A)
  if(dimension == 2L){
    Edge2$new(A, B)
  }else if(dimension == 3L){
    Edge3$new(A, B)
  }
}

#' @title IEdge2 or IEdge3
#' @noRd
newIEdge <- function(O, direction){
  dimension <- length(O)
  if(dimension == 2L){
    IEdge2$new(O, direction)
  }else if(dimension == 3L){
    IEdge3$new(O, direction)
  }
}
