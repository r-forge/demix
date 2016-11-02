###
### $Id$
###


##-----------------------------------------------------------------------------
## Generic typeless check
is.scalar <- function(x) {
    length(x) == 1L
}


##-----------------------------------------------------------------------------
is.scalar.character <- function(x) {
    is.character(x) && length(x) == 1
}


##-----------------------------------------------------------------------------
is.scalar.numeric <- function(x) {
    is.numeric(x) && length(x) == 1
}


##-----------------------------------------------------------------------------
is.scalar.integer <- function(x) {
    is.integer(x) && length(x) == 1
}


##-----------------------------------------------------------------------------
is.scalar.logical <- function(x) {
    is.logical(x) && length(x) == 1
}

