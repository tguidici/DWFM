#' simulated data for examples
#'
#' A list of lists containing simulated data from a 3x1 experimental design (\code{X[[1]]}), and a 2x2 experimental design (\code{X[[2]]})
#'
#' @format A list of length 2, containing lists of length 3 and 4, respectively
#' \describe{
#'     \item{X[[1]]}{data from a 3x1 experimental design (\code{X[[1]][[1]]} through \code{X[[1]][[3]]}.}
#'     \item{X[[2]]}{data from a 2x2 experimental design (\code{X[[2]][[1]]} through \code{X[[2]][[4]]}.}
#' }
"X"

#' Q used to generate simulated data
#'
#' true Q values used in data generation
#'
#' @format A matrix of size c(30,2)
"Q_vals"

#' B^k values used to generate simulated data
#'
#'  a list of length 2, containing matrices of sizes c(30,3) and c(30,4), respectively, containing true values used in data generation
#'
#' @format A list
#' \describe{
#'     \item{B_vals[[1]]}{a matrix of dimensions c(30,3), normalized with (ID1).}
#'      \item{B_vals[[2]]}{a matrix of dimensions c(30,4), normalized with (ID2).}
#' }
"B_vals"


#' Lambda^k values used to generate simulated data
#'
#'  a list of length 2, containing arrays of sizes c(30,2,3) and c(30,2,4), respectively, containing true values used in data generation.
#'  \code{Lambda[[1]][,,i]} is equivalent to \code{diag(B_vals[[1]][,i])\%\*\%Q_vals}
#'
#' @format A list
#' \describe{
#'     \item{Lambda_vals[[1]]}{an array of dimensions c(30,2,3).}
#'      \item{Lambda_vals[[2]]}{an array of dimensions c(30,2,4).}
#' }
"Lambda_vals"


#' a palette for visualizating reconstructed covariance matrices
#'
#' a vector of length 15, with
#' @format a vector of length 15
"bluered_palette"
