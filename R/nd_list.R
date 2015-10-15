#' North Dakota counties from 1890, 1900, and 1930.
#'
#' A list of three \code{SpatialPolygonsDataFrame} objects representing county
#' boundaries in North Dakota in 1890, 1900, and 1930 (polygons represent
#' boundaries as of October 31 of each year in question).  Data were drawn from
#' the Atlas of Historic County Boundaries Project (AHCBP).  Each
#' \code{SpatialPolygonsDataFrame} includes a single variable, ID_NUM, used as
#' a unique identifier by the AHCBP.
#'
#' @docType data
#' @name nd_list
#' @usage data(nd_list)
#' @format A list containing three \code{SpatialPolygonsDataFrame} objects.
#' The contents of the corresponding data slots are as follows:
#' \itemize{
#'   \item y1890.  58 rows and 1 column (ID_NUM).
#'   \item y1900.  40 rows and 1 column (ID_NUM).
#'   \item y1930.  53 rows and 1 column (ID_NUM).
#' }
#' @source \url{http://publications.newberry.org/ahcbp/}
NULL
