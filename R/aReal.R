#'Areal interpolation
#'
#'\code{aReal} uses areal interpolation to project a set of source values on to 
#'a set of target observations.

#'@param src A \code{SpatialPolygonsDataFrame} from which to project.
#'
#'@param trg A target geography on which to project values.
#'
#'@param vals A vector containing the names of the columns in the \code{data}
#'slot of the \code{src} object to be used as source values.
#'
#'@details Using the method proposed by Markoff and Shapiro (1973), \code{aReal} 
#'projects source values on to a target geography.  Source values are assumed
#'to be counts.
#'
#'@return A \code{data.frame} containing the projected values along with
#'polygon IDs drawn from \code{trg}, along with information on the structure
#'of the common geographies within which values were interpolated.
#'
#'@references Markoff, J. and Shapiro, G.  1973.  "The Linkage of Data
#'Describing Overlapping Geographical Units."  \emph{Historical Methods 
#'Newsletter} 7: 34--46.
#'
#'@examples
#'#load list containing partitions
#'#data(nd_list)
#'#aReal(nd_list[[1]], nd_list[[2]], nd_list[[1]]$wheat)

#'@rdname aReal
aReal <- function(src, trg, vals) {
  #generate W matrix
  src_names <- sapply(slot(src, 'polygons'), function(x) slot(x, 'ID'))
  trg_names <- sapply(slot(trg, 'polygons'), function(x) slot(x, 'ID'))
  int_shp <- gIntersection(src, trg, byid = TRUE, drop_not_poly = TRUE)
  areas <- gArea(int_shp, byid = TRUE)
  ids <- do.call(rbind, strsplit(names(areas), ' '))
  W <- matrix(0, nrow = length(src_names), ncol = length(trg_names))
  dimnames(W) <- list(src_names, trg_names)
  for (i in 1:nrow(ids)) {
    W[ids[i, 1], ids[i, 2]] <- areas[i] 
  }
  W <- W / rowSums(W)
  
  #identify clusters using W
  g_mat <- t(W) %*% W
  diag(g_mat) <- 0
  g_mat[g_mat > 0] <- 1
  clst <- igraph::clusters(igraph::graph.adjacency(g_mat))
  cluster <- clst$membership
  csize <- clst$csize[cluster]
  
  #build data.frame
  df <- do.call(data.frame, lapply(vals, function(x) tWy(W, src@data[, x])))
  names(df) <- vals
  df <- df %>%
    dplyr::mutate(id = rownames(df),
                  cluster = cluster,
                  csize = csize)
  df
}

#'@rdname aReal
tWy <- function(W, y) {
  y[is.na(y)] <- 0
  t(W) %*% y
}

