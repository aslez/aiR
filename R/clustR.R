#'Identify common geographies
#'
#'\code{clustR} identifies common geographies based on patterns of group overlap.
#'Input can take one of two forms: partitions or intersections.  A set of
#'partitions is represented by a list of two or more 
#'\code{SpatialPolygonsDataFrame} objects, each of which is composed of a set 
#'of areal units (e.g., counties, census tracts).  To identify common 
#'geographies \code{clustR} calculates the intersection of these partitions.  
#'This can be very slow depending on the number of observations and the 
#'resolution of the underlying boundary files.  Alternatively, intersections 
#'can be calculated using a dedicated GIS (e.g., ArcGIS, QGIS) and then passesd 
#'to \code{clustR} (recommended).  Intersections are represented as a single 
#'object.  This can be either a \code{SpatialPolygonsDataFrame} object or a
#'\code{data.frame} object.  \code{mp_shp} and \code{mp_int} are internal 
#'helper functions used to construct properly formatted membership profiles.

#'@param x Either a list of two or more \code{SpatialPolygonsDataFrame} objects 
#'or a single object depicting the intersection between partitions.  
#'Intersections can represented using either an \code{SpatialPolygonsDataFrame} 
#'or a \code{data.frame}.  When intersections are represented as a 
#'\code{SpatialPolygonsDataFrame}, the area of each intersection is calculated 
#'on the fly.  When intersections are represented as a \code{data.frame}, the 
#'area is included in the \code{data.frame} in question.  WARNING!!  PROCEED 
#'WITH EXTREME CAUTION WHEN WORKING WITH ORIGINAL PARTITIONS.  CHANGES IN THE 
#'RGEOS INTERSECTION ROUTINE ARE CAUSING POLYGONS TO BE DROPPED, LEADING TO 
#'INCORRECT CLUSTERS. 
#'
#'@param nid A character vector containing the column names used to identify
#'groups within each partition.  This is only required when starting with 
#'intersections as opposed to paritions.  When starting with partitions, 
#'\code{clustR} will default to polygon IDs unless otherwise indicated through
#'the \code{nid} argument.  Groups should be uniquely identified 
#'within partitions.
#'
#'@param area A string containing the name of the column containing data on the
#'area of overlap between groups.  This is only required when using a single
#'\code{data.frame} object.
#'
#'@param thresh A number between 0 and 1 used to drop ties resulting from
#'spurious polygons.  This value represents the area of group overlap expressed
#'as a proportion of the area of the smallest overlapping unit in question.
#'
#'@details \code{clustR} assigns areal units to common geographies using the
#'method outlined by Slez, O'Connell, and Curtis (2014), who show that
#'identifying the common geographies associated with a set of \eqn{k} partitions
#'is identical to identifying the components of a \eqn{k}-uniform
#'\eqn{k}-partite hypergraph.  Each edge in the hypergraph represents a
#'membership profile depicting the intersection between areal units.
#'
#'@return A \code{data.frame} depicting the relationship between hyperedges and
#'components.  Each hyperedge consists of a membership profile containing the
#'name of one group from each partition.  Each component refers to a common
#'geography.
#'
#'@references Slez, A., O'Connell, H.A., and Curtis, K.J.  2014.  "A Note on the
#'Identification of Common Geographies."
#'
#'@examples
#'#load list containing partitions
#'data(nd_list)
#'clustR(nd_list)
#'
#'#load data frame containing intersections
#'data(example)
#'
#'#add placeholder for area (real areas not needed)
#'example$AREA <- 1
#'clustR(example, nid = c("A", "B"), area = "AREA")
#'
#'#load data frame containing intersections
#'data(south_df)
#'clustR(south_df, nid = c("ID1860", "ID2000"), area = "AREA")


clustR <- function(x, nid = NULL, area = NULL, thresh = .05) {
  #Step 1: Construct Intersections
  if(class(x) == "list") mp <- mp_shp(x, nid)
  else mp <- mp_int(x, nid, area)

  #Step 2: Construct Affiliation Matrix
  df <- mp %>%
    tidyr::gather(part, group, -area, -edge) %>%
    dplyr::group_by(part, group) %>%
    dplyr::mutate(std_area = area / sum(area)) %>%
    dplyr::group_by(edge) %>%
    dplyr::mutate(tie = (max(std_area) >= thresh) * 1) %>%
    dplyr::filter(tie == 1) %>%
    dplyr::mutate(vertex = paste(part, group, sep = "_"))
  A_df <- data.frame(model.matrix(~ factor(edge) - 1, data = df),
                     vertex = df$vertex)
  A_df <- A_df %>%
    dplyr::group_by(vertex) %>%
    dplyr::summarise_all(dplyr::funs(sum)) %>%
    dplyr::select(-vertex)

  #Step 3: Identify Components
  E <- t(as.matrix(A_df)) %*% as.matrix(A_df)
  dimnames(E) <- list(unique(df$edge), unique(df$edge))
  g <- igraph::graph.adjacency(E)
  clst_df <- data.frame(edge = as.numeric(igraph::V(g)$name),
                        comp = igraph::clusters(g)$members)
  out_df <- df %>%
    dplyr::ungroup() %>%
    dplyr::left_join(clst_df, by = "edge") %>%
    dplyr::select(edge, part, group, comp) %>%
    tidyr::spread(part, group)
  out_df
}

#"@rdname clustR
mp_shp <- function(x, nid) {
  if (!all(sapply(x, class) == "SpatialPolygonsDataFrame")) stop("If x is list, each object in x should be a SpatialPolygonsDataFrame.")
  np <- length(x)
  axb <- rgeos::gIntersection(x[[1]], x[[2]],
                       byid = TRUE, drop_lower_td = TRUE)
  if (np > 2) {
    for (i in 3:length(x)) {
      axb <- rgeos::gIntersection(axb, x[[i]],
                           byid = TRUE, drop_lower_td = TRUE)
    }
  }
  area <- rgeos::gArea(axb, byid = TRUE)
  mp <- data.frame(mp = names(area), edge = 1:length(area), area)
  df <- tidyr::separate(mp, col = mp, into = paste0("p", 1:np))
  if (!is.null(nid)) {
    for (i in 1:np) {
      df[, paste0("p", i)] <- x[[i]]@data[df[, paste0("p", i)], nid]
    }
  }
  df
}

#"@rdname clustR
mp_int <- function(x, nid, area) {
  if (!class(x) %in% c("SpatialPolygonsDataFrame", "data.frame"))
    stop("Class of object in x is not valid.")
  if (class(x) == "SpatialPolygonsDataFrame") {
    if (!is.null(area)) warning("Polygon areas are being calculated directly.  Area field is ignored.")
    new_area <- rgeos::gArea(x, byid = TRUE)
  }
  if (class(x) == "data.frame") {
    if (is.null(area)) stop("Missing area field.")
    if (length(grep(area, names(x))) == 0) stop("Area field not found.")
    if (length(grep(area, names(x))) != 0) new_area <- x[, area]
  }
  mp <- data.frame(x[, nid], edge = 1:NROW(x), new_area)
  names(mp) <- c(nid, "edge", "area")
  mp
}
