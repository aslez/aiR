% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/south_df.R
\docType{data}
\name{south_df}
\alias{south_df}
\title{Overlap between 1860 and 2000 counties.}
\format{A data frame with 3826 rows and 3 variables:
\itemize{
  \item ID1860. NHGIS identifier (GISJOIN) for 1860.
  \item ID2000. NHGIS identifier (GISJOIN) for 2000.
  \item AREA.  Area of overlap measured in square miles.
}}
\source{
\url{https://www.nhgis.org/}
}
\usage{
data(south_df)
}
\description{
A dataset containing data on group overlap based on the intersection of
United States counties from 1860 and 2000.  The original geographic data
were drawn from the National Historical Geographic Information System
(NHGIS).  Intersections and areas were calculated using QGIS.
}
