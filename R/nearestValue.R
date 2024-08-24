#' @title nearestValue: Match species location data to environmental raster layers
#' @description Obtain the nearest covariate value at each of the species locations.
#' @param pts The species location data
#' @param r The raster file to extract the covariates at.
#' @param ... Extra arguments for \link[FNN]{knnx.index}.
#' @import FNN
#' @export

nearestValue <- function(pts, r,...) {
  # Convert SpatRaster to points
  rxy <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
  # Get the nearest raster cell index for each point
  ind <- FNN::knnx.index(rxy[, c("x", "y")], pts, k = 1, ...)[, 1]
  # Return the raster values corresponding to the nearest points
  rxy[ind, 3]
}