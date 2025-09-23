get_coordinate_of_max_value <- function(r) {
  cell <- which.max(terra::values(r))
  return(terra::xyFromCell(r, cell))
}

get_mode_distance <- function(r, x) {
  # Calculate performance metric 2 from Truszkowski et al.
  # Distance from the position of maximum probability to the
  # true sample location.
  x_geo <- sf::st_as_sf(vect(x, crs=crs(r)))
  x_mode <- get_coordinate_of_max_value(r)
  x_mode_geo <- sf::st_as_sf(vect(x_mode, crs=crs(r)))
  return(as.numeric(sf::st_distance(x_geo, x_mode_geo)[1, 1]/1000))
}

get_mean_absolute_error <- function(r, x) {
  # Calculate performance metric 3 from Truszkowski et al.
  # x is a 1x2 matrix coordinate, r is a *probability 
  # density* raster.

  # Calculate distance raster
  pt <- terra::vect(x, crs=crs(r))
  d <- terra::distance(r, pt, unit="km")
  
  return(sum(terra::values(r * d)))
}

get_area_scored_higher <- function(r, x) {
  # Calculate performance metric 4 from Truszkowski et al.
  # r is a p-value or probability density raster.
  # x is a 1x2 matrix coordinate indicating the true harvest location.
  # This metric calculates the raster area where r > r[x].
  pt <- terra::vect(x, crs=crs(r))
  v  <- terra::extract(r, pt)
  
  r_gt <- (r > v[1, 2]) * terra::cellSize(r, unit="km")
  return(sum(terra::values(r_gt)))
}
