get_mean_pairwise_distance <- function(d) {
  d_copy <- as.matrix(d) # deep copy
  diag(d_copy) <- NA
  return(mean(d_copy, na.rm=TRUE))
}

get_upper_whisker <- function(x) {
  q3 <- quantile(x, 0.75)
  q1 <- quantile(x, 0.25)
  iqr <- q3 - q1
  
  return(min(max(x), q3 + 1.5*iqr))
}

get_out_of_fold_distance <- function(x, folds) {
  # Calculate distance between data points
  # in x after Meyer and Pebesma.
  # x is assumed to be a nx2 coordinate matrix
  
  # Get distance matrix
  x_geo <- sf::st_as_sf(data.frame(x), coords=c(1, 2))
  x_geo <- sf::st_set_crs(x_geo, 4326)
  d <- sf::st_distance(x_geo)
  
  # Overall mean pairwise distance
  meandist <- get_mean_pairwise_distance(d)
  
  # Calculate what cells are in the same fold.
  same_fold <- outer(folds, folds, FUN=function(x, y) x == y)
  
  # Set distance to infinity if in the same fold
  d[same_fold] <- Inf
  
  # Calculate the minimum distance to another point that is not
  # in the same fold.
  mindist <- apply(d, 1, min)
  
  return(mindist)
}
