
# Function to calculate great-circle distances between two points on the Earth's surface
haversine_distance <- function(lon1, lat1, lon2, lat2) {
  # Convert latitude and longitude from degrees to radians
  lon1 <- lon1 * pi / 180
  lat1 <- lat1 * pi / 180
  lon2 <- lon2 * pi / 180
  lat2 <- lat2 * pi / 180
  
  # Haversine formula to calculate great-circle distance
  dlon <- lon2 - lon1 
  dlat <- lat2 - lat1 
  a <- sin(dlat/2)^2 + cos(lat1) * cos(lat2) * sin(dlon/2)^2 
  c <- 2 * asin(sqrt(a)) 
  r <- 6371  # Radius of the Earth in kilometers
  return(c * r)
}

space_mat<-function(d){ # d has lon, lat, and SCCS_ID as columns
  # Calculate distances between points using the haversine_distance function
  distances <- matrix(NA, nrow = nrow(d), ncol = nrow(d))
  for (i in 1:nrow(d)) {
    for (j in 1:nrow(d)) {
      distances[i, j] <- haversine_distance(d$lon[i], d$lat[i], d$lon[j], d$lat[j])
    }
  }
  
  # Normalize distances to get correlations and set diagonal to 1
  correlation_matrix <- 1 / distances
  diag(correlation_matrix) <- 1
  
  rownames(correlation_matrix) <- colnames(correlation_matrix) <- d$SCCS_ID
  # correlation_matrix will have perfect correlations (1) on the diagonal, and the off-diagonal elements represent the correlations between different locations based on their great-circle distances.
  return(correlation_matrix)
}


space_mat_varc<-function(d){ # d has lon, lat, and SCCS_ID as columns

  # Calculate distances between points using the haversine_distance function
  distances <- matrix(NA, nrow = nrow(d), ncol = nrow(d))
  for (i in 1:nrow(d)) {
    for (j in 1:nrow(d)) {
      distances[i, j] <- haversine_distance(d$lon[i], d$lat[i], d$lon[j], d$lat[j])
    }
  }
  
  kappa = 2 # smoothness parameter as recommended by Dinnage et al. (2020)
  phi = c(1, 1.15) # Sigma parameter. First value is not used.
  
  lower_triangle <- distances[lower.tri(distances, diag = FALSE)]
  smatrix3<-varcov.spatial(dists.lowertri=lower_triangle, cov.pars = phi, kappa = kappa)$varcov # this function can also just take a distance matrix
  rownames(smatrix3) <- colnames(smatrix3) <- d$SCCS_ID
  smatrix3<-smatrix3/max(smatrix3)
  return(smatrix3)
}