#########################

# SPATIAL COVARIANCE PLOT 

#########################

# Load necessary library
library(dplyr)

# Set seed for reproducibility
set.seed(123)

# Create simulated data
df <- data.frame(
  lon = runif(700, min = -35, max = 35), 
  lat = runif(700, min = -40, max = 40), 
  SCCS_ID = sample(1:700, 700, replace = TRUE)
)

# load scripts
source("scripts/HS_varcov.spatial.3D.R")
source("scripts/space_correl.R")

# Function to calculate distances

space<-function(d){ # d has lon, lat, and SCCS_ID as columns
  # Calculate distances between points using the haversine_distance function
  distances <- matrix(NA, nrow = nrow(d), ncol = nrow(d))
  for (i in 1:nrow(d)) {
    for (j in 1:nrow(d)) {
      distances[i, j] <- haversine_distance(d$lon[i], d$lat[i], d$lon[j], d$lat[j])
    }
  }
  
  # Normalize distances to get correlations and set diagonal to 1
  correlation_matrix <-  distances
  
  rownames(correlation_matrix) <- colnames(correlation_matrix) <- d$SCCS_ID
  # correlation_matrix will have perfect correlations (1) on the diagonal, and the off-diagonal elements represent the correlations between different locations based on their great-circle distances.
  return(correlation_matrix)
}

# Calculate distances
spatial_mat <- space(df)
Distances <- as.list(spatial_mat)

# Calculate covariances using various parameterisations of a Matérn covariance function 
cov1 <- varcov.spatial.3D(df[,c("lon", "lat")], cov.pars = c(1, 1.15), kappa = 2)$varcov
cov2 <- varcov.spatial.3D(df[,c("lon", "lat")], cov.pars = c(1, 2), kappa = 2)$varcov
cov3 <- varcov.spatial.3D(df[,c("lon", "lat")], cov.pars = c(1, 3), kappa = 2.5)$varcov
cov4 <- varcov.spatial.3D(df[,c("lon", "lat")], cov.pars = c(1, 5), kappa = 5)$varcov

# Calculate inverse spatial function
n <- 1:20000
reverse <- 1/n

# PLOT
plot(n, reverse,
     xlab = "Haversine distance (km)", ylab = "Covariance", 
     col="#999999FF", type = "l", lwd = 2,
     xlim = c(0, 10000))
points(Distances, cov1, pch = 19, cex = 0.1, col = "#E41A1CFF")
points(Distances, cov2, pch = 19, cex = 0.1, col = "#377EB8FF")
points(Distances, cov3, pch = 19, cex = 0.1, col = "#4DAF4AFF")
points(Distances, cov4, pch = 19, cex = 0.1, col = "#984EA3FF")

# Note: legend added manually in Inkscape
