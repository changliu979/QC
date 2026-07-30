# This code is written for readability rather than efficiency.
#
# Radar-assisted rainfall QC: Radar noise detect
#
# radar_day includes:
#   east, north, radar_value (24 hour radar accumulations), rain_grid (daily rainfall grids)
#

rm(list = ls())
gc()

library(RANN)
library(sf)

radius <- 2000
tolerance <- 100
n_permutations <- 100000

# Identify radar indicated rainfall areas.
close_function <- function(df, threshold, radius = 2000, tol = 100) {
  coords <- as.matrix(df[, c("east", "north")])
  
  nn_result <- nn2(data = coords,query = coords, searchtype = "radius",radius = radius + tol,k = 25)
  
  df$cluster_count_2km <- 0
  
  for (i in seq_len(nrow(df))) {
    neighbour_idx <- nn_result$nn.idx[i, ]
    neighbour_dist <- nn_result$nn.dists[i, ]
    
    keep <- which(neighbour_idx > 0 &!is.na(neighbour_dist) & neighbour_dist < (radius + tol))
    
    if (length(keep) == 0) next
    
    idx <- neighbour_idx[keep]
    local_mean <- mean(df$radar_value[idx], na.rm = TRUE)
    
    if (local_mean >= threshold) {
      df$cluster_count_2km[idx] <- df$cluster_count_2km[idx] + 1
    }
  }
  
  df
}

# merge clusters
onecluster_function <- function(df, dist_target = 2000, tol = 100) {
  df$id <- seq_len(nrow(df))
  df$cluster_num_2km <- 0
  
  cluster_df <- df[df$cluster_count_2km > 0, ]
  if (nrow(cluster_df) == 0) return(df)
  
  cluster_df$cluster_num_2km <- 0
  cluster_number <- 1
  
  while (any(cluster_df$cluster_num_2km == 0)) {
    seed_id <- cluster_df$id[cluster_df$cluster_num_2km == 0][1]
    cluster_ids <- seed_id
    growing <- TRUE
    
    while (growing) {
      current <- cluster_df[cluster_df$id %in% cluster_ids, ]
      remaining <- cluster_df[!cluster_df$id %in% cluster_ids, ]
      
      if (nrow(remaining) == 0) break
      
      nn <- nn2(data = remaining[, c("east", "north")],query = current[, c("east", "north")],k = min(nrow(remaining), nrow(current)))
      
      close_idx <- as.vector(nn$nn.idx)[as.vector(nn$nn.dists) < (dist_target + tol)]
      close_idx <- close_idx[close_idx > 0]
      
      if (length(close_idx) == 0) break
      
      new_ids <- unique(remaining$id[close_idx])
      previous_size <- length(cluster_ids)
      cluster_ids <- unique(c(cluster_ids, new_ids))
      growing <- length(cluster_ids) > previous_size
    }
    
    cluster_df$cluster_num_2km[cluster_df$id %in% cluster_ids] <- cluster_number
    
    cluster_number <- cluster_number + 1
  }
  
  df$cluster_num_2km[match(cluster_df$id, df$id)] <- cluster_df$cluster_num_2km
  
  df
}

# permuation test
run_permutation_test <- function(df,dist_target = 2000,tol = 100,n_perms = 100000) {
  df$id <- seq_len(nrow(df))
  
  cluster_ids <- sort(unique(df$cluster_num_2km[df$cluster_num_2km > 0]))
  if (length(cluster_ids) == 0) return(data.frame())
  
  spatial_df <- st_as_sf(df,coords = c("east", "north"))
  
  all_clustered_ids <- df$id[df$cluster_num_2km > 0]
  
  results <- lapply(cluster_ids, function(k) {
    cluster_points <- df[df$cluster_num_2km == k, ]
    cluster_sf <- spatial_df[spatial_df$cluster_num_2km == k, ]
    
    distance_matrix <- st_distance(cluster_sf, spatial_df)
    
    neighbour_position <- which(distance_matrix >= (dist_target - tol) &distance_matrix <= (dist_target + tol),arr.ind = TRUE)
    
    neighbour_ids <- if (nrow(neighbour_position) == 0) {
      integer(0)
    } else {
      unique(spatial_df$id[neighbour_position[, 2]])
    }
    
    neighbour_ids <- setdiff(neighbour_ids, all_clustered_ids)
    neighbour_points <- df[df$id %in% neighbour_ids, ]
    
    inside_values <- cluster_points$rain_grid
    outside_values <- neighbour_points$rain_grid
    
    if (length(outside_values) == 0 || length(inside_values) == 0) {
      p_value <- NA_real_
    } else {
      observed_difference <-mean(outside_values) - mean(inside_values)
      
      combined_values <- c(outside_values, inside_values)
      n_outside <- length(outside_values)
      
      permuted_difference <- numeric(n_perms)
      
      for (s in seq_len(n_perms)) {
        shuffled <- sample(combined_values)
        permuted_difference[s] <-mean(shuffled[seq_len(n_outside)]) -mean(shuffled[(n_outside + 1):length(shuffled)])
      }
      p_value <-(sum(abs(permuted_difference) >= abs(observed_difference)) + 1) /(length(permuted_difference) + 1)
    }
    
    inside_rain <- mean(cluster_points$rain_grid)
    outside_rain <- if (nrow(neighbour_points) > 0) {
      mean(neighbour_points$rain_grid)
    } else {
      NA_real_
    }
    
    inside_radar <- mean(cluster_points$radar_value)
    outside_radar <- if (nrow(neighbour_points) > 0) {
      mean(neighbour_points$radar_value)
    } else {
      NA_real_
    }
    
    data.frame(
      cluster = k,
      size = nrow(cluster_points),
      mean_east = mean(cluster_points$east),
      mean_north = mean(cluster_points$north),
      p_value = p_value,
      rain_inside = inside_rain,
      rain_outside = outside_rain,
      radar_inside = inside_radar,
      radar_outside = outside_radar,
      rain_ratio = inside_rain / outside_rain,
      radar_ratio = inside_radar / outside_radar
    )
  })
  
  do.call(rbind, results)
}


# Upper boxplot whisker used as the threshold
radar_threshold <- boxplot.stats(radar_day$radar_value)$stats[5]

# Identify radar indicated high rainfall area and clusters
radar_day <- close_function(
  radar_day,
  threshold = radar_threshold,
  radius = radius,
  tol = tolerance
)

radar_day <- onecluster_function(
  radar_day,
  dist_target = radius,
  tol = tolerance
)

# permutation test
output <- run_permutation_test(
  radar_day,
  dist_target = radius,
  tol = tolerance,
  n_perms = n_permutations
)


# noise detect
output$noise <- "Y"
n_tests <- nrow(output) # Bonferroni correction 
not_noise <- !is.na(output$p_value) & output$p_value < (0.05 / n_tests) & output$rain_ratio > 1 & output$radar_ratio < 10
output$noise[not_noise] <- "N"

# merge with radar_day
result <- merge(radar_day,output[, c("cluster", "noise")],
  by.x = "cluster_num_2km",
  by.y = "cluster",
  all.x = TRUE
)

result$noise[is.na(result$noise)] <- "N"
result <- result[,c("east","north","radar_value","rain_grid","cluster_num_2km","noise")]

colnames(result) <- c("east","north","radar_value","rain_grid","cluster","noise")

result
