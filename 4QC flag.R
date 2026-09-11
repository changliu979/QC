# This code is written for readability rather than efficiency.
#
# Radar-assisted rainfall QC
#
# Required station_day columns:
#   stno, rr, east, north, elev, qc_flag
#
# Required radar_day columns:
#   east, north, radar_value, noise
#
# noise:
#   "Y" = radar noise, "N" = radar value

rm(list = ls())
gc()

library(RANN)
library(gstat)
library(sp)
library(dplyr)

# calculate radar9.
n_radar <- 9
radar_distance_limit <- (1500 * 1.05) * sqrt(2)
powers <- c(1, 1.5, 2, 2.5, 3)
model_list <- list(
  "east+north+elev" =rr ~ east + north + elev,
  "east+north+elev+radar9" =rr ~ east + north + elev + radar9)

# Station data after station-based QC.
station_day <- station_day[!station_day$qc_flag %in% c("qc_fail", "sampling"),]

# Radar-grid result after radar-noise detection.
radar_day <- radar_day[radar_day$noise == "N", ]

# New var with 3x3 mean rainfall
add_radar9 <- function(station_data,radar_grid,n_neighbours = 9,distance_limit) {

  nearest_radar <- nn2(
    data = as.matrix(radar_grid[, c("east", "north")]),
    query = as.matrix(station_data[, c("east", "north")]),k = 10)
  
  station_data$radar9 <- NA_real_
  
  for (i in seq_len(nrow(station_data))) {
    
    neighbour_idx <- nearest_radar$nn.idx[i, ]
    neighbour_dist <- nearest_radar$nn.dists[i, ]
    
    keep <- neighbour_idx > 0  &neighbour_dist <= distance_limit
    
    if (!any(keep)) {
      next
    }
    
    radar_values <- radar_grid$radar_value[ neighbour_idx[keep]]
    
    radar_values <- radar_values[is.finite(radar_values)]
    
    if (length(radar_values) == 0) {
      next
    }
    
    station_data$radar9[i] <- mean(radar_values)
  }
  
  station_data
}


station_radar <- add_radar9(station_data = station_day,radar_grid = radar_day,n_neighbours = n_radar,distance_limit = radar_distance_limit)
station_radar <- station_radar[complete.cases(station_radar),]


# LOOCV with NRMSE
LOO_function <- function(formula,idw_power,station_data) {
  
  predictions <- rep(NA_real_,nrow(station_data))
  
  for (s in seq_len(nrow(station_data))) {
    
    train <- station_data[-s, ]
    test <- station_data[s, , drop = FALSE]
    
    regression_model <- lm(formula,data = train)
    
    trend_prediction <- predict(regression_model,newdata = test)
    
    train$residual <- residuals(regression_model)
    
    train_sp <- train
    coordinates(train_sp) <- ~ east + north
    
    idw_model <- gstat(id = "residual",formula = residual ~ 1,data = train_sp,set = list(idp = idw_power),nmax = nrow(train_sp))
    
    test_sp <- test
    coordinates(test_sp) <- ~ east + north
    
    residual_prediction <- predict(idw_model, newdata = test_sp)$residual.pred
    
    combined_prediction <-trend_prediction +residual_prediction
    
    if (is.finite(combined_prediction)) {
      predictions[s] <- max(combined_prediction,0)
    }
  }
  
  predictions
}


model_results <- list()
result_number <- 1

mean_rainfall <- mean(station_radar$rr,na.rm = TRUE)

for (model_name in names(model_list)) {
  
  for (idw_power in powers) {
    
    loo_prediction <- LOO_function(formula = model_list[[model_name]],idw_power = idw_power,station_data = station_radar)
    
    rmse <- sqrt(mean((loo_prediction - station_radar$rr)^2,na.rm = TRUE))
    
    nrmse <- if (is.finite(mean_rainfall) &&mean_rainfall > 0) { 
      rmse / mean_rainfall
    } else {
      NA_real_
    }
    
    model_results[[result_number]] <- data.frame(
      model = model_name,
      idw_power = idw_power,
      rmse = rmse,
      nrmse = nrmse
    )
    
    result_number <- result_number + 1
  }
}

model_results <- do.call(rbind,model_results)

best_index <- which.min(model_results$rmse)

best_result <- model_results[best_index,]

best_model <- as.character(best_result$model)

best_power <- best_result$idw_power


model_results
best_result



best_prediction <- LOO_function(formula = model_list[[best_model]],idw_power = best_power,station_data = station_radar)
station_radar$pred <- round(best_prediction,1)
output <- merge(station_day,station_radar[, c("stno",'radar9', "pred")], by = "stno",all.x = TRUE)
output$model <- as.character(best_result$model)


output <- output %>%
  mutate(
    diff = ifelse(model == "east+north+elev+radar9" &!is.na(pred), abs(rr - pred), NA_real_),
    log_diff = ifelse(model == "east+north+elev+radar9" & !is.na(pred),abs(log(rr+1) - log(pred+1)), NA_real_),
    
    final_flag = case_when( qc_flag %in% c("sampling", "qc_fail") ~ "fail",
      model == "east+north+elev" & qc_flag == "qc_ok" ~ "ok",
      model == "east+north+elev" & qc_flag != "qc_ok" ~ "fail",
      model == "east+north+elev+radar9" & qc_flag != "qc_ok" & diff > 5.1 & log_diff > 0.5 ~ "fail",
      model == "east+north+elev+radar9" &  !is.na(diff) & !is.na(log_diff) ~ "ok",
      TRUE ~ "fail")
  )

output