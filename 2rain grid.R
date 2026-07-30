# This code is written for readability rather than efficiency.
#
# Radar-assisted rainfall QC: provisional rainfall grids
#
# station_day is from station-based qc
#
# radar-day is from raw radar dataset including: east, north, elev and radar_value for each grid

rm(list = ls())
gc()

library(gstat)
library(sp)

# Candidate IDW powers.
powers <- c(1, 1.5, 2, 2.5, 3)

# Spatial regression model
rainfall_model <- rr ~ east + north + elev

# Leave-one-out prediction for one IDW power.
LOO_function <- function(station_data,formula,idw_power) {
  
  predictions <- numeric(nrow(station_data))
  
  for (s in seq_len(nrow(station_data))) {
    train <- station_data[-s, ]
    test <- station_data[s, , drop = FALSE]

    regression_model <- lm(formula,data = train)
    
    trend_prediction <- predict(regression_model,newdata = test)
    
    train$residual <- residuals(regression_model)
    
    train_sp <- train
    coordinates(train_sp) <- ~ east + north
    
    test_sp <- test
    coordinates(test_sp) <- ~ east + north
    
    idw_model <- gstat(id = "residual",formula = residual ~ 1,data = train_sp,set = list(idp = idw_power),nmax = nrow(train_sp))
    
    residual_prediction <- predict(idw_model,newdata = test_sp)$residual.pred
    
    predictions[s] <- max((trend_prediction +residual_prediction),0)
  }
  
  predictions
}


# Generate the rainfall grid using the selected IDW power.
predict_rainfall_grid <- function(station_data,radar_grid,formula,idw_power) {
  
  regression_model <- lm(formula,data = station_data)
  trend_prediction <- predict(regression_model,newdata = radar_grid)
  
  station_data$residual <- residuals(regression_model)
  station_sp <- station_data
  coordinates(station_sp) <- ~ east + north
  
  idw_model <- gstat(id = "residual",formula = residual ~ 1,data = station_sp,set = list(idp = idw_power),nmax = nrow(station_sp))
  
  radar_sp <- radar_grid
  coordinates(radar_sp) <- ~ east + north
  
  residual_prediction <- predict(idw_model,newdata = radar_sp)$residual.pred
  
  radar_grid$rain_grid <- pmax(as.numeric(trend_prediction) +as.numeric(residual_prediction),0)
  
  radar_grid
}


# Evaluate each candidate IDW power using leave-one-out.
rmse <- numeric(length(powers))

for (j in seq_along(powers)) {
  
  loo_prediction <- LOO_function(station_data = station_day,formula = rainfall_model,idw_power = powers[j])
  rmse[j] <- sqrt(mean((loo_prediction - station_day$rr)^2, na.rm = TRUE))
  
}


best_power <- powers[which.min(rmse)]

power_result <- data.frame(idw_power = powers,rmse = rmse)

radar_day <- predict_rainfall_grid(station_data = station_day,radar_grid = radar_day,formula = rainfall_model,idw_power = best_power)

output <- radar_day[,c("east","north","radar_value","rain_grid")]

output