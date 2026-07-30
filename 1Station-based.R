# This code is written for readability rather than efficiency.

# Station-based QC
#
# dff contain:
#   stno (unique station number), rr (station rainfall), east, north, elev, type (station type), diff_rr_nn (abs rainfall difference between near neighbour station), flag
#
# type 1 and type 2 stations are used as the reference stations (synoptic stations).

n_boot <- 30000
reference_types <- c("1", "2")
residual_limit <- 3
rainfall_difference_limit <- 5
reference_difference_limit <- 50
qc_fail_limit <- 250

df <- dff
df$type <- as.character(df$type)

# initial check
df$qc_flag <- "qc_ok"
valid_sampling_flags <- c(0, 1, 3, 4, 5, 7, 10, 12) # flag from raw data
df$qc_flag[is.na(df$rr) |is.na(df$flag) |!df$flag %in% valid_sampling_flags] <- "sampling"
df$qc_flag[!is.na(df$rr) & df$rr > qc_fail_limit] <- "qc_fail"

dff <- df[df$qc_flag == "qc_ok" & complete.cases(df[, c( "rr", "east", "north","elev", "type", "diff_rr_nn")]),]

# 1. Regression check
fit_all <- lm(rr ~ east + north + elev, data = dff)
dff$residual <- residuals(fit_all)

residual_sd <- sd(dff$residual)
dff$regression_flag <- abs(dff$residual) > residual_limit * residual_sd

dff$reference_reject <- dff$type %in% reference_types & dff$regression_flag & dff$diff_rr_nn > reference_difference_limit

reference_data <- dff[dff$type %in% reference_types & !dff$reference_reject,]

check_data <- dff[!dff$type %in% reference_types, ]

# 2. Bootstrap check
bootstrap_coef <- matrix(NA_real_, nrow = n_boot,ncol = 4,dimnames = list(NULL, c("intercept", "east", "north", "elev")))

for (i in seq_len(n_boot)) {
  
  bootstrap_sample <- do.call(
    rbind,
    lapply(reference_types, function(station_type) {
      
      group <- reference_data[reference_data$type == station_type, ]
      group[sample(seq_len(nrow(group)), replace = TRUE), ]
    })
  )
  
  fit_boot <- lm(rr ~ east + north + elev, data = bootstrap_sample)
  bootstrap_coef[i, ] <- coef(fit_boot)
}

check_data$bootstrap_flag <- FALSE

for (i in seq_len(nrow(check_data))) {
  
  predicted_rr <-bootstrap_coef[, "intercept"] + check_data$east[i] * bootstrap_coef[, "east"] + check_data$north[i] * bootstrap_coef[, "north"] + check_data$elev[i] * bootstrap_coef[, "elev"]
  
  check_data$bootstrap_flag[i] <- check_data$rr[i] < min(predicted_rr, na.rm = TRUE) | check_data$rr[i] > max(predicted_rr, na.rm = TRUE)
}

dff$bootstrap_flag <- FALSE
dff$bootstrap_flag[match(check_data$stno, dff$stno)] <-check_data$bootstrap_flag

# 3. Flag
dff$bootstrap_flag <- dff$bootstrap_flag & dff$diff_rr_nn >= rainfall_difference_limit

dff$qc_flag <- "qc_ok"
dff$qc_flag[dff$regression_flag & !dff$bootstrap_flag] <- "regression"
dff$qc_flag[!dff$regression_flag & dff$bootstrap_flag] <- "bootstrap"
dff$qc_flag[dff$regression_flag & dff$bootstrap_flag] <- "regression & bootstrap"
dff$qc_flag[dff$reference_reject] <- "regression_reject"


model_results <- dff[,c("stno","qc_flag")]

model_rows <- match(dff$stno, df$stno)
df$qc_flag[model_rows] <- dff$qc_flag

output <- df[ ,c("stno", "rr", "east", "north", "elev", "type", "diff_rr_nn", "qc_flag")]

output
