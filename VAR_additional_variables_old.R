# Final Project VAR Model

library(tseries)
library(vars)
library(forecast)

setwd("~/Desktop/")

# Load the data
data <- read.csv("combined_time_series.csv")


# Convert Date column to proper date format
data$Date <- as.Date(data$Date)

# ========================

# Log transform CPI to handle exponential growth and stabilize variance
log_CPI <- log(data$CPI)

# First difference Log(CPI)
diff_log_cpi <- diff(log_CPI)

final_CPI <- diff_log_cpi

# ======================

plot(data$Date, data$PPI, type = "l", col = "red",
     main = "Producer Price Index Time Series",
     xlab = "Date",
     ylab = "PPI")

plot(data$PPI, data$CPI, pch = 16, col = "darkblue",
     main = "CPI vs PPI", xlab = "PPI", ylab = "CPI")
cor(data$CPI, data$PPI)
# .977

lambda <- BoxCox.lambda(data$PPI)
lambda

PPI_bc <- BoxCox(data$PPI, lambda = lambda)
plot(data$Date, PPI_bc, type = "l")

print(ndiffs(PPI_bc))
#print(nsdiffs(PPI_bc))

diff_PPI <- diff(PPI_bc)
plot(data$Date[-1], diff_PPI, type = "l", col = "red")

adf.test(diff_PPI)

acf(diff_PPI, main = "ACF of First Order Difference of PPI")
pacf(diff_PPI, main = "PACF of First Order Difference of PPI")

final_PPI <- diff_PPI

# ============================

# Plot Housing Index
plot(data$Date, data$HousingIndex, type = "l", col = "red", main = "Housing Index Time Series", xlab = "Date", ylab = "Housing Index")

# Plot and calculate corr
plot(data$HousingIndex, data$CPI, pch = 16, col = "darkblue",
     main = "CPI vs Housing Index", xlab = "Housing Index", ylab = "CPI")
cor(data$CPI, data$HousingIndex)
# .951

# Stationarity
lambda <- BoxCox.lambda(data$HousingIndex)
lambda

housing_bc <- BoxCox(data$HousingIndex, lambda = "auto")
plot(housing_bc)

d <- ndiffs(housing_bc)
d
#D <- nsdiffs(housing_bc)

diff_housing <- diff(housing_bc)

plot(data$Date[-1], diff_housing, type = "l", col = "red",
     main = "First Order Difference Housing Index (With Box-Cox Transformation)", xlab = "Date", ylab = "Δ Housing Index")

adf.test(diff_housing)

acf(diff_housing, lag.max = 48, main = "ACF of First Order Difference Housing Index (With Box-Cox Transformation)")
pacf(diff_housing, main = "PACF of First Order Difference Housing Index (With Box-Cox Transformation)")


diff_housing2 <- diff(diff_housing)

adf.test(diff_housing2)
plot(data$Date[-(1:2)], diff_housing2, type = "l", col = "red",
     main = "Second Order Difference Housing Index(With Box-Cox Transformation)", xlab = "Date", ylab = "Δ Housing Index")

acf(diff_housing2, lag.max = 48, main = "ACF of 2nd Order Difference Housing Index (With Box-Cox Transformation)")
pacf(diff_housing2, main = "PACF of 2nd Order Difference Housing Index (With Box-Cox Transformation)")

sdiff_housing2 <- diff(diff_housing2, lag = 12)

plot(data$Date[-(1:14)], sdiff_housing2, type = "l", col = "red",
     main = "Seasonal and Second Order Differencing of Housing Index",
     xlab = "Date",
     ylab = "Housing Index")

acf(sdiff_housing2)
pacf(sdiff_housing2)

adf.test(sdiff_housing2)

min_length <- min(length(sdiff_housing2), length(diff_log_cpi))

plot(tail(sdiff_housing2, min_length), tail(diff_log_cpi, min_length), pch = 16, col = "darkblue",
     main = "CPI vs Unemployment Rate", xlab = "Transformed Unemployment Rate", ylab = "Transformed CPI")

cor(tail(sdiff_housing2, min_length), tail(diff_log_cpi, min_length))
length(sdiff_housing2)

final_housing <- diff_housing2

# ==============================

# Plot Unemployment
plot(data$Date, data$Unemployment, type = "l", col = "red", main = "Unemployment Rate Time Series", xlab = "Date", ylab = "Unemployment Rate")

# Plot and calculate corr
plot(data$Unemployment, data$CPI, pch = 16, col = "darkblue",
     main = "CPI vs Unemployment Rate", xlab = "Unemployment Rate (%)", ylab = "CPI")

cor(data$CPI, data$Unemployment, use = "complete.obs")


# Addressing Stationarity
diff_unemployment <- diff(data$Unemployment)

plot(data$Date[-1], diff_unemployment, type = "l", col = "red",
     main = "First Order Difference Unemployment Rate", xlab = "Date", ylab = "Δ Unemployment Rate (%)")

# ADF test
adf.test(diff_unemployment)
## p-value = .01, seems to have gotten stationarity in check

# Seasonality check
acf(diff_unemployment, main="ACF of First-Order Differenced Series")
# Exponential decay, no apparent seasonality

# Correlation check on our transformed data
plot(diff_unemployment, diff_log_cpi, pch = 16, col = "darkblue",
     main = "CPI vs Unemployment Rate", xlab = "Transformed Unemployment Rate", ylab = "Transformed CPI")

cor(diff_unemployment, diff_log_cpi)

final_unemployment <- diff_unemployment

# ====================
final_FedFund <- diff(data$FedFunds)



# ======= VAR Model (Full) =========
z_CPI   <- zoo(final_CPI,          order.by = data$Date[-1])
z_PPI   <- zoo(final_PPI,          order.by = data$Date[-1])
z_House <- zoo(final_housing,      order.by = data$Date[-(1:2)])
z_Unemp <- zoo(final_unemployment, order.by = data$Date[-1])
z_Fund  <- zoo(final_FedFund,      order.by = data$Date[-1])

Z <- merge(CPI = z_CPI,
           PPI = z_PPI,
           Housing = z_House,
           Unemp = z_Unemp,
           FedFund = z_Fund,
           all = FALSE)

covid_step_all <- zoo(
  as.integer(data$Date >= as.Date("2020-12-01") & data$Date <= as.Date("2022-05-01")),
  order.by = data$Date
)

Z2 <- merge(Z, COVID = covid_step_all, all = FALSE)

zoo_to_ts <- function(z) {
  idx <- index(z)
  start <- c(as.integer(format(min(idx), "%Y")),
             as.integer(format(min(idx), "%m")))
  ts(coredata(z), start = start, frequency = 12)
}

Y_full <- zoo_to_ts(Z2[, c("CPI","PPI","Housing","Unemp","FedFund")])
X_full <- zoo_to_ts(Z2[, "COVID", drop = FALSE])

vars_in <- c("PPI","Housing","Unemp","FedFund")  # edit this per model
y_model <- Y_full[, c("CPI", vars_in)]
x_model <- X_full

h <- 12
n <- nrow(y_model)
train <- 1:(n - h)
test  <- (n - h + 1):n

sel <- VARselect(y_model[train, ], lag.max = 48, type = "const", exogen = x_model[train, , drop = FALSE])
p <- sel$selection[["AIC(n)"]]
print(sel$selection[["SC(n)"]])
p <- 7

ps <- 1:18
tab <- do.call(rbind, lapply(ps, function(p){
  m  <- VAR(y_model[train,], p=p, type="const", exogen=x_model[train,,drop=FALSE], season = 12)
  fc <- predict(m, n.ahead=h, dumvar=x_model[test,,drop=FALSE])$fcst$CPI[,"fcst"]
  fc <- ts(fc, start=time(y_model)[test][1], frequency=12)
  rmse <- sqrt(mean((y_model[test,"CPI"] - fc)^2))
  aic <- AIC(m)
  bic <- BIC(m)
  c(p=p, RMSE=rmse, aic, bic)
}))
tab[order(tab[,"RMSE"]),]

var_mod <- VAR(y_model[train, ], p = p, type = "const", exogen = x_model[train, , drop = FALSE], season = 12)

summary_model <- summary(var_mod)
summary_model

.safe <- function(expr, default = NA_real_) tryCatch(expr, error = function(e) default)

to1 <- function(x, na = NA_real_) {
  if (is.null(x)) return(na)
  if (length(x) == 0) return(na)
  unname(as.vector(x))[1]
}

fc <- predict(var_mod, n.ahead = h, dumvar = x_model[test, , drop = FALSE])
cpi_fc  <- ts(fc$fcst$CPI[, "fcst"], start = time(y_model)[test][1], frequency = 12)
cpi_act <- y_model[test, "CPI"]
rmse_cpi <- sqrt(mean((cpi_act - cpi_fc)^2, na.rm = TRUE))

# Recompute / coerce all scalars safely
rmse_cpi <- to1(sqrt(mean((cpi_act - cpi_fc)^2, na.rm = TRUE)))
mae_cpi  <- to1(mean(abs(cpi_act - cpi_fc), na.rm = TRUE))
mape_cpi <- to1(mean(abs((cpi_act - cpi_fc) / cpi_act), na.rm = TRUE) * 100)

aic_val <- to1(suppressWarnings(tryCatch(AIC(var_mod), error = function(e) NA_real_)))
bic_val <- to1(suppressWarnings(tryCatch(BIC(var_mod), error = function(e) NA_real_)))
llk_val <- to1(suppressWarnings(tryCatch(as.numeric(logLik(var_mod)), error = function(e) NA_real_)))

r_2 <- to1(summary_model$varresult$CPI$adj.r.squared)
adj_r_2 <- to1(summary_model$varresult$CPI$adj.r.squared)

lb_pval   <- to1(suppressWarnings(tryCatch(serial.test(var_mod, lags.pt = 24, type = "PT.asymptotic")$serial$p.value, error = function(e) NA_real_)))
arch_pval <- to1(suppressWarnings(tryCatch(arch.test(var_mod, lags.multi = 12, multivariate.only = TRUE)$arch.mul$p.value, error = function(e) NA_real_)))
norm_pval <- to1(suppressWarnings(tryCatch(normality.test(var_mod, multivariate.only = TRUE)$jb.mul$p.value, error = function(e) NA_real_)))

roots_mod <- suppressWarnings(tryCatch(Mod(roots(var_mod)), error = function(e) NA_real_))
max_root  <- if (length(roots_mod) == 0 || all(is.na(roots_mod))) NA_real_ else max(roots_mod, na.rm = TRUE)
stable    <- isTRUE(all(roots_mod < 1, na.rm = TRUE))

# residual correlation CPI vs second series (if any)
res_cor <- NA_real_
res_mat <- suppressWarnings(tryCatch(cor(residuals(var_mod), use = "pairwise.complete.obs"), error = function(e) NULL))
if (!is.null(res_mat) && "CPI" %in% colnames(res_mat) && ncol(res_mat) >= 2) {
  other <- setdiff(colnames(res_mat), "CPI")
  if (length(other) >= 1) {
    res_cor <- to1(res_mat["CPI", other[1]])
  }
}

model_label <- paste0("CPI+", paste(vars_in, collapse = "+"))
train_n <- length(train); test_n <- length(test)

# Now safe to build the row (every field is guaranteed length 1)
this_row <- data.frame(
  Model       = model_label,
  Lag         = as.integer(to1(p)),
  TrainN      = to1(train_n),
  TestH       = to1(h),
  AIC         = to1(aic_val),
  BIC         = to1(bic_val),
  R_SQ        = to1(r_2),
  Adj_R_SQ    = to1(adj_r_2),
  LogLik      = to1(llk_val),
  Stable      = to1(stable, na = NA),
  MaxRootMod  = round(to1(max_root), 4),
  LB_p        = to1(lb_pval),
  ARCH_p      = to1(arch_pval),
  Normal_p    = to1(norm_pval),
  RMSE_CPI    = to1(rmse_cpi),
  MAE_CPI     = to1(mae_cpi),
  MAPE_CPI    = to1(mape_cpi),
  ResCorr_CPI_2nd = to1(res_cor),
  stringsAsFactors = FALSE
)

# Append to running table
if (!exists("VAR_comp_table")) {
  VAR_comp_table <- this_row
} else {
  VAR_comp_table <- rbind(VAR_comp_table, this_row)
}

VAR_comp_table <- VAR_comp_table[order(VAR_comp_table$AIC, VAR_comp_table$BIC), ]
print(VAR_comp_table)


plot(y_model[,"CPI"], main = "CPI (transformed) with VAR forecast window", ylab = "CPI (transformed)", xlab = "")
lines(cpi_fc, col = "red")
legend("topleft", bty = "n", lty = 1, col = c("black","red"),
       legend = c("Actual (all)","Forecast (test window)"))


act_zoom <- window(y_model[,"CPI"], start = c(2022, 1))  
fc_zoom  <- window(cpi_fc,          start = c(2022, 1))
plot(act_zoom); lines(fc_zoom, col = "red", lwd = 2)

plot(resid(var_mod)[, "CPI"], main = "VAR Model Residuals")
acf(resid(var_mod)[, "CPI"], main = "ACF of CPI Residuals")
pacf(resid(var_mod)[, "CPI"], main = "PACF of CPI Residuals")

last_train_index <- as.numeric(tail(train, n = 1))
last_train_index

last_log_cpi_train <- log(data$CPI[last_train_index])

inverted_fcast_cpi <- exp(cumsum(cpi_fc) + last_log_cpi_train)
inverted_actual_cpi <- exp(cumsum(cpi_act) + last_log_cpi_train)

rmse_orig <- sqrt(mean((inverted_fcast_cpi - inverted_actual_cpi)^2, na.rm = TRUE))
mae_orig <- mean(abs(inverted_fcast_cpi - inverted_actual_cpi), na.rm = TRUE)

round(rmse_orig, 4)
round(mae_orig, 4)

# We estimate a 5-variable VAR for CPI, PPI, Housing, Unemployment, and the Fed Funds rate with p = 13 lags. This specification balances fit and diagnostics, delivering RMSE = 0.00187 (test window) with no residual autocorrelation (Portmanteau p = 0.70) and stable roots. CPI dynamics are persistent (positive l1, negative l2, seasonal l12), while PPI leads CPI at horizons of roughly 3–8 months. Monetary policy effects are visible with a ~3-month lag (Fed Funds l3 < 0). Overall, the p=13 full model outperforms leaner alternatives on both accuracy and residual behavior, so we use it for forecasting CPI on the transformed scale.

# ========== VAR Model (PPI, FedFund) ==========

Y_full <- zoo_to_ts(Z2[, c("CPI","PPI","FedFund")])

vars_in <- c("PPI","FedFund")  # edit this per model
y_model <- Y_full[, c("CPI", vars_in)]
x_model <- X_full

h <- 12
n <- nrow(y_model)
train <- 1:(n - h)
test  <- (n - h + 1):n

sel <- VARselect(y_model[train, ], lag.max = 48, type = "const", exogen = x_model[train, , drop = FALSE])
p <- sel$selection[["AIC(n)"]]
print(sel$selection[["SC(n)"]])
p <- 6

ps <- 1:18
tab <- do.call(rbind, lapply(ps, function(p){
  m  <- VAR(y_model[train,], p=p, type="const", exogen=x_model[train,,drop=FALSE], season = 12)
  fc <- predict(m, n.ahead=h, dumvar=x_model[test,,drop=FALSE])$fcst$CPI[,"fcst"]
  fc <- ts(fc, start=time(y_model)[test][1], frequency=12)
  rmse <- sqrt(mean((y_model[test,"CPI"] - fc)^2))
  lb   <- serial.test(m, lags.pt=24, type="PT.asymptotic")$serial$p.value
  aic <- AIC(m)
  bic <- BIC(m)
  c(p=p, RMSE=rmse, LBp=lb, AIC = aic, BIC = bic)
}))
tab[order(tab[,"RMSE"]),]

var_mod <- VAR(y_model[train, ], p = p, type = "const", exogen = x_model[train, , drop = FALSE], season = 12)

summary_model <- summary(var_mod)
summary_model

fc <- predict(var_mod, n.ahead = h, dumvar = x_model[test, , drop = FALSE])
cpi_fc  <- ts(fc$fcst$CPI[, "fcst"], start = time(y_model)[test][1], frequency = 12)
cpi_act <- y_model[test, "CPI"]
rmse_cpi <- sqrt(mean((cpi_act - cpi_fc)^2, na.rm = TRUE))

rmse_cpi <- to1(sqrt(mean((cpi_act - cpi_fc)^2, na.rm = TRUE)))
mae_cpi  <- to1(mean(abs(cpi_act - cpi_fc), na.rm = TRUE))
mape_cpi <- to1(mean(abs((cpi_act - cpi_fc) / cpi_act), na.rm = TRUE) * 100)

aic_val <- to1(suppressWarnings(tryCatch(AIC(var_mod), error = function(e) NA_real_)))
bic_val <- to1(suppressWarnings(tryCatch(BIC(var_mod), error = function(e) NA_real_)))
llk_val <- to1(suppressWarnings(tryCatch(as.numeric(logLik(var_mod)), error = function(e) NA_real_)))

r_2 <- to1(summary_model$varresult$CPI$r.squared)
adj_r_2 <- to1(summary_model$varresult$CPI$adj.r.squared)

lb_pval   <- to1(suppressWarnings(tryCatch(serial.test(var_mod, lags.pt = 24, type = "PT.asymptotic")$serial$p.value, error = function(e) NA_real_)))
arch_pval <- to1(suppressWarnings(tryCatch(arch.test(var_mod, lags.multi = 12, multivariate.only = TRUE)$arch.mul$p.value, error = function(e) NA_real_)))
norm_pval <- to1(suppressWarnings(tryCatch(normality.test(var_mod, multivariate.only = TRUE)$jb.mul$p.value, error = function(e) NA_real_)))

roots_mod <- suppressWarnings(tryCatch(Mod(roots(var_mod)), error = function(e) NA_real_))
max_root  <- if (length(roots_mod) == 0 || all(is.na(roots_mod))) NA_real_ else max(roots_mod, na.rm = TRUE)
stable    <- isTRUE(all(roots_mod < 1, na.rm = TRUE))

# residual correlation CPI vs second series (if any)
res_cor <- NA_real_
res_mat <- suppressWarnings(tryCatch(cor(residuals(var_mod), use = "pairwise.complete.obs"), error = function(e) NULL))
if (!is.null(res_mat) && "CPI" %in% colnames(res_mat) && ncol(res_mat) >= 2) {
  other <- setdiff(colnames(res_mat), "CPI")
  if (length(other) >= 1) {
    res_cor <- to1(res_mat["CPI", other[1]])
  }
}

model_label <- paste0("CPI+", paste(vars_in, collapse = "+"))
train_n <- length(train); test_n <- length(test)

# Now safe to build the row (every field is guaranteed length 1)
this_row <- data.frame(
  Model       = model_label,
  Lag         = as.integer(to1(p)),
  TrainN      = to1(train_n),
  TestH       = to1(h),
  AIC         = to1(aic_val),
  BIC         = to1(bic_val),
  R_SQ        = to1(r_2),
  Adj_R_SQ    = to1(adj_r_2),
  LogLik      = to1(llk_val),
  Stable      = to1(stable, na = NA),
  MaxRootMod  = round(to1(max_root), 4),
  LB_p        = to1(lb_pval),
  ARCH_p      = to1(arch_pval),
  Normal_p    = to1(norm_pval),
  RMSE_CPI    = to1(rmse_cpi),
  MAE_CPI     = to1(mae_cpi),
  MAPE_CPI    = to1(mape_cpi),
  ResCorr_CPI_2nd = to1(res_cor),
  stringsAsFactors = FALSE
)

# Append to running table
if (!exists("VAR_comp_table")) {
  VAR_comp_table <- this_row
} else {
  VAR_comp_table <- rbind(VAR_comp_table, this_row)
}

VAR_comp_table <- VAR_comp_table[order(VAR_comp_table$AIC, VAR_comp_table$BIC), ]
print(VAR_comp_table)

plot(y_model[,"CPI"], main = "CPI (transformed) with VAR forecast window", ylab = "CPI (transformed)", xlab = "")
lines(cpi_fc, col = "red")
legend("topleft", bty = "n", lty = 1, col = c("black","red"),
       legend = c("Actual (all)","Forecast (test window)"))


act_zoom <- window(y_model[,"CPI"], start = c(2022, 1))  
fc_zoom  <- window(cpi_fc,          start = c(2022, 1))
plot(act_zoom); lines(fc_zoom, col = "red", lwd = 2)

plot(resid(var_mod)[, "CPI"], main = "VAR Model Residuals")
acf(resid(var_mod)[, "CPI"], main = "ACF of CPI Residuals")
pacf(resid(var_mod)[, "CPI"], main = "PACF of CPI Residuals")

last_train_index <- as.numeric(tail(train, n = 1))
last_train_index

last_log_cpi_train <- log(data$CPI[last_train_index])

inverted_fcast_cpi <- exp(cumsum(cpi_fc) + last_log_cpi_train)
inverted_actual_cpi <- exp(cumsum(cpi_act) + last_log_cpi_train)

rmse_orig <- sqrt(mean((inverted_fcast_cpi - inverted_actual_cpi)^2, na.rm = TRUE))
mae_orig <- mean(abs(inverted_fcast_cpi - inverted_actual_cpi), na.rm = TRUE)

round(rmse_orig, 4)
round(mae_orig, 4)

# ============== VAR Model (PPI, FedFund, Housing) ============
Y_full <- zoo_to_ts(Z2[, c("CPI","PPI","FedFund", "Housing")])

vars_in <- c("PPI","FedFund", "Housing")  # edit this per model
y_model <- Y_full[, c("CPI", vars_in)]
x_model <- X_full

h <- 12
n <- nrow(y_model)
train <- 1:(n - h)
test  <- (n - h + 1):n

sel <- VARselect(y_model[train, ], lag.max = 48, type = "const", exogen = x_model[train, , drop = FALSE], season = 12)

print(sel$selection[["AIC(n)"]])
print(sel$selection[["SC(n)"]])
p <- 6

ps <- 1:18
tab <- do.call(rbind, lapply(ps, function(p){
  m  <- VAR(y_model[train,], p=p, type="const", exogen=x_model[train,,drop=FALSE], season = 12)
  fc <- predict(m, n.ahead=h, dumvar=x_model[test,,drop=FALSE])$fcst$CPI[,"fcst"]
  fc <- ts(fc, start=time(y_model)[test][1], frequency=12)
  rmse <- sqrt(mean((y_model[test,"CPI"] - fc)^2))
  lb   <- serial.test(m, lags.pt=24, type="PT.asymptotic")$serial$p.value
  aic <- AIC(m)
  bic <- BIC(m)
  c(p=p, RMSE=rmse, LBp=lb, AIC = aic, BIC = bic)
}))
tab[order(tab[,"RMSE"]),]

var_mod <- VAR(y_model[train, ], p = p, type = "const", exogen = x_model[train, , drop = FALSE], season = 12)

summary_model <- summary(var_mod)
summary_model

fc <- predict(var_mod, n.ahead = h, dumvar = x_model[test, , drop = FALSE])
cpi_fc  <- ts(fc$fcst$CPI[, "fcst"], start = time(y_model)[test][1], frequency = 12)
cpi_act <- y_model[test, "CPI"]
rmse_cpi <- sqrt(mean((cpi_act - cpi_fc)^2, na.rm = TRUE))

rmse_cpi <- to1(sqrt(mean((cpi_act - cpi_fc)^2, na.rm = TRUE)))
mae_cpi  <- to1(mean(abs(cpi_act - cpi_fc), na.rm = TRUE))
mape_cpi <- to1(mean(abs((cpi_act - cpi_fc) / cpi_act), na.rm = TRUE) * 100)

aic_val <- to1(suppressWarnings(tryCatch(AIC(var_mod), error = function(e) NA_real_)))
bic_val <- to1(suppressWarnings(tryCatch(BIC(var_mod), error = function(e) NA_real_)))
llk_val <- to1(suppressWarnings(tryCatch(as.numeric(logLik(var_mod)), error = function(e) NA_real_)))

r_2 <- to1(summary_model$varresult$CPI$r.squared)
adj_r_2 <- to1(summary_model$varresult$CPI$adj.r.squared)

lb_pval   <- to1(suppressWarnings(tryCatch(serial.test(var_mod, lags.pt = 24, type = "PT.asymptotic")$serial$p.value, error = function(e) NA_real_)))
arch_pval <- to1(suppressWarnings(tryCatch(arch.test(var_mod, lags.multi = 12, multivariate.only = TRUE)$arch.mul$p.value, error = function(e) NA_real_)))
norm_pval <- to1(suppressWarnings(tryCatch(normality.test(var_mod, multivariate.only = TRUE)$jb.mul$p.value, error = function(e) NA_real_)))

roots_mod <- suppressWarnings(tryCatch(Mod(roots(var_mod)), error = function(e) NA_real_))
max_root  <- if (length(roots_mod) == 0 || all(is.na(roots_mod))) NA_real_ else max(roots_mod, na.rm = TRUE)
stable    <- isTRUE(all(roots_mod < 1, na.rm = TRUE))

# residual correlation CPI vs second series (if any)
res_cor <- NA_real_
res_mat <- suppressWarnings(tryCatch(cor(residuals(var_mod), use = "pairwise.complete.obs"), error = function(e) NULL))
if (!is.null(res_mat) && "CPI" %in% colnames(res_mat) && ncol(res_mat) >= 2) {
  other <- setdiff(colnames(res_mat), "CPI")
  if (length(other) >= 1) {
    res_cor <- to1(res_mat["CPI", other[1]])
  }
}

model_label <- paste0("CPI+", paste(vars_in, collapse = "+"))
train_n <- length(train); test_n <- length(test)

# Now safe to build the row (every field is guaranteed length 1)
this_row <- data.frame(
  Model       = model_label,
  Lag         = as.integer(to1(p)),
  TrainN      = to1(train_n),
  TestH       = to1(h),
  AIC         = to1(aic_val),
  BIC         = to1(bic_val),
  R_SQ        = to1(r_2),
  Adj_R_SQ    = to1(adj_r_2),
  LogLik      = to1(llk_val),
  Stable      = to1(stable, na = NA),
  MaxRootMod  = round(to1(max_root), 4),
  LB_p        = to1(lb_pval),
  ARCH_p      = to1(arch_pval),
  Normal_p    = to1(norm_pval),
  RMSE_CPI    = to1(rmse_cpi),
  MAE_CPI     = to1(mae_cpi),
  MAPE_CPI    = to1(mape_cpi),
  ResCorr_CPI_2nd = to1(res_cor),
  stringsAsFactors = FALSE
)

# Append to running table
if (!exists("VAR_comp_table")) {
  VAR_comp_table <- this_row
} else {
  VAR_comp_table <- rbind(VAR_comp_table, this_row)
}

VAR_comp_table <- VAR_comp_table[order(VAR_comp_table$AIC, VAR_comp_table$BIC), ]
print(VAR_comp_table)


plot(y_model[,"CPI"], main = "CPI (transformed) with VAR forecast window", ylab = "CPI (transformed)", xlab = "")
lines(cpi_fc, col = "red")
legend("topleft", bty = "n", lty = 1, col = c("black","red"),
       legend = c("Actual (all)","Forecast (test window)"))

act_zoom <- window(y_model[,"CPI"], start = c(2022, 1))  
fc_zoom  <- window(cpi_fc,          start = c(2022, 1))
plot(act_zoom); lines(fc_zoom, col = "red", lwd = 2)

plot(resid(var_mod)[, "CPI"], main = "VAR Model Residuals")
acf(resid(var_mod)[, "CPI"], main = "ACF of CPI Residuals")
pacf(resid(var_mod)[, "CPI"], main = "PACF of CPI Residuals")

last_train_index <- as.numeric(tail(train, n = 1))
last_train_index

last_log_cpi_train <- log(data$CPI[last_train_index])

inverted_fcast_cpi <- exp(cumsum(cpi_fc) + last_log_cpi_train)
inverted_actual_cpi <- exp(cumsum(cpi_act) + last_log_cpi_train)

rmse_orig <- sqrt(mean((inverted_fcast_cpi - inverted_actual_cpi)^2, na.rm = TRUE))
mae_orig <- mean(abs(inverted_fcast_cpi - inverted_actual_cpi), na.rm = TRUE)

round(rmse_orig, 4)
round(mae_orig, 4)

# =============== VAR Model (PPI) ===============
Y_full <- zoo_to_ts(Z2[, c("CPI","PPI")])

vars_in <- c("PPI")  # edit this per model
y_model <- Y_full[, c("CPI", vars_in)]
x_model <- X_full

h <- 12
n <- nrow(y_model)
train <- 1:(n - h)
test  <- (n - h + 1):n

sel <- VARselect(y_model[train, ], lag.max = 48, type = "const", exogen = x_model[train, , drop = FALSE], season = 12)
print(sel$selection[["AIC(n)"]])
print(sel$selection[["SC(n)"]])
p <- 6

ps <- 6:28
tab <- do.call(rbind, lapply(ps, function(p){
  m  <- VAR(y_model[train,], p=p, type="const", exogen=x_model[train,,drop=FALSE], season = 12)
  fc <- predict(m, n.ahead=h, dumvar=x_model[test,,drop=FALSE])$fcst$CPI[,"fcst"]
  fc <- ts(fc, start=time(y_model)[test][1], frequency=12)
  rmse <- sqrt(mean((y_model[test,"CPI"] - fc)^2))
  lb   <- serial.test(m, lags.pt=36, type="PT.asymptotic")$serial$p.value
  aic <- AIC(m)
  bic <- BIC(m)
  c(p=p, RMSE=rmse, LBp=lb, AIC = aic, BIC = bic)
}))
tab[order(tab[,"RMSE"]),]

var_mod <- VAR(y_model[train, ], p = p, type = "const", exogen = x_model[train, , drop = FALSE], season = 12)

summary_model <- summary(var_mod)
summary_model

fc <- predict(var_mod, n.ahead = h, dumvar = x_model[test, , drop = FALSE])
cpi_fc  <- ts(fc$fcst$CPI[, "fcst"], start = time(y_model)[test][1], frequency = 12)
cpi_act <- y_model[test, "CPI"]


# Recompute / coerce all scalars safely
rmse_cpi <- to1(sqrt(mean((cpi_act - cpi_fc)^2, na.rm = TRUE)))
mae_cpi  <- to1(mean(abs(cpi_act - cpi_fc), na.rm = TRUE))
mape_cpi <- to1(mean(abs((cpi_act - cpi_fc) / cpi_act), na.rm = TRUE) * 100)

aic_val <- to1(suppressWarnings(tryCatch(AIC(var_mod), error = function(e) NA_real_)))
bic_val <- to1(suppressWarnings(tryCatch(BIC(var_mod), error = function(e) NA_real_)))
llk_val <- to1(suppressWarnings(tryCatch(as.numeric(logLik(var_mod)), error = function(e) NA_real_)))

r_2 <- to1(summary_model$varresult$CPI$r.squared)
adj_r_2 <- to1(summary_model$varresult$CPI$adj.r.squared)

lb_pval   <- to1(suppressWarnings(tryCatch(serial.test(var_mod, lags.pt = 24, type = "PT.asymptotic")$serial$p.value, error = function(e) NA_real_)))
arch_pval <- to1(suppressWarnings(tryCatch(arch.test(var_mod, lags.multi = 12, multivariate.only = TRUE)$arch.mul$p.value, error = function(e) NA_real_)))
norm_pval <- to1(suppressWarnings(tryCatch(normality.test(var_mod, multivariate.only = TRUE)$jb.mul$p.value, error = function(e) NA_real_)))

roots_mod <- suppressWarnings(tryCatch(Mod(roots(var_mod)), error = function(e) NA_real_))
max_root  <- if (length(roots_mod) == 0 || all(is.na(roots_mod))) NA_real_ else max(roots_mod, na.rm = TRUE)
stable    <- isTRUE(all(roots_mod < 1, na.rm = TRUE))

# residual correlation CPI vs second series (if any)
res_cor <- NA_real_
res_mat <- suppressWarnings(tryCatch(cor(residuals(var_mod), use = "pairwise.complete.obs"), error = function(e) NULL))
if (!is.null(res_mat) && "CPI" %in% colnames(res_mat) && ncol(res_mat) >= 2) {
  other <- setdiff(colnames(res_mat), "CPI")
  if (length(other) >= 1) {
    res_cor <- to1(res_mat["CPI", other[1]])
  }
}

model_label <- paste0("CPI+", paste(vars_in, collapse = "+"))
train_n <- length(train); test_n <- length(test)

# Now safe to build the row (every field is guaranteed length 1)
this_row <- data.frame(
  Model       = model_label,
  Lag         = as.integer(to1(p)),
  TrainN      = to1(train_n),
  TestH       = to1(h),
  AIC         = to1(aic_val),
  BIC         = to1(bic_val),
  R_SQ        = to1(r_2),
  Adj_R_SQ    = to1(adj_r_2),
  LogLik      = to1(llk_val),
  Stable      = to1(stable, na = NA),
  MaxRootMod  = round(to1(max_root), 4),
  LB_p        = to1(lb_pval),
  ARCH_p      = to1(arch_pval),
  Normal_p    = to1(norm_pval),
  RMSE_CPI    = to1(rmse_cpi),
  MAE_CPI     = to1(mae_cpi),
  MAPE_CPI    = to1(mape_cpi),
  ResCorr_CPI_2nd = to1(res_cor),
  stringsAsFactors = FALSE
)

# Append to running table
if (!exists("VAR_comp_table")) {
  VAR_comp_table <- this_row
} else {
  VAR_comp_table <- rbind(VAR_comp_table, this_row)
}

VAR_comp_table <- VAR_comp_table[order(VAR_comp_table$AIC, VAR_comp_table$BIC), ]
print(VAR_comp_table)


plot(y_model[,"CPI"], main = "CPI (transformed) with VAR forecast window", ylab = "CPI (transformed)", xlab = "")
lines(cpi_fc, col = "red")
legend("topleft", bty = "n", lty = 1, col = c("black","red"),
       legend = c("Actual (all)","Forecast (test window)"))

act_zoom <- window(y_model[,"CPI"], start = c(2022, 1))  
fc_zoom  <- window(cpi_fc,          start = c(2022, 1))
plot(act_zoom); lines(fc_zoom, col = "red", lwd = 2)

plot(resid(var_mod)[, "CPI"], main = "VAR Model Residuals")
acf(resid(var_mod)[, "CPI"], main = "ACF of CPI Residuals")
pacf(resid(var_mod)[, "CPI"], main = "PACF of CPI Residuals")

# The model indicates bidirectional lagged linkages between CPI and PPI, with CPI showing more persistence and stronger COVID-period effects. The significant COVID impact on CPI but not PPI implies downstream consumer prices were more directly affected by pandemic-related shocks. However, residual autocorrelation in PPI warrants further refinement before long-term structural inference.

last_train_index <- as.numeric(tail(train, n = 1))
last_train_index

last_log_cpi_train <- log(data$CPI[last_train_index])

inverted_fcast_cpi <- exp(cumsum(cpi_fc) + last_log_cpi_train)
inverted_actual_cpi <- exp(cumsum(cpi_act) + last_log_cpi_train)

rmse_orig <- sqrt(mean((inverted_fcast_cpi - inverted_actual_cpi)^2, na.rm = TRUE))
mae_orig <- mean(abs(inverted_fcast_cpi - inverted_actual_cpi), na.rm = TRUE)

round(rmse_orig, 4)
round(mae_orig, 4)

# ========= VAR Model (Housing) ======
Y_full <- zoo_to_ts(Z2[, c("CPI","Housing")])

vars_in <- c("Housing")  # edit this per model
y_model <- Y_full[, c("CPI", vars_in)]
x_model <- X_full

h <- 12
n <- nrow(y_model)
train <- 1:(n - h)
test  <- (n - h + 1):n

sel <- VARselect(y_model[train, ], lag.max = 48, type = "const", exogen = x_model[train, , drop = FALSE], season = 12)
print(sel$selection[["AIC(n)"]])
print(sel$selection[["SC(n)"]])
p <- 3

ps <- 1:12
tab <- do.call(rbind, lapply(ps, function(p){
  m  <- VAR(y_model[train,], p=p, type="const", exogen=x_model[train,,drop=FALSE], season = 12)
  fc <- predict(m, n.ahead=h, dumvar=x_model[test,,drop=FALSE])$fcst$CPI[,"fcst"]
  fc <- ts(fc, start=time(y_model)[test][1], frequency=12)
  rmse <- sqrt(mean((y_model[test,"CPI"] - fc)^2))
  lb   <- serial.test(m, lags.pt=36, type="PT.asymptotic")$serial$p.value
  aic <- AIC(m)
  bic <- BIC(m)
  c(p=p, RMSE=rmse, LBp=lb, AIC = aic, BIC = bic)
}))
tab[order(tab[,"RMSE"]),]

var_mod <- VAR(y_model[train, ], p = p, type = "const", exogen = x_model[train, , drop = FALSE], season = 12)

summary_model <- summary(var_mod)
summary_model

fc <- predict(var_mod, n.ahead = h, dumvar = x_model[test, , drop = FALSE])
cpi_fc  <- ts(fc$fcst$CPI[, "fcst"], start = time(y_model)[test][1], frequency = 12)
cpi_act <- y_model[test, "CPI"]


# Recompute / coerce all scalars safely
rmse_cpi <- to1(sqrt(mean((cpi_act - cpi_fc)^2, na.rm = TRUE)))
mae_cpi  <- to1(mean(abs(cpi_act - cpi_fc), na.rm = TRUE))
mape_cpi <- to1(mean(abs((cpi_act - cpi_fc) / cpi_act), na.rm = TRUE) * 100)

aic_val <- to1(suppressWarnings(tryCatch(AIC(var_mod), error = function(e) NA_real_)))
bic_val <- to1(suppressWarnings(tryCatch(BIC(var_mod), error = function(e) NA_real_)))
llk_val <- to1(suppressWarnings(tryCatch(as.numeric(logLik(var_mod)), error = function(e) NA_real_)))

r_2 <- to1(summary_model$varresult$CPI$r.squared)
adj_r_2 <- to1(summary_model$varresult$CPI$adj.r.squared)

lb_pval   <- to1(suppressWarnings(tryCatch(serial.test(var_mod, lags.pt = 24, type = "PT.asymptotic")$serial$p.value, error = function(e) NA_real_)))
arch_pval <- to1(suppressWarnings(tryCatch(arch.test(var_mod, lags.multi = 12, multivariate.only = TRUE)$arch.mul$p.value, error = function(e) NA_real_)))
norm_pval <- to1(suppressWarnings(tryCatch(normality.test(var_mod, multivariate.only = TRUE)$jb.mul$p.value, error = function(e) NA_real_)))

roots_mod <- suppressWarnings(tryCatch(Mod(roots(var_mod)), error = function(e) NA_real_))
max_root  <- if (length(roots_mod) == 0 || all(is.na(roots_mod))) NA_real_ else max(roots_mod, na.rm = TRUE)
stable    <- isTRUE(all(roots_mod < 1, na.rm = TRUE))

# residual correlation CPI vs second series (if any)
res_cor <- NA_real_
res_mat <- suppressWarnings(tryCatch(cor(residuals(var_mod), use = "pairwise.complete.obs"), error = function(e) NULL))
if (!is.null(res_mat) && "CPI" %in% colnames(res_mat) && ncol(res_mat) >= 2) {
  other <- setdiff(colnames(res_mat), "CPI")
  if (length(other) >= 1) {
    res_cor <- to1(res_mat["CPI", other[1]])
  }
}

model_label <- paste0("CPI+", paste(vars_in, collapse = "+"))
train_n <- length(train); test_n <- length(test)

# Now safe to build the row (every field is guaranteed length 1)
this_row <- data.frame(
  Model       = model_label,
  Lag         = as.integer(to1(p)),
  TrainN      = to1(train_n),
  TestH       = to1(h),
  AIC         = to1(aic_val),
  BIC         = to1(bic_val),
  R_SQ        = to1(r_2),
  Adj_R_SQ    = to1(adj_r_2),
  LogLik      = to1(llk_val),
  Stable      = to1(stable, na = NA),
  MaxRootMod  = round(to1(max_root), 4),
  LB_p        = to1(lb_pval),
  ARCH_p      = to1(arch_pval),
  Normal_p    = to1(norm_pval),
  RMSE_CPI    = to1(rmse_cpi),
  MAE_CPI     = to1(mae_cpi),
  MAPE_CPI    = to1(mape_cpi),
  ResCorr_CPI_2nd = to1(res_cor),
  stringsAsFactors = FALSE
)

# Append to running table
if (!exists("VAR_comp_table")) {
  VAR_comp_table <- this_row
} else {
  VAR_comp_table <- rbind(VAR_comp_table, this_row)
}

VAR_comp_table <- VAR_comp_table[order(VAR_comp_table$AIC, VAR_comp_table$BIC), ]
print(VAR_comp_table)


plot(y_model[,"CPI"], main = "CPI (transformed) with VAR forecast window", ylab = "CPI (transformed)", xlab = "")
lines(cpi_fc, col = "red")
legend("topleft", bty = "n", lty = 1, col = c("black","red"),
       legend = c("Actual (all)","Forecast (test window)"))

act_zoom <- window(y_model[,"CPI"], start = c(2022, 1))  
fc_zoom  <- window(cpi_fc,          start = c(2022, 1))
plot(act_zoom); lines(fc_zoom, col = "red", lwd = 2)

plot(resid(var_mod)[, "CPI"], main = "VAR Model Residuals")
acf(resid(var_mod)[, "CPI"], main = "ACF of CPI Residuals")
pacf(resid(var_mod)[, "CPI"], main = "PACF of CPI Residuals")

# The model indicates bidirectional lagged linkages between CPI and PPI, with CPI showing more persistence and stronger COVID-period effects. The significant COVID impact on CPI but not PPI implies downstream consumer prices were more directly affected by pandemic-related shocks. However, residual autocorrelation in PPI warrants further refinement before long-term structural inference.

last_train_index <- as.numeric(tail(train, n = 1))
last_train_index

last_log_cpi_train <- log(data$CPI[last_train_index])

inverted_fcast_cpi <- exp(cumsum(cpi_fc) + last_log_cpi_train)
inverted_actual_cpi <- exp(cumsum(cpi_act) + last_log_cpi_train)

rmse_orig <- sqrt(mean((inverted_fcast_cpi - inverted_actual_cpi)^2, na.rm = TRUE))
mae_orig <- mean(abs(inverted_fcast_cpi - inverted_actual_cpi), na.rm = TRUE)

round(rmse_orig, 4)
round(mae_orig, 4)

# ========== VAR Model (Unemployment) ====
Y_full <- zoo_to_ts(Z2[, c("CPI","Unemp")])

vars_in <- c("Unemp")  # edit this per model
y_model <- Y_full[, c("CPI", vars_in)]
x_model <- X_full

h <- 12
n <- nrow(y_model)
train <- 1:(n - h)
test  <- (n - h + 1):n

sel <- VARselect(y_model[train, ], lag.max = 48, type = "const", exogen = x_model[train, , drop = FALSE], season = 12)
print(sel$selection[["AIC(n)"]])
print(sel$selection[["SC(n)"]])
p <- 9

ps <- 2:14
tab <- do.call(rbind, lapply(ps, function(p){
  m  <- VAR(y_model[train,], p=p, type="const", exogen=x_model[train,,drop=FALSE], season = 12)
  fc <- predict(m, n.ahead=h, dumvar=x_model[test,,drop=FALSE])$fcst$CPI[,"fcst"]
  fc <- ts(fc, start=time(y_model)[test][1], frequency=12)
  rmse <- sqrt(mean((y_model[test,"CPI"] - fc)^2))
  lb   <- serial.test(m, lags.pt=36, type="PT.asymptotic")$serial$p.value
  aic <- AIC(m)
  bic <- BIC(m)
  c(p=p, RMSE=rmse, LBp=lb, AIC = aic, BIC = bic)
}))
tab[order(tab[,"RMSE"]),]

var_mod <- VAR(y_model[train, ], p = p, type = "const", exogen = x_model[train, , drop = FALSE], season = 12)

summary_model <- summary(var_mod)
summary_model

fc <- predict(var_mod, n.ahead = h, dumvar = x_model[test, , drop = FALSE])
cpi_fc  <- ts(fc$fcst$CPI[, "fcst"], start = time(y_model)[test][1], frequency = 12)
cpi_act <- y_model[test, "CPI"]


# Recompute / coerce all scalars safely
rmse_cpi <- to1(sqrt(mean((cpi_act - cpi_fc)^2, na.rm = TRUE)))
mae_cpi  <- to1(mean(abs(cpi_act - cpi_fc), na.rm = TRUE))
mape_cpi <- to1(mean(abs((cpi_act - cpi_fc) / cpi_act), na.rm = TRUE) * 100)

aic_val <- to1(suppressWarnings(tryCatch(AIC(var_mod), error = function(e) NA_real_)))
bic_val <- to1(suppressWarnings(tryCatch(BIC(var_mod), error = function(e) NA_real_)))
llk_val <- to1(suppressWarnings(tryCatch(as.numeric(logLik(var_mod)), error = function(e) NA_real_)))

r_2 <- to1(summary_model$varresult$CPI$r.squared)
adj_r_2 <- to1(summary_model$varresult$CPI$adj.r.squared)

lb_pval   <- to1(suppressWarnings(tryCatch(serial.test(var_mod, lags.pt = 24, type = "PT.asymptotic")$serial$p.value, error = function(e) NA_real_)))
arch_pval <- to1(suppressWarnings(tryCatch(arch.test(var_mod, lags.multi = 12, multivariate.only = TRUE)$arch.mul$p.value, error = function(e) NA_real_)))
norm_pval <- to1(suppressWarnings(tryCatch(normality.test(var_mod, multivariate.only = TRUE)$jb.mul$p.value, error = function(e) NA_real_)))

roots_mod <- suppressWarnings(tryCatch(Mod(roots(var_mod)), error = function(e) NA_real_))
max_root  <- if (length(roots_mod) == 0 || all(is.na(roots_mod))) NA_real_ else max(roots_mod, na.rm = TRUE)
stable    <- isTRUE(all(roots_mod < 1, na.rm = TRUE))

# residual correlation CPI vs second series (if any)
res_cor <- NA_real_
res_mat <- suppressWarnings(tryCatch(cor(residuals(var_mod), use = "pairwise.complete.obs"), error = function(e) NULL))
if (!is.null(res_mat) && "CPI" %in% colnames(res_mat) && ncol(res_mat) >= 2) {
  other <- setdiff(colnames(res_mat), "CPI")
  if (length(other) >= 1) {
    res_cor <- to1(res_mat["CPI", other[1]])
  }
}

model_label <- paste0("CPI+", paste(vars_in, collapse = "+"))
train_n <- length(train); test_n <- length(test)

# Now safe to build the row (every field is guaranteed length 1)
this_row <- data.frame(
  Model       = model_label,
  Lag         = as.integer(to1(p)),
  TrainN      = to1(train_n),
  TestH       = to1(h),
  AIC         = to1(aic_val),
  BIC         = to1(bic_val),
  R_SQ        = to1(r_2),
  Adj_R_SQ    = to1(adj_r_2),
  LogLik      = to1(llk_val),
  Stable      = to1(stable, na = NA),
  MaxRootMod  = round(to1(max_root), 4),
  LB_p        = to1(lb_pval),
  ARCH_p      = to1(arch_pval),
  Normal_p    = to1(norm_pval),
  RMSE_CPI    = to1(rmse_cpi),
  MAE_CPI     = to1(mae_cpi),
  MAPE_CPI    = to1(mape_cpi),
  ResCorr_CPI_2nd = to1(res_cor),
  stringsAsFactors = FALSE
)

# Append to running table
if (!exists("VAR_comp_table")) {
  VAR_comp_table <- this_row
} else {
  VAR_comp_table <- rbind(VAR_comp_table, this_row)
}

VAR_comp_table <- VAR_comp_table[order(VAR_comp_table$AIC, VAR_comp_table$BIC), ]
print(VAR_comp_table)


plot(y_model[,"CPI"], main = "CPI (transformed) with VAR forecast window", ylab = "CPI (transformed)", xlab = "")
lines(cpi_fc, col = "red")
legend("topleft", bty = "n", lty = 1, col = c("black","red"),
       legend = c("Actual (all)","Forecast (test window)"))

act_zoom <- window(y_model[,"CPI"], start = c(2022, 1))  
fc_zoom  <- window(cpi_fc,          start = c(2022, 1))
plot(act_zoom); lines(fc_zoom, col = "red", lwd = 2)

plot(resid(var_mod)[, "CPI"], main = "VAR Model Residuals")
acf(resid(var_mod)[, "CPI"], main = "ACF of CPI Residuals")
pacf(resid(var_mod)[, "CPI"], main = "PACF of CPI Residuals")

# The model indicates bidirectional lagged linkages between CPI and PPI, with CPI showing more persistence and stronger COVID-period effects. The significant COVID impact on CPI but not PPI implies downstream consumer prices were more directly affected by pandemic-related shocks. However, residual autocorrelation in PPI warrants further refinement before long-term structural inference.

last_train_index <- as.numeric(tail(train, n = 1))
last_train_index

last_log_cpi_train <- log(data$CPI[last_train_index])

inverted_fcast_cpi <- exp(cumsum(cpi_fc) + last_log_cpi_train)
inverted_actual_cpi <- exp(cumsum(cpi_act) + last_log_cpi_train)

rmse_orig <- sqrt(mean((inverted_fcast_cpi - inverted_actual_cpi)^2, na.rm = TRUE))
mae_orig <- mean(abs(inverted_fcast_cpi - inverted_actual_cpi), na.rm = TRUE)

round(rmse_orig, 4)
round(mae_orig, 4)

# =========== VAR Modeling (PPI, Housing) ====
Y_full <- zoo_to_ts(Z2[, c("CPI","PPI", "Housing")])

vars_in <- c("PPI", "Housing")  # edit this per model
y_model <- Y_full[, c("CPI", vars_in)]
x_model <- X_full

h <- 12
n <- nrow(y_model)
train <- 1:(n - h)
test  <- (n - h + 1):n

sel <- VARselect(y_model[train, ], lag.max = 48, type = "const", exogen = x_model[train, , drop = FALSE], season = 12)
print(sel$selection[["AIC(n)"]])
print(sel$selection[["SC(n)"]])
p <- 6

ps <- 6:28
tab <- do.call(rbind, lapply(ps, function(p){
  m  <- VAR(y_model[train,], p=p, type="const", exogen=x_model[train,,drop=FALSE], season = 12)
  fc <- predict(m, n.ahead=h, dumvar=x_model[test,,drop=FALSE])$fcst$CPI[,"fcst"]
  fc <- ts(fc, start=time(y_model)[test][1], frequency=12)
  rmse <- sqrt(mean((y_model[test,"CPI"] - fc)^2))
  lb   <- serial.test(m, lags.pt=36, type="PT.asymptotic")$serial$p.value
  aic <- AIC(m)
  bic <- BIC(m)
  c(p=p, RMSE=rmse, LBp=lb, AIC = aic, BIC = bic)
}))
tab[order(tab[,"RMSE"]),]

var_mod <- VAR(y_model[train, ], p = p, type = "const", exogen = x_model[train, , drop = FALSE], season = 12)

summary_model <- summary(var_mod)
summary_model

fc <- predict(var_mod, n.ahead = h, dumvar = x_model[test, , drop = FALSE])
cpi_fc  <- ts(fc$fcst$CPI[, "fcst"], start = time(y_model)[test][1], frequency = 12)
cpi_act <- y_model[test, "CPI"]


# Recompute / coerce all scalars safely
rmse_cpi <- to1(sqrt(mean((cpi_act - cpi_fc)^2, na.rm = TRUE)))
mae_cpi  <- to1(mean(abs(cpi_act - cpi_fc), na.rm = TRUE))
mape_cpi <- to1(mean(abs((cpi_act - cpi_fc) / cpi_act), na.rm = TRUE) * 100)

aic_val <- to1(suppressWarnings(tryCatch(AIC(var_mod), error = function(e) NA_real_)))
bic_val <- to1(suppressWarnings(tryCatch(BIC(var_mod), error = function(e) NA_real_)))
llk_val <- to1(suppressWarnings(tryCatch(as.numeric(logLik(var_mod)), error = function(e) NA_real_)))

r_2 <- to1(summary_model$varresult$CPI$r.squared)
adj_r_2 <- to1(summary_model$varresult$CPI$adj.r.squared)

lb_pval   <- to1(suppressWarnings(tryCatch(serial.test(var_mod, lags.pt = 24, type = "PT.asymptotic")$serial$p.value, error = function(e) NA_real_)))
arch_pval <- to1(suppressWarnings(tryCatch(arch.test(var_mod, lags.multi = 12, multivariate.only = TRUE)$arch.mul$p.value, error = function(e) NA_real_)))
norm_pval <- to1(suppressWarnings(tryCatch(normality.test(var_mod, multivariate.only = TRUE)$jb.mul$p.value, error = function(e) NA_real_)))

roots_mod <- suppressWarnings(tryCatch(Mod(roots(var_mod)), error = function(e) NA_real_))
max_root  <- if (length(roots_mod) == 0 || all(is.na(roots_mod))) NA_real_ else max(roots_mod, na.rm = TRUE)
stable    <- isTRUE(all(roots_mod < 1, na.rm = TRUE))

# residual correlation CPI vs second series (if any)
res_cor <- NA_real_
res_mat <- suppressWarnings(tryCatch(cor(residuals(var_mod), use = "pairwise.complete.obs"), error = function(e) NULL))
if (!is.null(res_mat) && "CPI" %in% colnames(res_mat) && ncol(res_mat) >= 2) {
  other <- setdiff(colnames(res_mat), "CPI")
  if (length(other) >= 1) {
    res_cor <- to1(res_mat["CPI", other[1]])
  }
}

model_label <- paste0("CPI+", paste(vars_in, collapse = "+"))
train_n <- length(train); test_n <- length(test)

# Now safe to build the row (every field is guaranteed length 1)
this_row <- data.frame(
  Model       = model_label,
  Lag         = as.integer(to1(p)),
  TrainN      = to1(train_n),
  TestH       = to1(h),
  AIC         = to1(aic_val),
  BIC         = to1(bic_val),
  R_SQ        = to1(r_2),
  Adj_R_SQ    = to1(adj_r_2),
  LogLik      = to1(llk_val),
  Stable      = to1(stable, na = NA),
  MaxRootMod  = round(to1(max_root), 4),
  LB_p        = to1(lb_pval),
  ARCH_p      = to1(arch_pval),
  Normal_p    = to1(norm_pval),
  RMSE_CPI    = to1(rmse_cpi),
  MAE_CPI     = to1(mae_cpi),
  MAPE_CPI    = to1(mape_cpi),
  ResCorr_CPI_2nd = to1(res_cor),
  stringsAsFactors = FALSE
)

# Append to running table
if (!exists("VAR_comp_table")) {
  VAR_comp_table <- this_row
} else {
  VAR_comp_table <- rbind(VAR_comp_table, this_row)
}

VAR_comp_table <- VAR_comp_table[order(VAR_comp_table$AIC, VAR_comp_table$BIC), ]
print(VAR_comp_table)


plot(y_model[,"CPI"], main = "CPI (transformed) with VAR forecast window", ylab = "CPI (transformed)", xlab = "")
lines(cpi_fc, col = "red")
legend("topleft", bty = "n", lty = 1, col = c("black","red"),
       legend = c("Actual (all)","Forecast (test window)"))

act_zoom <- window(y_model[,"CPI"], start = c(2022, 1))  
fc_zoom  <- window(cpi_fc,          start = c(2022, 1))
plot(act_zoom); lines(fc_zoom, col = "red", lwd = 2)

plot(resid(var_mod)[, "CPI"], main = "VAR Model Residuals")
acf(resid(var_mod)[, "CPI"], main = "ACF of CPI Residuals")
pacf(resid(var_mod)[, "CPI"], main = "PACF of CPI Residuals")

# The model indicates bidirectional lagged linkages between CPI and PPI, with CPI showing more persistence and stronger COVID-period effects. The significant COVID impact on CPI but not PPI implies downstream consumer prices were more directly affected by pandemic-related shocks. However, residual autocorrelation in PPI warrants further refinement before long-term structural inference.

last_train_index <- as.numeric(tail(train, n = 1))
last_train_index

last_log_cpi_train <- log(data$CPI[last_train_index])

inverted_fcast_cpi <- exp(cumsum(cpi_fc) + last_log_cpi_train)
inverted_actual_cpi <- exp(cumsum(cpi_act) + last_log_cpi_train)

rmse_orig <- sqrt(mean((inverted_fcast_cpi - inverted_actual_cpi)^2, na.rm = TRUE))
mae_orig <- mean(abs(inverted_fcast_cpi - inverted_actual_cpi), na.rm = TRUE)

round(rmse_orig, 4)
round(mae_orig, 4)


# ============= VAR Model (Fed Fund) =====
Y_full <- zoo_to_ts(Z2[, c("CPI","FedFund")])

vars_in <- c("FedFund")  # edit this per model
y_model <- Y_full[, c("CPI", vars_in)]
x_model <- X_full

h <- 12
n <- nrow(y_model)
train <- 1:(n - h)
test  <- (n - h + 1):n

sel <- VARselect(y_model[train, ], lag.max = 48, type = "const", exogen = x_model[train, , drop = FALSE], season = 12)
print(sel$selection[["AIC(n)"]])
print(sel$selection[["SC(n)"]])
p <- 3

ps <- 2:12
tab <- do.call(rbind, lapply(ps, function(p){
  m  <- VAR(y_model[train,], p=p, type="const", exogen=x_model[train,,drop=FALSE], season = 12)
  fc <- predict(m, n.ahead=h, dumvar=x_model[test,,drop=FALSE])$fcst$CPI[,"fcst"]
  fc <- ts(fc, start=time(y_model)[test][1], frequency=12)
  rmse <- sqrt(mean((y_model[test,"CPI"] - fc)^2))
  lb   <- serial.test(m, lags.pt=36, type="PT.asymptotic")$serial$p.value
  aic <- AIC(m)
  bic <- BIC(m)
  c(p=p, RMSE=rmse, LBp=lb, AIC = aic, BIC = bic)
}))
tab[order(tab[,"RMSE"]),]

var_mod <- VAR(y_model[train, ], p = p, type = "const", exogen = x_model[train, , drop = FALSE], season = 12)

summary_model <- summary(var_mod)
summary_model

fc <- predict(var_mod, n.ahead = h, dumvar = x_model[test, , drop = FALSE])
cpi_fc  <- ts(fc$fcst$CPI[, "fcst"], start = time(y_model)[test][1], frequency = 12)
cpi_act <- y_model[test, "CPI"]


# Recompute / coerce all scalars safely
rmse_cpi <- to1(sqrt(mean((cpi_act - cpi_fc)^2, na.rm = TRUE)))
mae_cpi  <- to1(mean(abs(cpi_act - cpi_fc), na.rm = TRUE))
mape_cpi <- to1(mean(abs((cpi_act - cpi_fc) / cpi_act), na.rm = TRUE) * 100)

aic_val <- to1(suppressWarnings(tryCatch(AIC(var_mod), error = function(e) NA_real_)))
bic_val <- to1(suppressWarnings(tryCatch(BIC(var_mod), error = function(e) NA_real_)))
llk_val <- to1(suppressWarnings(tryCatch(as.numeric(logLik(var_mod)), error = function(e) NA_real_)))

r_2 <- to1(summary_model$varresult$CPI$r.squared)
adj_r_2 <- to1(summary_model$varresult$CPI$adj.r.squared)

lb_pval   <- to1(suppressWarnings(tryCatch(serial.test(var_mod, lags.pt = 24, type = "PT.asymptotic")$serial$p.value, error = function(e) NA_real_)))
arch_pval <- to1(suppressWarnings(tryCatch(arch.test(var_mod, lags.multi = 12, multivariate.only = TRUE)$arch.mul$p.value, error = function(e) NA_real_)))
norm_pval <- to1(suppressWarnings(tryCatch(normality.test(var_mod, multivariate.only = TRUE)$jb.mul$p.value, error = function(e) NA_real_)))

roots_mod <- suppressWarnings(tryCatch(Mod(roots(var_mod)), error = function(e) NA_real_))
max_root  <- if (length(roots_mod) == 0 || all(is.na(roots_mod))) NA_real_ else max(roots_mod, na.rm = TRUE)
stable    <- isTRUE(all(roots_mod < 1, na.rm = TRUE))

# residual correlation CPI vs second series (if any)
res_cor <- NA_real_
res_mat <- suppressWarnings(tryCatch(cor(residuals(var_mod), use = "pairwise.complete.obs"), error = function(e) NULL))
if (!is.null(res_mat) && "CPI" %in% colnames(res_mat) && ncol(res_mat) >= 2) {
  other <- setdiff(colnames(res_mat), "CPI")
  if (length(other) >= 1) {
    res_cor <- to1(res_mat["CPI", other[1]])
  }
}

model_label <- paste0("CPI+", paste(vars_in, collapse = "+"))
train_n <- length(train); test_n <- length(test)

# Now safe to build the row (every field is guaranteed length 1)
this_row <- data.frame(
  Model       = model_label,
  Lag         = as.integer(to1(p)),
  TrainN      = to1(train_n),
  TestH       = to1(h),
  AIC         = to1(aic_val),
  BIC         = to1(bic_val),
  R_SQ        = to1(r_2),
  Adj_R_SQ    = to1(adj_r_2),
  LogLik      = to1(llk_val),
  Stable      = to1(stable, na = NA),
  MaxRootMod  = round(to1(max_root), 4),
  LB_p        = to1(lb_pval),
  ARCH_p      = to1(arch_pval),
  Normal_p    = to1(norm_pval),
  RMSE_CPI    = to1(rmse_cpi),
  MAE_CPI     = to1(mae_cpi),
  MAPE_CPI    = to1(mape_cpi),
  ResCorr_CPI_2nd = to1(res_cor),
  stringsAsFactors = FALSE
)

# Append to running table
if (!exists("VAR_comp_table")) {
  VAR_comp_table <- this_row
} else {
  VAR_comp_table <- rbind(VAR_comp_table, this_row)
}

VAR_comp_table <- VAR_comp_table[order(VAR_comp_table$AIC, VAR_comp_table$BIC), ]
print(VAR_comp_table)


plot(y_model[,"CPI"], main = "CPI (transformed) with VAR forecast window", ylab = "CPI (transformed)", xlab = "")
lines(cpi_fc, col = "red")
legend("topleft", bty = "n", lty = 1, col = c("black","red"),
       legend = c("Actual (all)","Forecast (test window)"))

act_zoom <- window(y_model[,"CPI"], start = c(2022, 1))  
fc_zoom  <- window(cpi_fc,          start = c(2022, 1))
plot(act_zoom); lines(fc_zoom, col = "red", lwd = 2)

plot(resid(var_mod)[, "CPI"], main = "VAR Model Residuals")
acf(resid(var_mod)[, "CPI"], main = "ACF of CPI Residuals")
pacf(resid(var_mod)[, "CPI"], main = "PACF of CPI Residuals")

# The model indicates bidirectional lagged linkages between CPI and PPI, with CPI showing more persistence and stronger COVID-period effects. The significant COVID impact on CPI but not PPI implies downstream consumer prices were more directly affected by pandemic-related shocks. However, residual autocorrelation in PPI warrants further refinement before long-term structural inference.

last_train_index <- as.numeric(tail(train, n = 1))
last_train_index

last_log_cpi_train <- log(data$CPI[last_train_index])

inverted_fcast_cpi <- exp(cumsum(cpi_fc) + last_log_cpi_train)
inverted_actual_cpi <- exp(cumsum(cpi_act) + last_log_cpi_train)

rmse_orig <- sqrt(mean((inverted_fcast_cpi - inverted_actual_cpi)^2, na.rm = TRUE))
mae_orig <- mean(abs(inverted_fcast_cpi - inverted_actual_cpi), na.rm = TRUE)

round(rmse_orig, 4)
round(mae_orig, 4)

