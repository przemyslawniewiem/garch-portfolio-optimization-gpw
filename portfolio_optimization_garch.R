getwd()

# PROJEKT EFID - CZĘŚĆ A

# Instalacja i ładowanie bibliotek
if (!require("quantmod")) install.packages("quantmod")
if (!require("PerformanceAnalytics")) install.packages("PerformanceAnalytics")
if (!require("quadprog")) install.packages("quadprog")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("gridExtra")) install.packages("gridExtra")

library(quantmod)
library(PerformanceAnalytics)
library(quadprog)
library(ggplot2)
library(gridExtra)

# Wczytywanie danych
file_names <- c("PZU.csv", "BNP.csv", "CDR.csv", "ORANGE.csv", "ORLEN.csv")
ticker_names <- c("PZU", "BNP", "CDR", "ORANGE", "ORLEN")

data_list <- list()
for (i in 1:length(file_names)) {
  if (!file.exists(file_names[i])) {
    stop(paste("Brak pliku:", file_names[i]))
  }
  df <- read.csv(file_names[i])
  xts_obj <- xts(df$Zamkniecie, order.by = as.Date(df$Data))
  colnames(xts_obj) <- ticker_names[i]
  data_list[[i]] <- xts_obj
}

all_prices <- na.omit(do.call(merge, data_list))
returns <- na.omit(Return.calculate(all_prices, method = "log"))

# Wstępna analiza danych
cat("\n>>> STATYSTYKI CEN ZAMKNIĘCIA:\n")
price_stats <- data.frame(
  Instrument = ticker_names,
  Min = apply(all_prices, 2, min),
  Średnia = apply(all_prices, 2, mean),
  Mediana = apply(all_prices, 2, median),
  Max = apply(all_prices, 2, max),
  Odch_std = apply(all_prices, 2, sd)
)
print(price_stats[, -1], digits = 2)

cat("\n>>> STATYSTYKI STÓP ZWROTU:\n")
return_stats <- data.frame(
  Instrument = ticker_names,
  Średnia_proc = apply(returns, 2, mean, na.rm = TRUE) * 100,
  Mediana_proc = apply(returns, 2, median, na.rm = TRUE) * 100,
  Odch_std_proc = apply(returns, 2, sd, na.rm = TRUE) * 100,
  Min_proc = apply(returns, 2, min, na.rm = TRUE) * 100,
  Max_proc = apply(returns, 2, max, na.rm = TRUE) * 100,
  Skośność = apply(returns, 2, skewness, na.rm = TRUE),
  Kurtoza = apply(returns, 2, kurtosis, na.rm = TRUE)
)
print(return_stats[, -1], digits = 4)

cat("\n>>> MACIERZ KORELACJI:\n")
cor_matrix <- cor(returns)
print(round(cor_matrix, 3))

# Wykresy danych wejściowych
plot.zoo(all_prices, 
         main = "Ceny zamknięcia instrumentów",
         xlab = "Data", 
         ylab = "Cena",
         col = 1:5,
         lwd = 1.5,
         plot.type = "single")
legend("topleft", legend = ticker_names, col = 1:5, lty = 1, lwd = 2, cex = 0.8, ncol = 2)
grid()

par(mfrow = c(3, 2), mar = c(3, 4, 2, 1))
for (i in 1:ncol(returns)) {
  plot(returns[, i], main = ticker_names[i], ylab = "Stopa zwrotu", col = "steelblue", cex.main = 0.9)
  abline(h = 0, col = "red", lty = 2)
  grid()
}
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2))

# Parametry okna kroczącego
window_len <- 500
step_len <- 5
n_obs <- nrow(returns)

# Funkcja wag MVP
get_mvp_weights <- function(sigma) {
  n <- ncol(sigma)
  Dmat <- 2 * sigma
  dvec <- rep(0, n)
  Amat <- cbind(rep(1, n), diag(n))
  bvec <- c(1, rep(0, n))
  sol <- solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat, bvec = bvec, meq = 1)
  weights <- sol$solution
  names(weights) <- colnames(sigma)
  return(weights)
}

# Konstrukcja portfela MVP
mvp_returns_list <- list()
weight_history <- list()
iteration_count <- 0

for (i in seq(window_len, n_obs - step_len, by = step_len)) {
  iteration_count <- iteration_count + 1
  train_rets <- returns[(i - window_len + 1):i, ]
  test_end <- min(i + step_len, n_obs)
  test_rets <- returns[(i + 1):test_end, ]
  sigma <- cov(train_rets)
  w_mvp <- get_mvp_weights(sigma)
  weight_history[[iteration_count]] <- data.frame(Date = index(test_rets)[1], t(w_mvp))
  port_ret <- Return.portfolio(test_rets, weights = w_mvp)
  mvp_returns_list[[length(mvp_returns_list) + 1]] <- port_ret
}

mvp_returns <- do.call(rbind, mvp_returns_list)
colnames(mvp_returns) <- "MVP"

weights_df <- do.call(rbind, weight_history)
avg_weights <- colMeans(weights_df[, -1])

cat("\n>>> ŚREDNIE WAGI PORTFELA MVP:\n")
weight_summary <- data.frame(
  Instrument = ticker_names,
  Średnia_waga = avg_weights,
  Min_waga = apply(weights_df[, -1], 2, min),
  Max_waga = apply(weights_df[, -1], 2, max)
)
print(weight_summary[, -1], digits = 4)

# Przygotowanie danych do porównania (2019-2025)
period <- "2019-01-01/2025-12-31"
mvp_final <- mvp_returns[period]
common_dates <- index(mvp_final)

ew_weights <- rep(1 / ncol(returns), ncol(returns))
ew_returns <- Return.portfolio(returns[common_dates], weights = ew_weights)
colnames(ew_returns) <- "EW"

single_returns <- returns[common_dates]
colnames(single_returns) <- paste0("SINGLE_", ticker_names)

comparison_matrix <- merge(mvp_final, ew_returns, single_returns)

# Funkcja statystyk
calc_stats <- function(ret_series) {
  ret_series <- na.omit(ret_series)
  mean_daily <- mean(ret_series, na.rm = TRUE)
  mean_annual <- mean_daily * 252
  sd_daily <- sd(ret_series, na.rm = TRUE)
  sd_annual <- sd_daily * sqrt(252)
  var_95 <- as.numeric(VaR(ret_series, p = 0.95, method = "historical"))
  es_95 <- as.numeric(ES(ret_series, p = 0.95, method = "historical"))
  sharpe <- mean_annual / sd_annual
  
  res <- c(
    "Śr. stopa zwrotu (dzienna) %" = mean_daily * 100,
    "Śr. stopa zwrotu (roczna) %" = mean_annual * 100,
    "Odch. std (dzienne) %" = sd_daily * 100,
    "Odch. std (roczne) %" = sd_annual * 100,
    "VaR 95% (hist)" = var_95,
    "ES 95% (hist)" = es_95,
    "Sharpe Ratio (roczny)" = sharpe
  )
  return(res)
}

# Podsumowanie dla całego okresu (2019-2025)
cat("\n>>> STATYSTYKI - CAŁY OKRES (2019-2025):\n")
total_summary <- t(apply(comparison_matrix, 2, calc_stats))
print(round(total_summary, 4))

cat("\n>>> RANKING STRATEGII (wg Sharpe Ratio):\n")
ranking <- data.frame(
  Strategia = rownames(total_summary),
  Sharpe = total_summary[, "Sharpe Ratio (roczny)"],
  Srednia_roczna = total_summary[, "Śr. stopa zwrotu (roczna) %"],
  Ryzyko_roczne = total_summary[, "Odch. std (roczne) %"]
)
ranking <- ranking[order(-ranking$Sharpe), ]
rownames(ranking) <- NULL
print(ranking[, -1], digits = 4)

# Statystyki w poszczególnych latach
years <- unique(format(index(comparison_matrix), "%Y"))
yearly_results <- list()

for (yr in years) {
  cat(paste0("\n>>> STATYSTYKI - ROK ", yr, ":\n"))
  yearly_data <- comparison_matrix[yr]
  yearly_summary <- t(apply(yearly_data, 2, calc_stats))
  print(round(yearly_summary, 4))
  yearly_results[[yr]] <- yearly_summary
}

# Wizualizacje
par(mfrow = c(2, 1), mar = c(4, 4, 3, 2))

# Cumulative returns
cum_mvp_ew <- cumprod(1 + comparison_matrix[, c("MVP", "EW")])
plot(index(cum_mvp_ew), cum_mvp_ew[, "MVP"],
     type = "l", lwd = 2.5, col = "darkblue",
     ylim = range(cum_mvp_ew),
     main = "Cumulative Returns",
     xlab = "", ylab = "Wartość portfela",
     xaxt = "n")
lines(index(cum_mvp_ew), cum_mvp_ew[, "EW"], lwd = 2.5, col = "darkred")
legend("topleft", legend = c("MVP", "EW"), col = c("darkblue", "darkred"), lty = 1, lwd = 2.5, cex = 0.8)
grid()

# Drawdown
dd_mvp <- PerformanceAnalytics:::Drawdowns(comparison_matrix[, "MVP"])
dd_ew <- PerformanceAnalytics:::Drawdowns(comparison_matrix[, "EW"])
plot(index(dd_mvp), dd_mvp,
     type = "l", lwd = 2, col = "darkblue",
     ylim = range(c(dd_mvp, dd_ew)),
     main = "Drawdown",
     xlab = "Data", ylab = "Drawdown")
lines(index(dd_ew), dd_ew, lwd = 2, col = "darkred")
abline(h = 0, col = "gray", lty = 2)
grid()

par(mfrow = c(1, 1), mar = c(5, 4, 4, 2))

# Krzywe kapitału - wszystkie strategie
cum_returns <- cumprod(1 + comparison_matrix)

plot(index(cum_returns), cum_returns[, "MVP"], 
     type = "l", lwd = 2.5, col = "darkblue",
     ylim = range(cum_returns),
     main = "Krzywe kapitału wszystkich strategii",
     xlab = "Data", ylab = "Wartość portfela")

lines(index(cum_returns), cum_returns[, "EW"], lwd = 2.5, col = "darkred")

for (i in 1:length(ticker_names)) {
  col_name <- paste0("SINGLE_", ticker_names[i])
  lines(index(cum_returns), cum_returns[, col_name], lwd = 1, lty = 2, col = i + 2)
}

legend("topleft", 
       legend = c("MVP", "EW", paste0("SINGLE_", ticker_names)),
       col = c("darkblue", "darkred", (1:length(ticker_names)) + 2),
       lty = c(1, 1, rep(2, length(ticker_names))),
       lwd = c(2.5, 2.5, rep(1, length(ticker_names))),
       cex = 0.7, ncol = 2)
grid()

# Sharpe Ratio
par(mfrow = c(1, 1), mar = c(8, 4, 4, 2))
sharpe_vals <- total_summary[, "Sharpe Ratio (roczny)"]
barplot(sharpe_vals,
        main = "Sharpe Ratio - porównanie strategii",
        ylab = "Sharpe Ratio",
        las = 2,
        col = ifelse(sharpe_vals > 0, "forestgreen", "firebrick"),
        ylim = c(min(sharpe_vals) * 1.2, max(sharpe_vals) * 1.2))
abline(h = 0, col = "black", lwd = 2)
grid()

# VaR i ES
par(mfrow = c(1, 2), mar = c(8, 4, 4, 2))
var_vals <- total_summary[, "VaR 95% (hist)"] * 100
barplot(abs(var_vals),
        main = "VaR 95%",
        ylab = "VaR (%)",
        las = 2,
        col = "coral")
grid()

es_vals <- total_summary[, "ES 95% (hist)"] * 100
barplot(abs(es_vals),
        main = "Expected Shortfall 95%",
        ylab = "ES (%)",
        las = 2,
        col = "darkred")
grid()

# Ewolucja wag MVP
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2))
weights_ts <- xts(weights_df[, -1], order.by = weights_df$Date)

plot.zoo(weights_ts,
         main = "Ewolucja wag portfela MVP",
         xlab = "Data",
         ylab = "Waga",
         col = 1:5,
         lwd = 1.5,
         plot.type = "single",
         ylim = c(0, 1))
legend("topright", legend = ticker_names, col = 1:5, lty = 1, lwd = 2, cex = 0.8)
grid()

# Stacked area chart wag
weights_mat <- as.matrix(weights_df[, -1])
dates_weights <- weights_df$Date
n_cols <- ncol(weights_mat)

# Definicja kolorów 
cols <- c("black", "red", "green", "blue", "cyan") 

# Obliczenie skumulowanych wag
cum_weights <- t(apply(weights_mat, 1, cumsum))

plot(dates_weights, cum_weights[, n_cols],
     type = "n",
     ylim = c(0, 1),
     main = "Struktura portfela MVP (stacked)",
     xlab = "Data",
     ylab = "Udział",
     xaxs = "i", yaxs = "i")

for (i in n_cols:1) {
  if (i == 1) {
    polygon(c(dates_weights, rev(dates_weights)),
            c(cum_weights[, i], rep(0, nrow(weights_mat))),
            col = cols[i], border = "white", lwd = 0.5)
  } else {
    polygon(c(dates_weights, rev(dates_weights)),
            c(cum_weights[, i], rev(cum_weights[, i-1])),
            col = cols[i], border = "white", lwd = 0.5)
  }
}

legend("topright", 
       legend = rev(ticker_names), 
       fill = rev(cols), 
       cex = 0.8, 
       bg = "white")


# ================================================================================
# PROJEKT EFID - CZĘŚĆ B (MODELE GARCH)
# ================================================================================

cat("\n\n")
cat("================================================================================\n")
cat("                         CZĘŚĆ B - MODELE GARCH                                 \n")
cat("================================================================================\n")

# ================================================================================
# 1. BIBLIOTEKI
# ================================================================================
if (!require("rugarch")) install.packages("rugarch")
if (!require("rmgarch")) install.packages("rmgarch")
if (!require("mgarchBEKK")) install.packages("mgarchBEKK")
if (!require("vars")) install.packages("vars")
if (!require("tseries")) install.packages("tseries")
if (!require("parallel")) install.packages("parallel")
if (!require("doParallel")) install.packages("doParallel")
if (!require("foreach")) install.packages("foreach")
if (!require("doSNOW")) install.packages("doSNOW")
if (!require("forecast")) install.packages("forecast")

library(forecast)
library(rugarch)
library(rmgarch)
library(mgarchBEKK)
library(vars)
library(tseries)
library(parallel)
library(doParallel)
library(foreach)
library(doSNOW)

# ================================================================================
# 2. PRZYGOTOWANIE DANYCH
# ================================================================================

cat("\n>>> Przygotowanie danych do analizy GARCH...\n")

# Konwersja danych do formatu macierzowego
ret_cols <- ticker_names
mat_rets <- as.matrix(returns)
colnames(mat_rets) <- ret_cols
vec_dates <- index(returns)
n_assets <- ncol(mat_rets)

# Parametry analizy
window_length <- 500
hold_period   <- 5
bekk_step     <- 50
start_target  <- as.Date("2019-01-01")

start_idx <- which(vec_dates >= start_target)[1]

# Walidacja
if (is.na(start_idx)) stop("Brak danych po 2019-01-01")
if (start_idx <= window_length) stop("Za mało danych historycznych przed 2019")

cat("    - Długość okna estymacji:", window_length, "dni\n")
cat("    - Częstotliwość rebalansowania:", hold_period, "dni\n")
cat("    - Częstotliwość estymacji BEKK:", bekk_step, "dni\n")
cat("    - Data startu:", as.character(start_target), "\n")

# Funkcja MVP
calc_mvp_weights <- function(Sigma) {
  inv_Sigma <- tryCatch(solve(Sigma), error = function(e) NULL)
  if (!is.null(inv_Sigma)) {
    ones <- rep(1, ncol(Sigma))
    numerator <- inv_Sigma %*% ones
    denominator <- as.numeric(t(ones) %*% inv_Sigma %*% ones)
    return(as.numeric(numerator / denominator))
  } else {
    return(rep(1/ncol(Sigma), ncol(Sigma)))
  }
}

# Funkcja statystyk
calc_stats_partB <- function(ret_series) {
  ret_series <- na.omit(ret_series)
  mean_daily <- mean(ret_series, na.rm = TRUE)
  mean_annual <- mean_daily * 252
  sd_daily <- sd(ret_series, na.rm = TRUE)
  sd_annual <- sd_daily * sqrt(252)
  
  var_95 <- as.numeric(PerformanceAnalytics::VaR(ret_series, p = 0.95, method = "historical"))
  es_95 <- as.numeric(PerformanceAnalytics::ES(ret_series, p = 0.95, method = "historical"))
  sharpe <- mean_annual / sd_annual
  
  res <- c(
    "Śr. dzienna %" = mean_daily * 100,
    "Śr. roczna %" = mean_annual * 100,
    "SD dzienna %" = sd_daily * 100,
    "SD roczna %" = sd_annual * 100,
    "VaR 95%" = var_95,
    "ES 95%" = es_95,
    "Sharpe" = sharpe
  )
  return(res)
}

# ================================================================================
# 3. MVP STATYCZNY (REKALKULACJA DLA SPÓJNOŚCI)
# ================================================================================

cat("\n>>> Obliczenia MVP Static...\n")

portfolio_mvp_static <- list()
current_idx <- start_idx

while ((current_idx + hold_period - 1) <= nrow(mat_rets)) {
  
  train_data  <- mat_rets[(current_idx - window_length):(current_idx - 1), ]
  Sigma <- cov(train_data) 
  w_mvp <- calc_mvp_weights(Sigma)
  
  test_end   <- min(current_idx + hold_period - 1, nrow(mat_rets))
  test_data  <- mat_rets[current_idx:test_end, , drop=FALSE]
  
  port_ret   <- test_data %*% w_mvp
  
  portfolio_mvp_static[[length(portfolio_mvp_static) + 1]] <- data.frame(
    Date = vec_dates[current_idx:test_end],
    MVP_Static = as.numeric(port_ret)
  )
  
  current_idx <- current_idx + hold_period
}

mvp_static_results <- do.call(rbind, portfolio_mvp_static)

cat("    ✓ MVP Static - gotowe\n")

# ================================================================================
# 4. DCC-GARCH (NORMALNY + T-STUDENTA)
# ================================================================================

cat("\n>>> Obliczenia DCC-GARCH...\n")

# Specyfikacje GARCH
uspec_n <- ugarchspec(
  mean.model = list(armaOrder = c(0,0), include.mean = FALSE),
  variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
  distribution.model = "norm"
)

uspec_t <- ugarchspec(
  mean.model = list(armaOrder = c(0,0), include.mean = FALSE),
  variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
  distribution.model = "std"
)

# Specyfikacje DCC
spec_dcc_n <- dccspec(
  uspec = multispec(replicate(n_assets, uspec_n)), 
  dccOrder = c(1,1), 
  distribution = "mvnorm"
)

spec_dcc_t <- dccspec(
  uspec = multispec(replicate(n_assets, uspec_t)), 
  dccOrder = c(1,1), 
  distribution = "mvt"
)

# Konfiguracja obliczeń równoległych
cl <- makeCluster(parallel::detectCores() - 1)
registerDoParallel(cl)

indices_dcc <- seq(start_idx, nrow(mat_rets) - hold_period + 1, by = hold_period)

cat("    - Liczba iteracji DCC:", length(indices_dcc), "\n")

results_dcc_list <- foreach(
  current_idx = indices_dcc, 
  .packages = c("rugarch", "rmgarch", "vars"), 
  .export = c("calc_mvp_weights", "spec_dcc_n", "spec_dcc_t", "mat_rets", 
              "window_length", "hold_period", "vec_dates", "n_assets")
) %dopar% {
  
  train_data  <- mat_rets[(current_idx - window_length):(current_idx - 1), ]
  
  # VAR -> Reszty
  var_fit <- tryCatch(
    VAR(train_data, p=1, type="const"), 
    error=function(e) NULL
  )
  residuals_var <- if(!is.null(var_fit)) residuals(var_fit) else train_data
  
  # Fallback
  static_cov <- cov(train_data)
  
  # DCC Normalny
  H_n <- tryCatch({
    fit <- dccfit(spec_dcc_n, data = residuals_var, solver = "solnp")
    if(inherits(fit, "DCCfit")) {
      rcov(dccforecast(fit, n.ahead = 1))[[1]][,,1]
    } else {
      static_cov
    }
  }, error = function(e) static_cov)
  
  # DCC t-Student
  H_t <- tryCatch({
    fit <- dccfit(spec_dcc_t, data = residuals_var, solver = "solnp")
    if(inherits(fit, "DCCfit")) {
      rcov(dccforecast(fit, n.ahead = 1))[[1]][,,1]
    } else {
      static_cov
    }
  }, error = function(e) static_cov)
  
  # Wagi MVP
  w_n <- calc_mvp_weights(H_n)
  w_t <- calc_mvp_weights(H_t)
  
  # Wyniki na okresie testowym
  test_end   <- min(current_idx + hold_period - 1, nrow(mat_rets))
  test_data  <- mat_rets[current_idx:test_end, , drop=FALSE]
  
  data.frame(
    Date = vec_dates[current_idx:test_end],
    DCC_Normal = as.numeric(test_data %*% w_n),
    DCC_Student = as.numeric(test_data %*% w_t)
  )
}

stopCluster(cl)

results_dcc_df <- do.call(rbind, results_dcc_list)

cat("    ✓ DCC-GARCH - gotowe\n")

# ================================================================================
# 5. BEKK-GARCH
# ================================================================================

cat("\n>>> Obliczenia BEKK-GARCH...\n")
cat("    [To może potrwać 1-2 godziny...]\n")

# Punkty estymacji (co 50 dni)
indices_bekk <- seq(start_idx, nrow(mat_rets) - bekk_step + 1, by = bekk_step)
iterations <- length(indices_bekk)

cat("    - Liczba okien estymacji BEKK:", iterations, "\n")

# Konfiguracja klastra
cores <- parallel::detectCores() - 1
cl <- makeCluster(cores)
registerDoSNOW(cl)

# Pasek postępu
pb <- txtProgressBar(max = iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

# Główna pętla równoległa BEKK
results_bekk_list <- foreach(
  i = 1:iterations, 
  .packages = c("mgarchBEKK", "vars"), 
  .export = c("calc_mvp_weights", "mat_rets", "vec_dates", 
              "indices_bekk", "window_length", "bekk_step", "hold_period"),
  .options.snow = opts
) %dopar% {
  
  current_idx <- indices_bekk[i]
  
  # Dane treningowe
  train_data  <- mat_rets[(current_idx - window_length):(current_idx - 1), ]
  
  # VAR -> Reszty
  var_fit <- tryCatch(
    VAR(train_data, p=1, type="const"), 
    error=function(e) NULL
  )
  res_var <- if(!is.null(var_fit)) residuals(var_fit) else train_data
  
  # Estymacja BEKK i prognoza na 50 dni
  H_forecasts <- tryCatch({
    fit_bekk <- BEKK(res_var, order = c(1, 1), method = "BFGS")
    predict(fit_bekk, n.ahead = bekk_step)$H
  }, error = function(e) NULL)
  
  # Fallback w razie błędu
  if(is.null(H_forecasts)) {
    static <- cov(train_data)
    H_forecasts <- replicate(bekk_step, static, simplify = FALSE)
  }
  
  # Rebalansowanie wewnątrz 50-dniowego okna prognozy
  days_available <- min(current_idx + bekk_step - 1, nrow(mat_rets)) - current_idx + 1
  results_chunk <- list()
  
  if (days_available > 0) {
    sub_indices <- seq(1, days_available, by = hold_period)
    
    for (k in sub_indices) {
      # Prognozowana macierz na dany dzień
      H_curr <- H_forecasts[[k]]
      
      # Wagi MVP
      w_bekk <- calc_mvp_weights(H_curr)
      
      # Zastosowanie wag na kolejne 5 dni
      sub_end <- min(k + hold_period - 1, days_available)
      real_idx_start <- current_idx + k - 1
      real_idx_end   <- current_idx + sub_end - 1
      
      chunk_data <- mat_rets[real_idx_start:real_idx_end, , drop=FALSE]
      chunk_dates <- vec_dates[real_idx_start:real_idx_end]
      
      results_chunk[[length(results_chunk) + 1]] <- data.frame(
        Date = chunk_dates,
        BEKK = as.numeric(chunk_data %*% w_bekk)
      )
    }
  }
  
  if(length(results_chunk) > 0) do.call(rbind, results_chunk) else NULL
}

close(pb)
stopCluster(cl)

results_bekk_df <- do.call(rbind, results_bekk_list)

cat("\n    ✓ BEKK-GARCH - gotowe\n")

# ================================================================================
# 6. AGREGACJA WYNIKÓW
# ================================================================================

cat("\n>>> Łączenie wyników...\n")

# Łączenie wszystkich strategii (BEZ PIPE OPERATOR)
final_df <- merge(mvp_static_results, results_dcc_df, by = "Date", all = TRUE)
final_df <- merge(final_df, results_bekk_df, by = "Date", all = TRUE)

# Dodanie EW i pojedynczych aktywów
assets_df <- data.frame(
  Date = vec_dates,
  EW = rowMeans(mat_rets),
  mat_rets
)

final_df <- merge(final_df, assets_df, by = "Date", all = TRUE)

# Usunięcie braków
final_df <- na.omit(final_df)

# Przygotowanie nazw kolumn dla pojedynczych aktywów
single_cols <- paste0("SINGLE_", ticker_names)
colnames(final_df)[(ncol(final_df) - n_assets + 1):ncol(final_df)] <- single_cols

# Konwersja do xts
all_strategies_B <- xts(
  final_df[, -1], 
  order.by = final_df$Date
)

mvp_models <- c("MVP_Static", "DCC_Normal", "DCC_Student", "BEKK")

cat("    ✓ Strategie gotowe:", paste(colnames(all_strategies_B), collapse = ", "), "\n")

# ================================================================================
# 7. STATYSTYKI - CAŁY OKRES
# ================================================================================

cat("\n================================================================================\n")
cat("                    STATYSTYKI - CAŁY OKRES (2019-2025)                        \n")
cat("================================================================================\n")

stats_overall_B <- t(apply(all_strategies_B, 2, calc_stats_partB))
print(round(stats_overall_B, 4))

cat("\n>>> RANKING STRATEGII (według Sharpe Ratio):\n")
ranking_B <- data.frame(
  Strategia = rownames(stats_overall_B),
  Sharpe = stats_overall_B[, "Sharpe"],
  Srednia_roczna = stats_overall_B[, "Śr. roczna %"],
  Ryzyko_roczne = stats_overall_B[, "SD roczna %"],
  VaR_95 = stats_overall_B[, "VaR 95%"],
  ES_95 = stats_overall_B[, "ES 95%"]
)
ranking_B <- ranking_B[order(-ranking_B$Sharpe), ]
rownames(ranking_B) <- NULL
print(ranking_B, digits = 4)

# ================================================================================
# 8. STATYSTYKI ROCZNE
# ================================================================================

cat("\n================================================================================\n")
cat("                          STATYSTYKI ROCZNE                                     \n")
cat("================================================================================\n")

years_B <- unique(format(index(all_strategies_B), "%Y"))
yearly_results_B <- list()

for (yr in years_B) {
  cat(paste0("\n>>> ROK ", yr, ":\n"))
  yearly_data_B <- all_strategies_B[yr]
  yearly_summary_B <- t(apply(yearly_data_B, 2, calc_stats_partB))
  print(round(yearly_summary_B, 4))
  yearly_results_B[[yr]] <- yearly_summary_B
}

# ================================================================================
# 9. TESTY STATYSTYCZNE
# ================================================================================

cat("\n================================================================================\n")
cat("                          TESTY STATYSTYCZNE                                    \n")
cat("================================================================================\n")

cat("\n>>> Test Diebold-Mariano (DM) - porównanie z MVP_Static:\n")
cat("    (H0: modele równie dobre; H1: model 2 lepszy)\n\n")

for (model in mvp_models[-1]) {
  e1 <- (all_strategies_B[, "MVP_Static"])^2
  e2 <- (all_strategies_B[, model])^2
  
  dm_test <- tryCatch({
    dm.test(e1, e2, alternative = "greater", h = 1)
  }, error = function(e) {
    return(NULL)
  })
  
  if (!is.null(dm_test)) {
    cat(sprintf("  %s vs MVP_Static:\n", model))
    cat(sprintf("    DM statystyka = %.4f, p-value = %.4f %s\n", 
                dm_test$statistic, 
                dm_test$p.value,
                ifelse(dm_test$p.value < 0.05, "***", 
                       ifelse(dm_test$p.value < 0.10, "*", ""))))
  }
}

cat("\n>>> Test Jarque-Bera (normalność rozkładu zwrotów):\n")
cat("    (H0: rozkład normalny; H1: rozkład nienormalny)\n\n")

for (model in mvp_models) {
  jb_test <- tryCatch({
    jarque.bera.test(as.numeric(all_strategies_B[, model]))
  }, error = function(e) {
    return(NULL)
  })
  
  if (!is.null(jb_test)) {
    cat(sprintf("  %s:\n", model))
    cat(sprintf("    JB statystyka = %.4f, p-value = %.4f %s\n", 
                jb_test$statistic, 
                jb_test$p.value,
                ifelse(jb_test$p.value < 0.05, "***", "")))
  }
}

# ================================================================================
# 10. WIZUALIZACJE
# ================================================================================

cat("\n>>> Generowanie wizualizacji...\n")

n_strat <- ncol(all_strategies_B)
colors_B <- rep("grey70", n_strat)
lwd_B <- rep(1, n_strat)
lty_B <- rep(2, n_strat)

# Kolorowanie modeli MVP
colors_B[colnames(all_strategies_B) == "MVP_Static"] <- "black"
colors_B[colnames(all_strategies_B) == "DCC_Normal"] <- "blue"
colors_B[colnames(all_strategies_B) == "DCC_Student"] <- "darkgreen"
colors_B[colnames(all_strategies_B) == "BEKK"] <- "purple"

lwd_B[colnames(all_strategies_B) %in% mvp_models] <- 3
lty_B[colnames(all_strategies_B) %in% mvp_models] <- 1

# EW
colors_B[colnames(all_strategies_B) == "EW"] <- "orange"
lwd_B[colnames(all_strategies_B) == "EW"] <- 2
lty_B[colnames(all_strategies_B) == "EW"] <- 1

# Wykres 1: Performance Summary
par(mfrow = c(2, 1), mar = c(4, 4, 3, 2))

# Cumulative Returns
cum_all <- cumprod(1 + all_strategies_B)
plot(index(cum_all), cum_all[, "MVP_Static"],
     type = "n",
     ylim = range(cum_all),
     main = "Cumulative Returns - Wszystkie strategie",
     xlab = "",
     ylab = "Wartość portfela")

for (i in 1:ncol(all_strategies_B)) {
  lines(index(cum_all), cum_all[, i],
        col = colors_B[i],
        lwd = lwd_B[i],
        lty = lty_B[i])
}

legend("topleft",
       legend = c("MVP Static", "DCC-Normal", "DCC-Student", "BEKK", "EW"),
       col = c("black", "blue", "darkgreen", "purple", "orange"),
       lty = 1,
       lwd = c(3, 3, 3, 3, 2),
       cex = 0.6,
       ncol = 2)
grid()

# Drawdown
dd_all <- PerformanceAnalytics:::Drawdowns(all_strategies_B[, "MVP_Static"])
plot(index(dd_all), dd_all,
     type = "n",
     ylim = range(sapply(1:length(mvp_models), function(i) 
       range(PerformanceAnalytics:::Drawdowns(all_strategies_B[, mvp_models[i]])))),
     main = "Drawdown - Modele MVP",
     xlab = "Data",
     ylab = "Drawdown")

for (i in 1:length(mvp_models)) {
  dd <- PerformanceAnalytics:::Drawdowns(all_strategies_B[, mvp_models[i]])
  lines(index(dd), dd,
        col = c("black", "blue", "darkgreen", "purple")[i],
        lwd = 2)
}

abline(h = 0, col = "gray", lty = 2)
grid()

par(mfrow = c(1, 1), mar = c(5, 4, 4, 2))

# Wykres 2: Kroczące ryzyko
safe_width <- 60

if (nrow(all_strategies_B) > safe_width) {
  chart.RollingPerformance(
    all_strategies_B,
    width = safe_width,
    FUN = "StdDev.annualized",
    scale = 252,
    main = paste0("Kroczące ryzyko (", safe_width, " dni)"),
    colorset = colors_B,
    lwd = lwd_B,
    lty = lty_B,
    legend.loc = "topleft",
    cex.legend = 0.7
  )
}

# Wykres 3: Sharpe Ratio
sharpe_mvp <- stats_overall_B[mvp_models, "Sharpe"]

par(mar = c(8, 4, 4, 2))
bp <- barplot(sharpe_mvp,
              main = "Sharpe Ratio - porównanie modeli MVP",
              ylab = "Sharpe Ratio",
              las = 2,
              col = c("black", "blue", "darkgreen", "purple"),
              ylim = c(0, max(sharpe_mvp) * 1.2),
              names.arg = c("Static", "DCC-Normal", "DCC-Student", "BEKK"))

text(bp, sharpe_mvp, labels = round(sharpe_mvp, 3), 
     pos = 3, cex = 0.9, font = 2)

abline(h = 0, col = "red", lwd = 2, lty = 2)
grid()
par(mar = c(5, 4, 4, 2))

# Wykres 4: VaR i ES
par(mfrow = c(1, 2), mar = c(8, 4, 4, 2))

var_mvp <- abs(stats_overall_B[mvp_models, "VaR 95%"] * 100)
bp1 <- barplot(var_mvp,
               main = "VaR 95% - modele MVP",
               ylab = "VaR (%)",
               las = 2,
               col = "coral",
               ylim = c(0, max(var_mvp) * 1.3),  # Zwiększone z 1.2 na 1.3
               names.arg = c("Static", "DCC-N", "DCC-S", "BEKK"),
               space = 0.5)  # Dodane odstępy między słupkami
text(bp1, var_mvp + 0.1, labels = round(var_mvp, 2), pos = 3, cex = 0.85, font = 2)  # Przesunięcie w górę
grid()

es_mvp <- abs(stats_overall_B[mvp_models, "ES 95%"] * 100)
bp2 <- barplot(es_mvp,
               main = "Expected Shortfall 95% - modele MVP",
               ylab = "ES (%)",
               las = 2,
               col = "darkred",
               ylim = c(0, max(es_mvp) * 1.3),  # Zwiększone z 1.2 na 1.3
               names.arg = c("Static", "DCC-N", "DCC-S", "BEKK"),
               space = 0.5)  # Dodane odstępy między słupkami
text(bp2, es_mvp + 0.1, labels = round(es_mvp, 2), pos = 3, cex = 0.85, font = 2)  # Przesunięcie w górę
grid()

par(mfrow = c(1, 1), mar = c(5, 4, 4, 2))

# Wykres 5: Krzywe kapitału - tylko MVP modele
mvp_only <- all_strategies_B[, mvp_models]
cum_mvp <- cumprod(1 + mvp_only)

plot(index(cum_mvp), cum_mvp[, "MVP_Static"],
     type = "l", lwd = 3, col = "black",
     ylim = range(cum_mvp),
     main = "Krzywe kapitału - porównanie modeli MVP",
     xlab = "Data",
     ylab = "Wartość portfela (start = 1)")
lines(index(cum_mvp), cum_mvp[, "DCC_Normal"], lwd = 3, col = "blue")
lines(index(cum_mvp), cum_mvp[, "DCC_Student"], lwd = 3, col = "darkgreen")
lines(index(cum_mvp), cum_mvp[, "BEKK"], lwd = 3, col = "purple")

legend("topleft",
       legend = c("MVP Static", "DCC-Normal", "DCC-Student", "BEKK"),
       col = c("black", "blue", "darkgreen", "purple"),
       lty = 1,
       lwd = 3,
       cex = 0.8)
grid()

# Wykres 6: Rozkład różnic zwrotów
par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))

for (model in mvp_models[-1]) {
  diff_ret <- all_strategies_B[, model] - all_strategies_B[, "MVP_Static"]
  
  if (sd(diff_ret, na.rm = TRUE) > 1e-10) {
    hist(diff_ret,
         breaks = 30,
         main = paste(gsub("_", " ", model), "vs Static"),
         xlab = "Różnica zwrotów",
         col = "lightblue",
         border = "white")
    abline(v = 0, col = "red", lwd = 2, lty = 2)
    abline(v = mean(diff_ret, na.rm = TRUE), col = "darkblue", lwd = 2)
    grid()
  } else {
    plot(1, type = "n", xlim = c(0, 1), ylim = c(0, 1),
         xlab = "", ylab = "", main = paste(model, "vs Static"))
    text(0.5, 0.5, "Modele identyczne", cex = 1.2)
  }
}

par(mfrow = c(1, 1), mar = c(5, 4, 4, 2))

# ================================================================================
# 11. PODSUMOWANIE KOŃCOWE
# ================================================================================

cat("\n================================================================================\n")
cat("                          PODSUMOWANIE CZĘŚCI B                                 \n")
cat("================================================================================\n\n")

cat(">>> PORÓWNANIE MODELI MVP:\n\n")

summary_table <- data.frame(
  Model = mvp_models,
  Sharpe = round(stats_overall_B[mvp_models, "Sharpe"], 4),
  Srednia_roczna = round(stats_overall_B[mvp_models, "Śr. roczna %"], 2),
  SD_roczna = round(stats_overall_B[mvp_models, "SD roczna %"], 2),
  VaR_95 = round(abs(stats_overall_B[mvp_models, "VaR 95%"]) * 100, 2),
  ES_95 = round(abs(stats_overall_B[mvp_models, "ES 95%"]) * 100, 2)
)

rownames(summary_table) <- NULL
print(summary_table)

cat("\n>>> GŁÓWNE WNIOSKI:\n")

best_model <- mvp_models[which.max(stats_overall_B[mvp_models, "Sharpe"])]
cat(sprintf("  1. Najwyższy Sharpe Ratio: %s (%.4f)\n", 
            best_model, 
            max(stats_overall_B[mvp_models, "Sharpe"])))

lowest_risk <- mvp_models[which.min(stats_overall_B[mvp_models, "SD roczna %"])]
cat(sprintf("  2. Najniższe ryzyko roczne: %s (%.2f%%)\n", 
            lowest_risk, 
            min(stats_overall_B[mvp_models, "SD roczna %"])))

lowest_es <- mvp_models[which.min(abs(stats_overall_B[mvp_models, "ES 95%"]))]
cat(sprintf("  3. Najniższy ES 95%%: %s (%.2f%%)\n", 
            lowest_es, 
            min(abs(stats_overall_B[mvp_models, "ES 95%"])) * 100))

