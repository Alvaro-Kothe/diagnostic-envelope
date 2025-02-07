library("MASS") # To fit negative binomial models

get_measures <- function(obj) {
  out <- summary(obj, dispersion = 1)$coefficients
  out[, 4] <- out[, 4] * 100 # Multply p-value by 100
  out
}

build_table_coefs <- function(models_list) {
  measures_models <- lapply(models_list, get_measures)
  cols_idx <- c(1, 2, 4)
  n <- nrow(measures_models[[1]])
  nmodels <- length(models_list)
  nmeasures <- length(cols_idx)
  cols <- nmodels * nmeasures
  out <- matrix(NA_real_, nrow = n, ncol = cols, dimnames = list(rownames(measures_models[[1]])))

  j <- 1
  for (measure_col in cols_idx) {
    for (model_measure in measures_models) {
      out[, j] <- model_measure[, measure_col]
      j <- j + 1
    }
  }
  out
}

# Load data
data_url <- "https://raw.githubusercontent.com/yairgoldy/BNT162b2_waning_immunity/2fe12b0ac5aeab12f02cd65c5a3b3c13d506bc83/pos_data_days11-31_7.csv" # nolint
data_positive <- read.csv(data_url)
colnames(data_positive) <- c("vp", "gender", "age", "pcr", "sector", "week", "n", "rate1k")
data_positive$cases <- with(data_positive, round(n * rate1k / 1000))
data_positive$vp <- factor(data_positive$vp, c(
  "JanB", "FebA", "FebB", "MarA", "MarB", "Apr", "May"
))
data_positive$sector <- factor(data_positive$sector, c(
  "General Jewish", "Arab", "ultra-Orthodox Jewish"
))

fit_poisson <- glm(
  cases ~ gender + age + pcr + sector + week + age:vp + offset(log(n)),
  data = data_positive, family = poisson(link = "log")
)

fit_negbin <- glm(
  cases ~ gender + age + pcr + sector + week + age:vp + offset(log(n)),
  data = data_positive,
  family = MASS::negative.binomial(theta = 20, link = "log")
)

fit_list <- list(fit_poisson, fit_negbin)
models_names <- c("Poisson", "NB20")

est_table <- build_table_coefs(fit_list)

est_table |>
  knitr::kable("latex", col.names = rep(models_names, 3L)) |>
  cat()
