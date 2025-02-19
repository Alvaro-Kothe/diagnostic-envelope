library("MASS")

data_url <- paste0(
  "https://raw.githubusercontent.com/",
  "yairgoldy/BNT162b2_waning_immunity/",
  "2fe12b0ac5aeab12f02cd65c5a3b3c13d506bc83/",
  "pos_data_days11-31_7.csv"
)
data_positive <- read.csv(data_url)
colnames(data_positive) <- c(
  "vp", "gender", "age", "pcr", "sector", "week", "n", "rate1k"
)
data_positive$cases <- with(data_positive, round(n * rate1k / 1000))
data_positive$vp <- factor(data_positive$vp, c(
  "JanB", "FebA", "FebB", "MarA", "MarB", "Apr", "May"
))
data_positive$sector <- factor(data_positive$sector, c(
  "General Jewish", "Arab", "ultra-Orthodox Jewish"
))

# Start: article code
residual_exclusion <- function(object, dispersion = 1) {
  deviance_residual <- residuals(object, type = "deviance")
  pearson_residual <- residuals(object, type = "pearson")
  h <- hatvalues(object)
  std_den <- sqrt(dispersion * (1 - h))
  standard_dev <- deviance_residual / std_den
  standard_pear <- pearson_residual / std_den
  residual_squared <-
    ((1 - h) * (standard_dev^2)) + (h * (standard_pear^2))
  sqrt(residual_squared)
}

discrepancy_envelope <- function(residual, lower, upper) {
  re <- ifelse(
    residual < lower, abs(residual - lower),
    ifelse(residual > upper, abs(residual - upper), 0.0)
  )

  sqrt(mean(re^2))
}

thetas <- c(5, 10, 20, 30, 50, 55.86, 100, 1000)
discrepancies <- numeric(length(thetas))

set.seed(505)
for (i in seq_along(thetas)) {
  fit <- glm(
    cases ~ gender + age + pcr + sector + week + age:vp +
      offset(log(n)),
    data = data_positive,
    family = MASS::negative.binomial(theta = thetas[[i]])
  )
  envel <- asympDiag::envelope(fit,
    nsim = 100, residual_fn = residual_exclusion,
    plot.it = FALSE
  )
  discrepancies[i] <-
    discrepancy_envelope(envel$observed, envel$lower, envel$upper)
}

cbind(thetas, discrepancies)
