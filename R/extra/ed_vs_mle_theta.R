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

discrepancy_envelope <- function(residual, lower, upper) {
  re <- ifelse(
    residual < lower, abs(residual - lower),
    ifelse(residual > upper, abs(residual - upper), 0.0)
  )

  sqrt(mean(re^2))
}

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

fit_negbin <- glm(
  cases ~ gender + age + pcr + sector + week + age:vp + offset(log(n)),
  data = data_positive,
  family = MASS::negative.binomial(theta = 20, link = "log")
)

new_responses <- simulate(fit_negbin, nsim = 10000L, seed = 505)[1:5000]

# This will be used to speed up the model refit
generator_params <- coef(fit_negbin)

env_discrepancy_theta <- function(theta, y, n_envelopes) {
  stopifnot(theta > 0)
  # HACK: can't use `y` directly in the formula, so change the data locally.
  # The data_positive outside of the function scope won't change.
  data_positive$cases <- y
  fit <- glm(
    cases ~ gender + age + pcr + sector + week + age:vp +
      offset(log(n)),
    data = data_positive,
    family = MASS::negative.binomial(theta = theta),
    start = generator_params
  )

  # It is necessary to define refit_fn because of limitation with
  # the current implementation of the default refit method.
  refit_fn <- function(.y, ...) {
    glm(
      .y ~ gender + age + pcr + sector + week + age:vp +
        offset(log(n)),
      data = data_positive,
      family = MASS::negative.binomial(theta = theta),
      start = generator_params
    )
  }

  # Compute several discrepancies from different envelope plots
  discrepancies <- replicate(n_envelopes, {
    envel <- asympDiag::envelope(fit,
      nsim = 100, residual_fn = residual_exclusion,
      refit_fn = refit_fn,
      plot.it = FALSE
    )

    discrepancy_envelope(envel$observed, envel$lower, envel$upper)
  })

  # Choose the median to avoid picking simulations that were too good/bad.
  result <- median(discrepancies)
  cat(sprintf("theta=%.5f\tED=%.8f\n", theta, result))
  result
}

model_metrics <- function(object, test_coefficients) {
  b <- coef(object)
  vc <- vcov(object, dispersion = 1)

  num_diff <- b - test_coefficients

  wald_joint_stat <- c(num_diff %*% solve(vc) %*% num_diff)
  wald_joint_pv <- 1 - pchisq(wald_joint_stat, length(test_coefficients))

  wald_uni_stat <- (num_diff^2) / diag(vc)
  wald_uni_pv <- 1 - pchisq(wald_uni_stat, 1)
  names(wald_uni_pv) <- paste0("wald_uni_pv_", names(wald_uni_pv))

  c(wald_joint_pv = wald_joint_pv, wald_uni_pv)
}

result_output <- "ed-ml-bn20-sim.csv"

# If the CSV file already exists, read the done simulation IDs; otherwise, create the file with a header.
if (file.exists(result_output)) {
  done_sims <- read.csv(result_output)$sim
} else {
  done_sims <- numeric()
}

todo_sims <- setdiff(seq_along(new_responses), done_sims)

for (i in todo_sims) {
  cat(sprintf("Sim: %d\n", i))

  # Ensure reproducible results
  set.seed(505)
  # Use GSS to find optimal theta.
  # I'm not looking for absolute precision.
  # If it finds a theta close to optimal is fine.
  theta_ed <- optimize(env_discrepancy_theta,
    y = new_responses[[i]], n_envelopes = 3,
    lower = 1, upper = 1000, tol = 1
  )$minimum

  negbin_ed <- glm(
    new_responses[[i]] ~ gender + age + pcr + sector + week + age:vp +
      offset(log(n)),
    data = data_positive,
    family = MASS::negative.binomial(theta = theta_ed),
    start = generator_params
  )

  negbin_ml <- glm.nb(
    new_responses[[i]] ~ gender + age + pcr + sector + week + age:vp +
      offset(log(n)),
    data = data_positive,
    start = generator_params
  )

  ed_metrics <- model_metrics(negbin_ed, test_coefficients = generator_params)
  names(ed_metrics) <- paste0("ed_", names(ed_metrics))
  ml_metrics <- model_metrics(negbin_ml, test_coefficients = generator_params)
  names(ml_metrics) <- paste0("ml_", names(ml_metrics))

  metrics <- c(
    ed_theta = theta_ed, ml_theta = negbin_ml$theta,
    ed_metrics, ml_metrics
  )
  csv_line <- c(sim = i, metrics)

  if (!file.exists(result_output)) {
    write.table(t(csv_line),
      sep = ",",
      file = result_output,
      row.names = FALSE,
      col.names = TRUE
    )
  } else {
    write.table(t(csv_line),
      sep = ",",
      file = result_output,
      row.names = FALSE,
      col.names = FALSE,
      append = TRUE
    )
  }
}

simulation_results <- read.csv(result_output)
significance_levels <- c(.01, .05, .10)

pv_ed_columns <- c("ed_wald_joint_pv")
pv_ml_columns <- c("ml_wald_joint_pv")

rejec_rates_ed <- sapply(
  significance_levels,
  function(alpha) mean(simulation_results[, pv_ed_columns] <= alpha)
)
rejec_rates_ml <- sapply(
  significance_levels,
  function(alpha) mean(simulation_results[, pv_ml_columns] <= alpha)
)
ed_row <- c(
  " " = "ED",
  "Bias" = mean(simulation_results$ed_theta) - 20,
  "MSE" = mean((simulation_results$ed_theta - 20)^2),
  rejec_rates_ed
)
ml_row <- c(
  " " = "ML",
  "Bias" = mean(simulation_results$ml_theta) - 20,
  "MSE" = mean((simulation_results$ml_theta - 20)^2),
  rejec_rates_ml
)

rbind(ed_row, ml_row) |>
  knitr::kable("latex", booktabs = TRUE, row.names = FALSE) |>
  print()
