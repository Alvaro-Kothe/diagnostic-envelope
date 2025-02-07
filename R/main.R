# Install the package asympDiag for the analysis
# if (!require("pak")) install.packages("pak")
# pak::pak("Alvaro-Kothe/asympDiag")

library("MASS") # To fit negative binomial models.

# Save location.
envelope_location <- file.path("figures", "envelope")
ecdf_location <- file.path("figures", "ecdf")

envelope_nsim <- 100L
wald_nsim <- 10000L
run_wald_sim <- tolower(Sys.getenv("RUN_WALD_SIM", "0")) %in% c("1", "t", "true")

# Load data.
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

# Define function to compute the deletion residual.
residual_exclusion <- function(object, dispersion = 1) {
  # We define the dispersion parameter as 1 because we are only
  # dealing with Poisson and negative binomial models.
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

# Let's fit the poisson model.
fit_poisson <- glm(
  cases ~ gender + age + pcr + sector + week + age:vp + offset(log(n)),
  data = data_positive, family = poisson(link = "log")
)

# Create envelope plot.
set.seed(505)
pdf(file.path(envelope_location, "poisson.pdf"))
asympDiag::envelope(fit_poisson,
  nsim = envelope_nsim, residual_fn = residual_exclusion
)
dev.off()

# Let's fit the negative binomial model.

# We will define the shape parameter as 20 instead of estimating it,
# because simulations using this dataset showed that this parameter
# tends to be overestimated.
shape_parameter <- 20
fit_negbin <- glm(
  cases ~ gender + age + pcr + sector + week + age:vp + offset(log(n)),
  data = data_positive,
  family = MASS::negative.binomial(theta = shape_parameter, link = "log")
)

# Create envelope plot.
set.seed(505)
pdf(file.path(envelope_location, "negbin.pdf"))
asympDiag::envelope(fit_negbin,
  nsim = envelope_nsim, residual_fn = residual_exclusion
)
dev.off()

# Verify the asymptotic approximation for Wald test.
if (run_wald_sim) {
  # First we consider that the NB20 model is well defined.
  nb20_data <- simulate(fit_negbin, nsim = wald_nsim, seed = 505)
  pdf(file.path(ecdf_location, "negbin.pdf"))
  asympDiag::simulate_wald_pvalues(
    fit_negbin,
    test_coefficients = coef(fit_negbin),
    responses = nb20_data, plot.it = FALSE,
    vcov_fn = function(obj) vcov(obj, dispersion = 1)
  ) |>
    # With these arguments, the plot will only show the joint Wald test;
    # don't show Kolmogorov-Smirnov test results;
    # display all rejection rates for signif. levels of 1%, 5% and 10%;
    # and don't show plot legend.
    plot(
      which = Inf, ks_test = FALSE,
      discrepancy_tol = -1, uniform_legend = FALSE
    )
  dev.off()
  # Secondly, we verify if the Poisson Wald tests are reliable when
  # the data comes from a negative binomial with shape parameter 20.
  pdf(file.path(ecdf_location, "resp-negbin-fit-poisson.pdf"))
  asympDiag::simulate_wald_pvalues(
    fit_poisson,
    test_coefficients = coef(fit_negbin),
    responses = nb20_data, plot.it = FALSE,
    vcov_fn = function(obj) vcov(obj, dispersion = 1)
  ) |>
    plot(
      which = Inf, ks_test = FALSE,
      discrepancy_tol = -1, uniform_legend = FALSE
    )
  dev.off()
}
