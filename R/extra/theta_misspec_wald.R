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

fit_negbin <- glm(
  cases ~ gender + age + pcr + sector + week + age:vp + offset(log(n)),
  data = data_positive,
  family = MASS::negative.binomial(theta = 20, link = "log")
)

set.seed(505)
# Start: article code
thetas <- c(1, 5, 10, 20, 30, 50, 60, 100, 250, 1000)
significance_levels <- c(0.01, 0.05, 0.10)
new_responses <- simulate(fit_negbin, nsim = 10000L)
rejection_rates <-
  matrix(NA, nrow = length(thetas), ncol = length(significance_levels))

for (i in seq_along(thetas)) {
  wald_pvalues <- asympDiag::simulate_wald_pvalues(
    fit_negbin,
    responses = new_responses,
    family = negative.binomial(thetas[[i]]),
    vcov_fn = function(obj) vcov(obj, dispersion = 1),
    plot.it = FALSE
  )$pvalues_joint

  rejection_rates[i, ] <- sapply(
    significance_levels,
    function(alpha) mean(wald_pvalues <= alpha)
  )
}

cbind(thetas, rejection_rates)
# End: article code

cbind(thetas, rejection_rates) |>
  knitr::kable("latex", booktabs = TRUE, linesep = "")
