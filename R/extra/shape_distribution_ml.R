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

fit_negbin_ml <- glm.nb(
  cases ~ gender + age + pcr + sector + week + age:vp + offset(log(n)),
  data = data_positive,
  maxit = 1000
)

set.seed(505)
thetas_ml <- asympDiag::parametric_bootstrap(fit_negbin_ml,
  statistic = function(model) model$theta,
  nsim = 10000L
)$result

hist(thetas_ml, breaks = 20)
abline(v = fit_negbin_ml$theta, col = "red")
