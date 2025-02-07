library("MASS") # To fit negative binomial models

wald_test <- function(est, vc, idx = seq_along(est), null = rep(0, length(idx))) {
  est <- est[idx]
  vc <- vc[idx, idx]
  p <- length(idx)
  num <- est - null
  estat <- num %*% solve(vc) %*% num
  pv <- 1 - pchisq(estat, p)
  c(estat = estat, pvalue = pv)
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

index_age_vp <- 11:28
compute_wald_stat <- function(object, index) {
  b <- coef(object)
  vc <- vcov(object, dispersion = 1)
  wald_test(b, vc, index)
}

fit_poisson <- glm(
  cases ~ gender + age + pcr + sector + week + age:vp + offset(log(n)),
  data = data_positive, family = poisson(link = "log")
)

fit_negbin <- glm(
  cases ~ gender + age + pcr + sector + week + age:vp + offset(log(n)),
  data = data_positive,
  family = MASS::negative.binomial(theta = 20, link = "log")
)

sapply(list("Poisson" = fit_poisson, "NB20" = fit_negbin), compute_wald_stat, index = index_age_vp)
