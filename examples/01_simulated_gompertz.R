# Run after installing GOLDBioAge: source("examples/01_simulated_gompertz.R")
library(GOLDBioAge)

# 1. Simulate baseline predictors and Gompertz event times (follow-up in years).
set.seed(20261106)
n <- 2000L
age <- runif(n, 40, 80)
marker_risk <- rnorm(n)
marker_protective <- rnorm(n)
noise <- rnorm(n)
shape <- 0.08
individual_rate <- exp(-9 + 0.08 * age + 0.4 * marker_risk - 0.3 * marker_protective)
event_time <- flexsurv::rgompertz(n, shape = shape, rate = individual_rate)
censor_time <- runif(n, 5, 20)
simulated <- data.frame(
  time = pmin(event_time, censor_time),
  status = as.integer(event_time <= censor_time),
  age, marker_risk, marker_protective, noise
)

# 2. Split participants before fitting or feature selection.
train_id <- sample.int(n, size = floor(0.8 * n))
train <- simulated[train_id, , drop = FALSE]
test <- simulated[-train_id, , drop = FALSE]
variables <- c("age", "marker_risk", "marker_protective", "noise")

# 3. Fit on the development sample and inspect weights in years.
sim_fit <- gold_bioage(train, variables)
print(sim_fit$coefficients)

# 4. Score held-out predictors; no follow-up or event status is passed to scoring.
test$GOLDbioage <- predict_gold_bioage(sim_fit, test[, variables, drop = FALSE])
test$age_difference <- test$GOLDbioage - test$age
print(c(n_train = nrow(train), n_test = nrow(test),
        events_train = sum(train$status), events_test = sum(test$status),
        mean_train_age_difference = mean(sim_fit$GOLDBioage$GOLDbioage - train$age),
        mean_test_age_difference = mean(test$age_difference)))
stopifnot(length(test$GOLDbioage) == nrow(test), all(is.finite(test$GOLDbioage)))

# 5. Optional selection is performed only within the development sample.
if (requireNamespace("glmnet", quietly = TRUE) && requireNamespace("ggplot2", quietly = TRUE)) {
  sim_lasso <- gold_bioage(train, variables, feature_selection = TRUE,
                           selection_method = "lasso", nfolds = 5, seed = 123)
  sim_elasticnet <- gold_bioage(train, variables, feature_selection = TRUE,
                                selection_method = "elasticnet", alpha = 0.5,
                                nfolds = 5, seed = 123)
  print(sim_lasso$selected_vars)
  print(sim_elasticnet$selected_vars)
  if (interactive()) print(sim_lasso$cv_plot)
  test$GOLDbioage_lasso <- predict_gold_bioage(sim_lasso, test)
}

# These are artificial data and an algorithm demonstration, not clinical validation.
