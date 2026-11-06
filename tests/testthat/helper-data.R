make_gompertz_data <- function() {
  set.seed(20261106)
  n <- 2000L
  age <- runif(n, 40, 80)
  marker_risk <- rnorm(n)
  marker_protective <- rnorm(n)
  noise <- rnorm(n)
  rate <- exp(-9 + 0.08 * age + 0.4 * marker_risk - 0.3 * marker_protective)
  event_time <- flexsurv::rgompertz(n, shape = 0.08, rate = rate)
  censor_time <- runif(n, 5, 20)
  data.frame(time = pmin(event_time, censor_time),
             status = as.integer(event_time <= censor_time),
             age, marker_risk, marker_protective, noise)
}

test_data <- make_gompertz_data()
test_vars <- c("age", "marker_risk", "marker_protective", "noise")
