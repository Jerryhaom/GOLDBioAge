.prepare_bioage_data <- function(d4, var) {
  if (!is.data.frame(d4)) stop("'d4' must be a data frame.")
  if (!is.character(var) || !length(var) || anyNA(var) ||
      anyDuplicated(var) || any(!nzchar(var))) {
    stop("'var' must contain unique, non-empty column names.")
  }
  if (!"age" %in% var) stop("'age' must be included in 'var'.")
  if (any(var %in% c("time", "status", "beta0")) ||
      any(grepl("`", var, fixed = TRUE))) {
    stop("'var' cannot include time, status, beta0, or names containing backticks.")
  }
  var <- c("age", setdiff(var, "age"))
  required <- c("time", "status", var)
  missing <- setdiff(required, names(d4))
  if (length(missing)) stop("Missing columns: ", paste(missing, collapse = ", "))
  if (anyDuplicated(names(d4))) stop("'d4' must have unique column names.")
  data <- as.data.frame(d4[, required, drop = FALSE])
  numeric_columns <- vapply(data, is.numeric, logical(1))
  numeric_columns["status"] <- is.numeric(data$status) || is.logical(data$status)
  if (!all(numeric_columns)) {
    stop("Model columns must be numeric; status may also be logical. Check: ",
         paste(names(data)[!numeric_columns], collapse = ", "))
  }
  if (any(vapply(data, function(x) any(is.infinite(x)), logical(1)))) {
    stop("Model columns must contain finite values or NA, not Inf.")
  }
  if (any(!is.na(data$status) & !data$status %in% c(0, 1))) {
    stop("'status' must use 0 = censored and 1 = event.")
  }
  if (any(data$time < 0, na.rm = TRUE) || any(data$age < 0, na.rm = TRUE)) {
    stop("'time' and 'age' must be non-negative.")
  }
  keep <- stats::complete.cases(data)
  data <- data[keep, , drop = FALSE]
  if (nrow(data) <= length(var) + 1L) stop("Too few complete observations to fit the model.")
  if (!any(data$status == 1)) stop("At least one observed event is required.")
  constant <- vapply(data[, var, drop = FALSE], function(x) length(unique(x)) < 2L,
                     logical(1))
  if (any(constant)) stop("Predictors must vary: ", paste(var[constant], collapse = ", "))
  if (any(!keep)) message("Omitted ", sum(!keep), " rows with missing model variables.")
  list(data = data, var = var, rows_used = which(keep), rows_omitted = which(!keep))
}

.bioage_formula <- function(var, offset_age = NULL) {
  terms <- paste0("`", var, "`")
  if (!is.null(offset_age)) {
    terms <- c(paste0("offset(", format(offset_age, digits = 17, scientific = TRUE),
                      " * age)"), terms)
  }
  stats::as.formula(paste("survival::Surv(time, status) ~", paste(terms, collapse = " + ")))
}

.with_bioage_seed <- function(seed, fn) {
  if (is.null(seed)) return(fn())
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  })
  set.seed(seed)
  fn()
}
