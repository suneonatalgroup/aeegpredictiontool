## rounding helper

roundC <- function(x, d = 0) {
  trimws(format(round(x, d), nsmall = d))
}


## ROC table

roc_table <- function(df, xv, outcome, step = 0.1) {
  x <- df[[xv]]
  b <- (df[[outcome]] %in% c(1, "Y", "Adverse"))
  if (length(unique(x)) == 2) {
    x1 <- max(x)
    x0 <- x1
  } else {
    x0 <- round(min(x), 1) - step
    x1 <- round(max(x), 1) + step
    if (step == 0.01) x0 <- 0
  }
  K <- seq(x1, x0, -step)
  l <- length(K)
  fp <- rep(NA, l)
  tn <- rep(NA, l)
  tn_ci <- rep(NA, l)
  tp <- rep(NA, l)
  tp_ci <- rep(NA, l)
  ppe <- rep(NA, l)
  ppe_ci <- rep(NA, l)
  npe <- rep(NA, l)
  npe_ci <- rep(NA, l)
  sum_b0 <- rep(NA, l)
  sum_b1 <- rep(NA, l)
  sum_t0 <- rep(NA, l)
  sum_t1 <- rep(NA, l)
  phi <- rep(NA, l)
  for (i in 1:l) {
    t <- (x >= K[i])
    fp[i]     <- sum( t & !b) / sum(!b)
    tn[i]     <- 1 - fp[i]
    tn_ci[i]  <- paste(roundC(binom.test(sum(!t & !b), sum(!b))$conf.int, 6), collapse = " - ")
    tp[i]     <- sum( t &  b) / sum( b)
    tp_ci[i]  <- paste(roundC(binom.test(sum( t &  b), sum( b))$conf.int, 6), collapse = " - ")
    ppe[i]    <- sum( b &  t) / sum( t)
    if (sum( t) > 0) ppe_ci[i] <- paste(roundC(binom.test(sum( b &  t), sum( t))$conf.int, 6), collapse = " - ")
    npe[i]    <- sum(!b & !t) / sum(!t)
    if (sum(!t) > 0) npe_ci[i] <- paste(roundC(binom.test(sum(!b & !t), sum(!t))$conf.int, 6), collapse = " - ")
    sum_b0[i] <- sum(!b)
    sum_b1[i] <- sum( b)
    sum_t0[i] <- sum(!t)
    sum_t1[i] <- sum( t)
    if (sd(t) > 0 && sd(b) > 0) {
      phi[i] <- cor(as.numeric(t), as.numeric(b), method = "pearson")
    }
  }
  data.table(K, fp, tp, d = tp - fp, tn, tn_ci, tp_ci,
             ppe, ppe_ci, npe, npe_ci,
             sum_b0, sum_b1, sum_t0, sum_t1, phi)
}


## ROC calculations

tab_opt_roc <- function(rt, ds, v1, v2, model = NULL, w = which.max(rt$d)) {
  ts0 <- rt$K[w]
  spec <- rt$tn[w]
  spec_ci <- rt$tn_ci[w]
  sens <- rt$tp[w]
  sens_ci <- rt$tp_ci[w]
  ppe <- rt$ppe[w]
  ppe_ci <- rt$ppe_ci[w]
  npe <- rt$npe[w]
  npe_ci <- rt$npe_ci[w]
  roc <- suppressMessages(roc(ds[[v1]], ds[[v2]]))
  aucci <- ci.auc(roc)
  if (!is.null(model)) nagr2 <- PseudoR2(model, which = "Nagelkerke") else nagr2 <- NA_real_
  phi <- rt$phi[w]
  dt <- data.table(" "            = c("AUC",
                                      "Nagelkerke R-sq",
                                      "Optimal threshold",
                                      "Sensitivity",
                                      "Specificity",
                                      "Pos Pred Value",
                                      "Neg Pred Value",
                                      "Phi coefficient"),
                   "Value"        = c(roundC(aucci[2], 6),
                                      roundC(nagr2, 6),
                                      ts0,roundC(sens, 6),
                                      roundC(spec, 6),
                                      roundC(ppe, 6),
                                      roundC(npe, 6),
                                      roundC(phi, 6)),
                   "Exact 95% CI" = c(paste(roundC(aucci[c(1, 3)], 6), collapse = " - "),
                                      "",
                                      "",
                                      sens_ci,
                                      spec_ci,
                                      ppe_ci,
                                      npe_ci,
                                      ""),
                   check.names = FALSE)
  if (nrow(rt) == 1) dt <- dt[4:8]
  if (is.null(model)) dt <- dt[` ` != "Nagelkerke R-sq"]
  print(dt)
  invisible(dt)
}


## Model evaluation

model_eval <- function(modelN, model_fit, dx) {
  cat("\nModel", modelN, "\n")
  dx[, paste0("m", modelN, "_pred") := predict(model_fit, type = c("response"))]
  rt <- roc_table(dx, paste0("m", modelN, "_pred"), "outcome", 0.01)
  opt <- tab_opt_roc(rt, dx, "outcome", paste0("m", modelN, "_pred"), model = model_fit)
  cat("\n")
  invisible(opt)
}


## Five-fold cross validation (https://github.com/ledell/cvAUC)

CV_AUC_CI <- function(df, print = TRUE) {
  train <- setDF(copy(df))
  y <- "outcome"
  V <- 5
  ## helpers
  cv_folds <- function(Y, V) {
    ## stratify by outcome
    w_Y0 <- which(Y == 0)
    w_Y1 <- which(Y == 1)
    Y0 <- split(sample(w_Y0), rep(1:V, length.out = length(w_Y0)))
    Y1 <- split(sample(w_Y1), rep(1:V, length.out = length(w_Y1)))
    folds <- list()
    for (v in 1:V) folds[[v]] <- c(Y0[[v]], Y1[[v]])
    folds
  }
  predict_logreg <- function(v, folds, train) {
    dx_train <- train[-folds[[v]], ]
    dx_test <- train[folds[[v]], ]
    model_fit <- glm(outcome ~ ., data = dx_train, family = "binomial")
    predict(model_fit, newdata = dx_test, type = c("response"))
  }
  ## create folds
  folds <- cv_folds(Y = train[, c(y)], V = V)
  ## generate CV predicted values
  predictions <- c()
  for (v in 1:V) predictions <- c(predictions, predict_logreg(v, folds, train))
  predictions[unlist(folds)] <- predictions
  names(predictions) <- as.character(seq_along(predictions))
  ## get CV AUC and 95% confidence interval
  res1 <- cvAUC(predictions = predictions, labels = train[, c(y)], folds = folds)
  res2 <- ci.cvAUC(predictions = predictions, labels = train[, c(y)],
                   folds = folds, confidence = 0.95)
  if (print) {
    cat("\nModel", modelN, "-", paste0(V, "-fold cross-validation"), "\n\n")
    cat("AUC mean:", round(res2$cvAUC, 4), "\n")
    cat("AUC 95% CI:", paste0("(", round(res2$ci[1], 4), ", ", round(res2$ci[2], 4), ")\n"))
    cat("AUC 95% CI radius:", paste0(round((res2$ci[2] - res2$ci[1]) / 2, 4), "\n"))
    cat("AUC range:", round(min(res1$fold.AUC), 4), "-", round(max(res1$fold.AUC), 4), "\n\n\n")
  }
  ## return
  invisible(data.table(auc = res2$cvAUC, ci_radius = (res2$ci[2] - res2$ci[1]) / 2))
}
