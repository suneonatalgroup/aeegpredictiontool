## Prediction of long-term neurodevelopmental outcome

## Input: `d` data.frame or data.table with the following variables:
## - sid: unique subject ID
## - outcome: 0 for good outcome, 1 for poor outcome
## - aEEG_p01: aEEG score (1, 2, 3, 4 or 5) for the  1st time interval (00-06h)
## - aEEG_p02: aEEG score (1, 2, 3, 4 or 5) for the  2nd time interval (06-12h)
## - aEEG_p03: aEEG score (1, 2, 3, 4 or 5) for the  3rd time interval (12-18h)
## - aEEG_p04: aEEG score (1, 2, 3, 4 or 5) for the  4th time interval (18-24h)
## - aEEG_p05: aEEG score (1, 2, 3, 4 or 5) for the  5th time interval (24-30h)
## - aEEG_p06: aEEG score (1, 2, 3, 4 or 5) for the  6th time interval (30-36h)
## - aEEG_p07: aEEG score (1, 2, 3, 4 or 5) for the  7th time interval (36-42h)
## - aEEG_p08: aEEG score (1, 2, 3, 4 or 5) for the  8th time interval (42-48h)
## - aEEG_p09: aEEG score (1, 2, 3, 4 or 5) for the  9th time interval (48-54h)
## - aEEG_p10: aEEG score (1, 2, 3, 4 or 5) for the 10th time interval (54-60h)
## - aEEG_p11: aEEG score (1, 2, 3, 4 or 5) for the 11th time interval (60-66h)
## - aEEG_p12: aEEG score (1, 2, 3, 4 or 5) for the 12th time interval (66-72h)
## - aEEG_p13: aEEG score (1, 2, 3, 4 or 5) for the 13th time interval (72-78h)
## - aEEG_p14: aEEG score (1, 2, 3, 4 or 5) for the 14th time interval (78-84h)
## - swc: sleep-wake cycling coded as 0 for absent, 1 for abnormal, 2 for normal

## Example data structure:

if (FALSE) {
  d <- data.table::fread("
    S01 0 1 2 1 2 3 1 2 1 1 1 2 1 2 1 2
    S02 1 3 4 3 2 4 1 2 2 3 2 1 2 3 2 0
    S03 0 2 1 1 2 1 2 1 1 2 2 2 2 1 1 1
    S04 1 4 2 1 4 3 1 2 1 1 2 2 4 1 3 1
    S05 1 3 4 3 2 4 1 2 3 2 2 2 2 3 1 0
    S06 1 3 2 3 2 4 1 2 1 3 2 1 2 3 2 2
    S07 0 2 1 1 2 1 2 1 2 1 2 1 2 1 1 2
    S08 0 1 2 1 1 3 1 2 1 3 2 1 4 1 3 1
    S09 1 1 2 1 4 3 1 2 1 2 2 1 4 1 1 2
    S10 1 3 1 5 2 2 1 2 5 3 5 2 1 2 1 1
    S11 0 1 3 1 2 3 1 3 1 1 1 2 2 2 2 2
    S12 1 3 4 3 2 4 1 3 2 3 2 1 3 3 1 1
    S13 0 2 2 1 2 1 2 2 1 2 2 2 3 1 4 2
    S14 1 4 3 1 4 3 1 3 1 1 2 2 5 1 1 0
    S15 1 3 5 3 2 4 1 3 3 2 2 2 3 3 3 0
    S16 1 3 3 3 2 4 1 3 1 3 2 1 3 3 1 2
    S17 0 2 2 1 2 1 2 2 2 1 2 1 3 1 4 1
    S18 0 1 3 1 1 3 1 3 1 3 2 1 5 1 1 2
    S19 1 1 4 1 4 3 1 3 1 2 2 1 5 1 2 1
    S20 1 3 2 5 2 2 1 3 5 3 5 2 2 2 1 2
    S21 1 1 3 1 2 3 1 3 2 1 1 2 2 2 4 0
    S22 0 3 4 3 2 4 1 3 1 3 2 1 3 3 4 2
    S23 1 2 2 1 2 1 2 2 2 2 2 2 3 1 3 0
    S24 0 4 3 1 4 3 1 3 2 1 2 2 5 1 2 0
    S25 0 3 5 3 2 4 1 3 1 2 2 2 3 3 2 1
    S26 0 3 3 3 2 4 1 3 2 3 2 1 3 3 3 2
    S27 1 2 2 1 2 1 2 2 3 1 2 1 3 1 3 0
    S28 1 1 3 1 1 3 1 3 2 3 2 1 5 1 4 2
    S29 0 1 4 1 4 3 1 3 3 2 2 1 5 1 5 1
    S30 0 3 2 5 2 2 1 3 2 3 5 2 2 2 3 2
  ")
  data.table::setnames(d, c(
    "sid", "outcome",
    paste0("aEEG_p", formatC(1:14, width = 2, flag = "0")),
    "swc"))
}


## load libraries and prepare data

library(data.table)
library(pROC)
library(DescTools)
library(cvAUC)

source("helper_functions.R")
source("prepare_data.R")



## All measurements model
## --- all 14 aEEG scores were used as separate predictors

modelN <- 1
vars <- c("outcome", paste0("aEEG_p", formatC(1:14, width = 2, flag = "0")))
dx <- na.omit(d[, ..vars])
model_fit <- glm(outcome ~ ., data = dx, family = "binomial")
model_eval(modelN, model_fit, dx)


## Best measurement model
## --- in our case the best aEEG score was found at 18 hours

modelN <- 2
vars <- c("outcome", "aEEG_p03")
dx <- na.omit(d[, ..vars])
model_fit <- glm(outcome ~ ., data = dx, family = "binomial")
model_eval(modelN, model_fit, dx)


## Summary measurement (mean) method
## --- mean aEEG score of the 14 periods

modelN <- 3
vars <- c("outcome", "aEEG_mean")
dx <- na.omit(d[, ..vars])
model_fit <- glm(outcome ~ ., data = dx, family = "binomial")
model_eval(modelN, model_fit, dx)


## Summary measurement (maximum) method
## --- worst aEEG score for each individual summarised

modelN <- 4
vars <- c("outcome", "aEEG_max")
dx <- na.omit(d[, ..vars])
model_fit <- glm(outcome ~ ., data = dx, family = "binomial")
model_eval(modelN, model_fit, dx)


## Changes between subsequent measurements method
## --- calculation is based on the difference between the consecutive aEEG scores


modelN <- 5
vars <- c("outcome",  paste0("d", formatC(0:13, width = 2, flag = "0")))
dx <- na.omit(d[, ..vars])
model_fit <- glm(outcome ~ ., data = dx, family = "binomial")
model_eval(modelN, model_fit, dx)


## Conditional measurements method
## --- calculation is based on the difference between the measured and
## --- estimated aEEG scores using regression from the previous measurements

modelN <- 6
vars <- c("outcome",  paste0("c", formatC(1:14, width = 2, flag = "0")))
dx <- na.omit(d[, ..vars])
model_fit <- glm(outcome ~ ., data = dx, family = "binomial")
model_eval(modelN, model_fit, dx)


## Growth curve method
## --- calculation is based on the mean and the slope coefficient of aEEG scores
## --- using a linear regression curve for each individual

modelN <- 7
vars <- c("outcome", "aEEG_mean", "aEEG_slope")
dx <- na.omit(d[, ..vars])
model_fit <- glm(outcome ~ ., data = dx, family = "binomial")
model_eval(modelN, model_fit, dx)


## Growth curve method extended with the `swc` variable (7B model)

modelN <- "7B"
vars <- c("outcome", "aEEG_mean", "aEEG_slope", "swc")
dx <- na.omit(d[, ..vars])
model_fit <- glm(outcome ~ ., data = dx, family = "binomial")
model_eval(modelN, model_fit, dx)


## 7B model: 5-fold cross validation

modelN <- "7B"
vars <- c("outcome", "aEEG_mean", "aEEG_slope", "swc")
dx <- na.omit(d[, ..vars])
set.seed(1)
CV_AUC_CI(dx)


## 7B model: 5-fold cross validation repeated 1000 times

modelN <- "7B"
vars <- c("outcome", "aEEG_mean", "aEEG_slope", "swc")
dx <- na.omit(d[, ..vars])
set.seed(1)
auc_rep <- list()
for (i in 1:1000) auc_rep[[i]] <- CV_AUC_CI(dx, print = FALSE)
auc_rep <- rbindlist(auc_rep)
cat("\nAUC mean from 1000 repetitions:", mean(auc_rep$auc), "\n\n")
