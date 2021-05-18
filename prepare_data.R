## keep only complete observations

setDT(d)
d[, sid := as.character(sid)]
d <- na.omit(d)

## both seizure types

d[, both_seizure_types := (electrographic_seizure & clinical_seizure)]

## aEEG mean

d[, aEEG_mean := 1 / 14 * (aEEG_p01 + aEEG_p02 + aEEG_p03 + aEEG_p04 +
                           aEEG_p05 + aEEG_p06 + aEEG_p07 + aEEG_p08 +
                           aEEG_p09 + aEEG_p10 + aEEG_p11 + aEEG_p12 +
                           aEEG_p13 + aEEG_p14)]

## aEEG maximum

d[, aEEG_max := pmax(aEEG_p01, aEEG_p02, aEEG_p03, aEEG_p04,
                     aEEG_p05, aEEG_p06, aEEG_p07, aEEG_p08,
                     aEEG_p09, aEEG_p10, aEEG_p11, aEEG_p12,
                     aEEG_p13, aEEG_p14)]

## aEEG differences

d[, `:=`(d00 = aEEG_p01,
         d01 = aEEG_p02 - aEEG_p01, d02 = aEEG_p03 - aEEG_p02,
         d03 = aEEG_p04 - aEEG_p03, d04 = aEEG_p05 - aEEG_p04,
         d05 = aEEG_p06 - aEEG_p05, d06 = aEEG_p07 - aEEG_p06,
         d07 = aEEG_p08 - aEEG_p07, d08 = aEEG_p09 - aEEG_p08,
         d09 = aEEG_p10 - aEEG_p09, d10 = aEEG_p11 - aEEG_p10,
         d11 = aEEG_p12 - aEEG_p11, d12 = aEEG_p13 - aEEG_p12,
         d13 = aEEG_p14 - aEEG_p13)]

## aEEG measured - predicted

d[, c01 := aEEG_p01]
d[, c02 := aEEG_p02 - predict(glm(aEEG_p02 ~ aEEG_p01, data = d))]
d[, c03 := aEEG_p03 - predict(glm(aEEG_p03 ~ aEEG_p01 + aEEG_p02, data = d))]
d[, c04 := aEEG_p04 - predict(glm(aEEG_p04 ~ aEEG_p01 + aEEG_p02 + aEEG_p03, data = d))]
d[, c05 := aEEG_p05 - predict(glm(aEEG_p05 ~ aEEG_p01 + aEEG_p02 + aEEG_p03 + aEEG_p04, data = d))]
d[, c06 := aEEG_p06 - predict(glm(aEEG_p06 ~ aEEG_p01 + aEEG_p02 + aEEG_p03 + aEEG_p04 +
                                             aEEG_p05, data = d))]
d[, c07 := aEEG_p07 - predict(glm(aEEG_p07 ~ aEEG_p01 + aEEG_p02 + aEEG_p03 + aEEG_p04 +
                                             aEEG_p05 + aEEG_p06, data = d))]
d[, c08 := aEEG_p08 - predict(glm(aEEG_p08 ~ aEEG_p01 + aEEG_p02 + aEEG_p03 + aEEG_p04 +
                                             aEEG_p05 + aEEG_p06 + aEEG_p07, data = d))]
d[, c09 := aEEG_p09 - predict(glm(aEEG_p09 ~ aEEG_p01 + aEEG_p02 + aEEG_p03 + aEEG_p04 +
                                             aEEG_p05 + aEEG_p06 + aEEG_p07 + aEEG_p08, data = d))]
d[, c10 := aEEG_p10 - predict(glm(aEEG_p10 ~ aEEG_p01 + aEEG_p02 + aEEG_p03 + aEEG_p04 +
                                             aEEG_p05 + aEEG_p06 + aEEG_p07 + aEEG_p08 + aEEG_p09, data = d))]
d[, c11 := aEEG_p11 - predict(glm(aEEG_p11 ~ aEEG_p01 + aEEG_p02 + aEEG_p03 + aEEG_p04 +
                                             aEEG_p05 + aEEG_p06 + aEEG_p07 + aEEG_p08 + aEEG_p09 +
                                             aEEG_p10, data = d))]
d[, c12 := aEEG_p12 - predict(glm(aEEG_p12 ~ aEEG_p01 + aEEG_p02 + aEEG_p03 + aEEG_p04 +
                                             aEEG_p05 + aEEG_p06 + aEEG_p07 + aEEG_p08 + aEEG_p09 +
                                             aEEG_p10 + aEEG_p11, data = d))]
d[, c13 := aEEG_p13 - predict(glm(aEEG_p13 ~ aEEG_p01 + aEEG_p02 + aEEG_p03 + aEEG_p04 +
                                             aEEG_p05 + aEEG_p06 + aEEG_p07 + aEEG_p08 + aEEG_p09 +
                                             aEEG_p10 + aEEG_p11 + aEEG_p12, data = d))]
d[, c14 := aEEG_p14 - predict(glm(aEEG_p14 ~ aEEG_p01 + aEEG_p02 + aEEG_p03 + aEEG_p04 +
                                             aEEG_p05 + aEEG_p06 + aEEG_p07 + aEEG_p08 + aEEG_p09 +
                                             aEEG_p10 + aEEG_p11 + aEEG_p12 + aEEG_p13, data = d))]

## aEEG slope coefficient

d_long <- melt(d, id.vars = c("sid", "outcome"),
                    measure.vars = paste0("aEEG_p", formatC(1:14, width = 2, flag = "0")),
                    na.rm = TRUE,
                    variable.name = "period",
                    value.name = "aeeg")
d_long[, periodN := as.numeric(substring(period, 7))]
slope <- list()
for (s in unique(d_long$sid)) {
  dx <- d_long[sid == s]
  model_fit <- glm(aeeg ~ periodN, data = dx)
  slope[[s]] <- data.table(slope = coef(model_fit)["periodN"])
}
slope <- rbindlist(slope, idcol = "sid")
slope[, sid := as.character(sid)]
d[slope, aEEG_slope := slope, on = "sid"]
