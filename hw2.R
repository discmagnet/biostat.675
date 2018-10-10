# Problem 1

library(haven)
library(dplyr)
library(ggfortify)
library(survival)
library(ggplot2)

anemia2 <- read_sas("~/WORKING_DIRECTORIES/biostat.675/anemia2.sas7bdat")

# (a) Carry out separate two-sided log rank and Wilcoxon tests in order
#     to determine which of the following predict time until GVHD.

#     List your results in a small table listing covariate, test, test
#     statistic, and p-value.

#     In carrying out each test, plot the Survival Function, Cumulative
#     Hazard Function, and Log Cumulative Hazard Function estimators.

#     (i) age group (use categories <= 19 and >= 20)

anemia2 <- mutate(anemia2, under20 = 1*(age <= 19))

fit_age <- survfit(Surv(obs_time,GVHD) ~ under20, data = anemia2)
plot01 <- autoplot(fit_age, conf.int = FALSE) +
  ggtitle("Survival Function by Age Group")
fit_age$surv <- -log(fit_age$surv)
plot02 <- autoplot(fit_age, conf.int = FALSE, surv.connect = FALSE) +
  ggtitle("Cumulative Hazard Function by Age Group") +
  ylab("-log(surv)")
fit_age$surv <- log(fit_age$surv)
plot03 <- autoplot(fit_age, conf.int = FALSE, surv.connect = FALSE) +
  ggtitle("Log Cumulative Hazard Function by Age Group") +
  ylab("log(-log(surv))")

plot01
plot02
plot03

#     (ii) LAF

fit_laf <- survfit(Surv(obs_time,GVHD) ~ LAF, data = anemia2)
plot04 <- autoplot(fit_laf, conf.int = FALSE) +
  ggtitle("Survival Function by LAF")
fit_laf$surv <- -log(fit_laf$surv)
plot05 <- autoplot(fit_laf, conf.int = FALSE, surv.connect = FALSE) +
  ggtitle("Cumulative Hazard Function by LAF") +
  ylab("-log(surv)")
fit_laf$surv <- log(fit_laf$surv)
plot06 <- autoplot(fit_laf, conf.int = FALSE, surv.connect = FALSE) +
  ggtitle("Log Cumulative Hazard Function by LAF") +
  ylab("log(-log(surv))")

plot04
plot05
plot06

#     (iii) use of both CSP and MTX (versus not)

fit_cm <- survfit(Surv(obs_time,GVHD) ~ CSP_MTX, data = anemia2)
plot07 <- autoplot(fit_cm, conf.int = FALSE) +
  ggtitle("Survival Function by CSP_MTX")
fit_cm$surv <- -log(fit_cm$surv)
plot08 <- autoplot(fit_cm, conf.int = FALSE, surv.connect = FALSE) +
  ggtitle("Cumulative Hazard Function by CSP_MTX") +
  ylab("-log(surv)")
fit_cm$surv <- log(fit_cm$surv)
plot09 <- autoplot(fit_cm, conf.int = FALSE, surv.connect = FALSE) +
  ggtitle("Log Cumulative Hazard Function by CSP_MTX") +
  ylab("log(-log(surv))")

plot07
plot08
plot09



anemia2 <- arrange(anemia2, obs_time)
anemia2$atRisk <- c(64:59,59,57,56,55,55,53,53,rep(51,4),47,47,45:24,rep(23,23))
anemia2 <- subset(anemia2, GVHD == 1)

chf <- vector()
chf[1] <- 1/anemia2$atRisk[1]
for(i in 2:nrow(anemia2)){
  chf[i] <- chf[i-1] + 1/anemia2$atRisk[i]
}
anemia2 <- mutate(anemia2,chf)