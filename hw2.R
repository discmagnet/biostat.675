# Problem 1 -----------------------------------------------------------------------

library(haven)
library(dplyr)
library(ggfortify)
library(survival)
library(ggplot2)
library(broom)

anemia2 <- read_sas("~/WORKING_DIRECTORIES/biostat.675/anemia2.sas7bdat")

# (a) Carry out separate two-sided log rank and Wilcoxon tests in order
#     to determine which of the following predict time until GVHD.

#     List your results in a small table listing covariate, test, test
#     statistic, and p-value.

#     In carrying out each test, plot the Survival Function, Cumulative
#     Hazard Function, and Log Cumulative Hazard Function estimators.

#     (i) age group (use categories <= 19 and >= 20)

anemia2 <- mutate(anemia2, under20 = 1*(age <= 19))

# Logrank Test based on Age Group (rho set to 0)
LRT_age <- survdiff(Surv(obs_time,GVHD) ~ under20, data = anemia2, rho = 0)
# Wilcoxon Test based on Age Group (rho set to 1)
WXT_age <- survdiff(Surv(obs_time,GVHD) ~ under20, data = anemia2, rho = 1)

# Plot of Survival Function
fit_age <- survfit(Surv(obs_time,GVHD) ~ under20, data = anemia2)
plot01 <- autoplot(fit_age, conf.int = FALSE) +
  ggtitle("Survival Function by Age Group")
# Plot of Cumulative Hazard Function
fit_age$surv <- -log(fit_age$surv) # transforms survival probabilities
plot02 <- autoplot(fit_age, conf.int = FALSE, surv.connect = FALSE) +
  ggtitle("Cumulative Hazard Function by Age Group") +
  ylab("-log(surv)")
# Plot of Log Cumulative Hazard Function
fit_age$surv <- log(fit_age$surv) # transforms cumulative hazard probabilities
plot03 <- autoplot(fit_age, conf.int = FALSE, surv.connect = FALSE) +
  ggtitle("Log Cumulative Hazard Function by Age Group") +
  ylab("log(-log(surv))")

plot01
plot02
plot03

#     (ii) LAF

# Logrank Test based on LAF (rho set to 0)
LRT_laf <- survdiff(Surv(obs_time,GVHD) ~ LAF, data = anemia2, rho = 0)
# Wilcoxon Test based on LAF (rho set to 1)
WXT_laf <- survdiff(Surv(obs_time,GVHD) ~ LAF, data = anemia2, rho = 1)

# Plot of Survival Function
fit_laf <- survfit(Surv(obs_time,GVHD) ~ LAF, data = anemia2)
plot04 <- autoplot(fit_laf, conf.int = FALSE) +
  ggtitle("Survival Function by LAF")
# Plot of Cumulative Hazard Function
fit_laf$surv <- -log(fit_laf$surv) # transforms survival probabilities
plot05 <- autoplot(fit_laf, conf.int = FALSE, surv.connect = FALSE) +
  ggtitle("Cumulative Hazard Function by LAF") +
  ylab("-log(surv)")
# Plot of Log Cumulative Hazard Function
fit_laf$surv <- log(fit_laf$surv) # transforms cumulative hazard probabilities
plot06 <- autoplot(fit_laf, conf.int = FALSE, surv.connect = FALSE) +
  ggtitle("Log Cumulative Hazard Function by LAF") +
  ylab("log(-log(surv))")

plot04
plot05
plot06

#     (iii) use of both CSP and MTX (versus not)

# Logrank Test based on CSP_MTX (rho set to 0)
LRT_csp <- survdiff(Surv(obs_time,GVHD) ~ CSP_MTX, data = anemia2, rho = 0)
# Wilcoxon Test based on CSP_MTX (rho set to 1)
WXT_csp <- survdiff(Surv(obs_time,GVHD) ~ CSP_MTX, data = anemia2, rho = 1)

# Plot of Survival Function
fit_cm <- survfit(Surv(obs_time,GVHD) ~ CSP_MTX, data = anemia2)
plot07 <- autoplot(fit_cm, conf.int = FALSE) +
  ggtitle("Survival Function by CSP_MTX")
# Plot of Cumulative Hazard Function
fit_cm$surv <- -log(fit_cm$surv) # transforms survival probabilities
plot08 <- autoplot(fit_cm, conf.int = FALSE, surv.connect = FALSE) +
  ggtitle("Cumulative Hazard Function by CSP_MTX") +
  ylab("-log(surv)")
# Plot of Log Cumulative Hazard Function
fit_cm$surv <- log(fit_cm$surv) # transforms cumulative hazard probabilities
plot09 <- autoplot(fit_cm, conf.int = FALSE, surv.connect = FALSE) +
  ggtitle("Log Cumulative Hazard Function by CSP_MTX") +
  ylab("log(-log(surv))")

plot07
plot08
plot09

# (c) Carry out a Logrank Test of CSP and MTX, but stratified on age group.

# Logrank Test based on CSP_MTX (rho set to 0)
LRT_csp_str01 <- survdiff(Surv(obs_time,GVHD) ~ CSP_MTX, data = anemia2, rho = 0,
                    subset = (under20 == 1))
# Logrank Test based on CSP_MTX (rho set to 0)
LRT_csp_str02 <- survdiff(Surv(obs_time,GVHD) ~ CSP_MTX, data = anemia2, rho = 0,
                    subset = (under20 == 0))

# Problem 2 -----------------------------------------------------------------------

anemia2 <- arrange(anemia2, obs_time)
anemia2$atRisk <- c(64:59,59,57,56,55,55,53,53,rep(51,4),47,47,45:24,rep(23,23))
anemia2 <- mutate(anemia2, hazard = GVHD/atRisk) # compute hazard function

b <- 4 # specify desired bandwidth
# plot smoothed hazard estimator for b < t < t_n - b
time <- c((10*b):(10*(max(anemia2$obs_time[anemia2$GVHD==1]-b))))/10

# compute smoothed hazard function
lambda <- vector()
for(i in 1:length(time)){
  s <- 0
  for(j in 1:nrow(anemia2)){
    x <- (time[i] - anemia2$obs_time[j])/b # compute (t-t_j)/b
    k <- 3/4*(1-x^2)*(abs(x) <= 1) # plug into kernel equation
    # (in our case, the Epanechnikov kernel function)
    h <- anemia2$hazard[j]
    # change in cumulative hazard is same as just the hazard
    s <- s + k*h 
  }
  lambda[i] <- s/b 
}

smoothed <- data.frame(time,lambda)

plot10 <- ggplot(data = smoothed, aes(x = time, y = lambda)) + 
  geom_line() +
  xlab("Time") +
  ylab("Hazard Rate") +
  ggtitle(paste("Smoothed Plot of the Hazard Function with Bandwidth = ", b))

plot10