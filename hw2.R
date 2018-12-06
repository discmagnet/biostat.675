# Homework 2 ------------------------------------------------------------------------------
#
# Problem 1 - Breast Cancer Study

hw1 <- read_delim("~/WORKING_DIRECTORIES/biostat.675/Breast_cancer_Table_1_2_Collet.txt",
                  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
colnames(hw1) <- c("i","Xi","Di","Zi")
library(dplyr)
library(tidyr)

# (a) Compute the Nelson-Aalen estimator of the Cumulative Hazard function.

hw1a <- arrange(hw1, Xi)
hw1a$atRisk <- c(45:24,24,22:5,5,3:1)
hw1b <- subset(hw1a, Di == 1)

chf <- vector()
chf[1] <- 1/hw1b$atRisk[1]
for(i in 2:nrow(hw1b)){
  chf[i] <- chf[i-1] + 1/hw1b$atRisk[i]
}
hw1c <- mutate(hw1b,chf)

# (b) Compute 95% confidence intervals, assuming the estimated Cumulative Hazard function
# follows a Normal distribution.

var_chf <- vector()
var_chf[1] <- 1/hw1b$atRisk[1]^2
for(i in 2:nrow(hw1b)){
  var_chf[i] <- var_chf[i-1] + 1/hw1b$atRisk[i]^2
}
hw1d <- mutate(hw1c,var_chf)

chf_CI_lb <- vector()
chf_CI_ub <- vector()
for(i in 1:nrow(hw1b)){
  chf_CI_lb[i] <- chf[i] - 1.96*sqrt(var_chf[i])
  chf_CI_ub[i] <- chf[i] + 1.96*sqrt(var_chf[i])
}
hw1e <- mutate(hw1d,chf_CI_lb,chf_CI_ub)

# (c) Compute 95% confidence intervals, assuming the log of the estimated Cumulative Hazard
# function follows a Normal distribution.

chf_CI_log_lb <- vector()
chf_CI_log_ub <- vector()
for(i in 1:nrow(hw1b)){
  chf_CI_log_lb[i] <- exp(log(chf[i]) - 1.96*sqrt(var_chf[i]/chf[i]^2))
  chf_CI_log_ub[i] <- exp(log(chf[i]) + 1.96*sqrt(var_chf[i]/chf[i]^2))
}
hw1f <- mutate(hw1e,chf_CI_log_lb,chf_CI_log_ub)

# (d) Estimate S(t) using the Nelson-Aalen procedure, then compute 95% confidence intervals
# for S(t), assuming the log of the estimated Cumulative Hazard function follows a
# Normal distribution.

surv_CI_log_lb <- vector()
surv_CI_log_ub <- vector()
surv <- vector()
for(i in 1:nrow(hw1b)){
  surv_CI_log_lb[i] <- exp(-exp(log(chf[i]) + 1.96*sqrt(var_chf[i]/chf[i]^2)))
  surv_CI_log_ub[i] <- exp(-exp(log(chf[i]) - 1.96*sqrt(var_chf[i]/chf[i]^2)))
  surv[i] <- exp(-exp(log(chf[i])))
}
hw1g <- mutate(hw1f,surv_CI_log_lb,surv_CI_log_ub,surv)

# (e) Where possible, estimate the median and quartiles of the survival distribution.

hw1g_25 <- subset(hw1g,surv <= 1 - 0.25)
t_25 <- min(hw1g_25$Xi)

hw1g_50 <- subset(hw1g,surv <= 1 - 0.5)
t_50 <- min(hw1g_50$Xi)

hw1g_75 <- subset(hw1g,surv <= 1 - 0.75)
t_75 <- min(hw1g_75$Xi)
