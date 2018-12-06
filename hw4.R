# Homework 4 --------------------------------------------------------------------
#
# Data are available on n = 43 bone marrow transplant (BMT) patients treated at
# the Ohio State University BMT Unit. All patients had either Hodgkin's disease
# (HOD) or non-Hodgkin's lymphoma (NHL). Patients were given either an allogenic
# transplant from a sibling donor or an autogenic transplant (i.e., their own
# marrow was cleansed and transplanted). Data are also available on Karnofsky
# score and time (months) until BMT. Note that, although time until transplant
# is measured in months, observation time is measured in days.

library(haven)
bmt <- read_sas("~/WORKING_DIRECTORIES/biostat.675/bmt_lymphoma_1.sas7bdat")

bmt$Wait_time_yrs <- bmt$Wait_time/12

# (a) Fit a Weibull regression model containing all main effects, with waiting
#     time measured in years. List the parameter estimates and corresponding
#     estimated standard errors.

library(survival)
library(broom)
model <- survreg(formula = Surv(obs_time,dead) ~ Karnof + 
                   Wait_time_yrs + factor(Tx_type) + factor(cancer),
                 data = bmt,
                 dist = "weibull")
est <- tidy(model)

# (b) Provide a frequency table for Karnofsky score.



# (c) Estimate and interpret the Karnofsky effect from a PH perspective; pay
#     attention to the table you produced in (b). Compute a 95% confidence 
#     interval.

# (d) Again, estimate and interpret the Karnofsky effect, but from an AFT
#     perspective; pay attention to the table you produced in (b). Compute a 95%
#     confidence interval.

# (e) Test H0: sigma == 1 versus H1: sigma =! 1 using the Wald test. Could the
#     exponential model be used for this data set?