# Homework 5 --------------------------------------------------------------------
#
# Problem 1 - Data on n = 64 patients with severe anemia are contained in the SAS
# file anemia2.sas7bdat located in the "Data Sets" folder. The failure time of
# interest is time until graft-versus-host-disease (GVHD), measured in days. The
# variate CSP.MTX is a treatment indicator (0 = no, 1 = yes), while LAF is an
# indicator for being assigned to an airflow isolation room. AGE is recorded in
# years.
library(haven)
anemia2 <- read_sas("~/WORKING_DIRECTORIES/biostat.675/datasets/anemia2.sas7bdat")
# (a) Fit a Cox model with CSP.MTX (Zi), LAF (Li) and AGE (Ai) as covariates.
#     Interpret each of the hazard ratios.
library(survival)
model <- coxph(Surv(obs_time,GVHD) ~ CSP_MTX + LAF + age, data = anemia2)
summary(model)

# (b) Refit the model, with Ai replaced by Ai/5. Compare each parameter estimate
#     to that from (a) and comment on their similarity or differences.
library(dplyr)
anemia2b <- mutate(anemia2, age5 = age/5)
model2 <- coxph(Surv(obs_time,GVHD) ~ CSP_MTX + LAF + age5, data = anemia2b)
summary(model2)

# (c) Is the treatment effect (i.e., effect of Zi) different for subjects of
#     different ages? Carry out an appropriate Wald test.
model3 <- coxph(Surv(obs_time,GVHD) ~ CSP_MTX + LAF + age + CSP_MTX*age, data = anemia2)
summary(model3)

# (d) Re-evaluate the hypothesis from (c), but this time use a likelihood ratio
#     test. Compare your results to that obtained through the Wald test and
#     comment.
like_full <- model3$loglik[2]
like_rdcd <- model$loglik[2]
LRT <- 2*(like_full - like_rdcd)
pval_LRT <- 1 - pchisq(LRT, 1)

# (f) Is the effect of age linear? Support your response empirically by fitting
#     an appropriate main effects model and providing the appropriate plot.
summary(anemia2$age)
anemia2c <- mutate(anemia2b,
                   age_cat1 = 1*(age > 0)*(age <= 14),
                   age_cat2 = 1*(age > 14)*(age <= 21),
                   age_cat3 = 1*(age > 21)*(age <= 27.25),
                   age_cat4 = 1*(age > 27.25)*(age <= 42))
model5 <- coxph(Surv(obs_time, GVHD) ~ CSP_MTX + LAF + age_cat2 + age_cat3 + 
                  age_cat4,
                data = anemia2c)
summary(model5)
betas <- c(0,model5$coefficients[3:5])
names(betas)[1] <- "age_cat1"
mean1 <- mean(anemia2c$age[anemia2c$age_cat1 == 1])
mean2 <- mean(anemia2c$age[anemia2c$age_cat2 == 1])
mean3 <- mean(anemia2c$age[anemia2c$age_cat3 == 1])
mean4 <- mean(anemia2c$age[anemia2c$age_cat4 == 1])
means <- c(mean1,mean2,mean3,mean4)
plot_data <- data.frame(betas,means)
library(ggplot2)
plot01 <- ggplot(plot_data, aes(x = means, y = betas)) +
  geom_point() +
  geom_line() +
  xlab("Age") +
  ylab("Age Effect Estimate (Betas)")
plot01

# Problem 3 - End-stage renal disease (ESRD; also referred to as ‘renal failure’)
# is increasing in many countries worldwide, including the United States and Canada.
# Due to the shortfall in the available donor organs, donor kidneys are now being
# transplanted which would in the past have been discarded; the so-called Expanded 
# Criteria Donor (ECD) kidneys. By definition, ECD kidneys are more likely to suffer 
# graft failure (GF), the condition wherein the transplanted kidney stops functioning
# sufficiently. A random sample of U.S. transplant recipients was assembled, in order
# to study the effects on the mortality hazard of ECD (vs non-ECD) kidneys and graft 
# failure (GF). Data are contained in the file “kidney-ECD-1.sas7bdat”, with fields:

#   IDNUM: patient ID number
#   ECD: equals 1 for an ECD kidney, and 0 for non-ECD
#   time-to-GF: time until graft failure (missing, if GF did not occur)
#   time-to-death: time until death (missing, if death not observed)
#   time-to-censor: potential time until censoring (non-missing for all patients)
#   AGE: age at transplant (years)
#   SEX
#   DIABETES: indicator that diabetes was the cause of renal failure
#   COMORBID: number of comorbid conditions (illnesses, not counting ESRD,
#   existing at the time of transplant)
kidney <- read_sas("~/WORKING_DIRECTORIES/biostat.675/kidney_ecd_1.sas7bdat")

# (a) Fit a model which contains only factors known at the time of transplant (t = 0).
#     List the factors that significantly predict death.
kidney2 <- mutate(kidney,
                  death = 1-is.na(time_to_death),
                  time_to_event = time_to_death)
for(i in 1:nrow(kidney2)){
  if(is.na(kidney2$time_to_event[i])){
    kidney2$time_to_event[i] <- kidney2$time_to_censor[i]
  }
}
model4 <- coxph(data = kidney2,
                formula = Surv(time_to_event, death) ~ age + 
                  male + diabetes + comorbid + ECD)
summary(model4)