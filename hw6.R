# Homework 6 --------------------------------------------------------------------------
#
# Problem 1 - End-stage renal disease (ESRD; also referred to as ‘renal failure’)
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
kidney_data <- read_sas("~/WORKING_DIRECTORIES/biostat.675/kidney_ecd_1.sas7bdat")

# FROM HOMEWORK 5.3(a)
#     Fit a model which contains only factors known at the time of transplant (t = 0).
#     List the factors that significantly predict death.
library(dplyr)
library(survival)
kidney <- mutate(kidney_data,
                 death = 1-is.na(time_to_death),
                 time_to_event = time_to_death)
for(i in 1:nrow(kidney)){
  if(is.na(kidney$time_to_event[i])){
    kidney$time_to_event[i] <- kidney$time_to_censor[i]
  }
}
model01 <- coxph(data = kidney,
                formula = Surv(time_to_event, death) ~ age + 
                  male + diabetes + comorbid + ECD)
summary(model01)

# (a) Fit a model with graft failure (GF) as a time-dependent covariate.

# put into time-dependent format
kidney$t1 <- 0
kidney$t2 <- 0
kidney$GF <- 0
kidney$death_new <- 0
kidA <- kidney
kidB <- kidney
kidC <- kidney
for(i in 1:nrow(kidney)){
  if(is.na(kidney$time_to_GF[i]) == F){
    kidA$t1[i] <- 0
    kidB$t1[i] <- kidney$time_to_GF[i]
    kidA$t2[i] <- kidney$time_to_GF[i]
    kidB$t2[i] <- kidney$time_to_event[i]
    kidA$GF[i] <- 0
    kidB$GF[i] <- 1
    kidA$death_new[i] <- 0
    kidB$death_new[i] <- kidney$death[i]
  } else{
    kidC$t1[i] <- 0
    kidC$t2[i] <- kidney$time_to_event[i]
    kidC$GF[i] <- 0
    kidC$death_new[i] <- kidney$death[i]
  }
}
kidA <- subset(kidA,is.na(kidney$time_to_GF) == F)
kidB <- subset(kidB,is.na(kidney$time_to_GF) == F)
kidC <- subset(kidC,is.na(kidney$time_to_GF) == T)
kidney_TD <- rbind(kidA,kidB,kidC)
kidney_TD <- arrange(kidney_TD, idnum)

# now fit the model
model02 <- coxph(data = kidney_TD,
                 formula = Surv(t1, t2, death_new) ~ age + 
                   male + diabetes + comorbid + ECD + GF)
summary(model02)


# Problem 2 - Asthma remains one of the most common chronic childhood illnesses, and
# a leading cause of hospital admissions. The file, asthma_1.sas7bdat contains data
# obtained from Alberta Health, a provincial health organization in Canada. Through
# linkages to several health administrative databases, a random sample of children
# born between January 1, 1995 and December 31, 1999 were retrospectively followed
# until their first physician visit for asthma or the end of the observation period
# (December 31, 1999), whichever occurred first.
library(haven)
asthma <- read_sas("~/WORKING_DIRECTORIES/biostat.675/asthma_1.sas7bdat")
# The file contains the following fields:
#   IDNUM: patient ID number
#   DT_BIRTH: date of birth
#   DT_ASTHMA: date of first physician visit reporting asthma
#   BWT: birth weight (in kg)
#   SEX
#   URBAN: indicator for living in an urban (as opposed to rural) area
#   RESP_DIST: indicator for experiencing respiratory distress at birth

# (a) Fit a model which assumes proportionality for all covariates. Code BWT as a
#     continuous covariate. Which factors appear to significantly affect asthma
#     incidence?

# (b) Repeat (a), but code BWT using an indicator for low birth weight (defined as
#     weighing <= 2.5 kg). Compare the parameter estimates with those from (a) and
#     comment on the similarities and/or differences.

# (c) Suppose, for part (c) only, that RESP_DIST was of no interest, except as an
#     adjustment covariate. Suppose also that you have no knowledge (an no desire
#     to learn) about the nature of theh non-proportionality. Fit an appropriate
#     model, and briefly defend your choice.

# (d) Fit a model which assumes that the RESP_DIST effect follows a year-specific
#     step function. Interpret the RESP_DIST effect, as estimated from this model.

# (e) Plot the age-specific RESP_DIST against the year mid-points. Describe the shape
#     of the plot and its implications (if any) for modelling the RESP_DIST effect.

# (f) Fit a model wherein the RESP_DIST regression coefficient is assumed to change
#     linearly with age (scaled to years). Interpret your parameter estimates.

# (g) Based on the model in (f), estimate the age at which children with and without
#     RESP_DIST have equal asthma hazard.
