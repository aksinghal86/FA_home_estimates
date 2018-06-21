##########################ATTORNEY-CLIENT PRIVILEGED############################

#' In-home formaldehyde estimates from Lumber Liquidators' laminate flooring
#' Plaintiffs: Thomas, Choe and Stein
#' asinghal 5/17/2018

#' NOTE: No intact emissions data were provided by the Plaintiff, i.e., 
#' only deconstructed board data were provided.  
#' Therefore, emissions specific to each house estimated using two approaches: 
#' 1) Using Phase I intact-board ER data
#' mean, median and 95% of ER (along with decay model in Sheehan et al. 2018) 
#' were used to estimate emissions at present time and at the time of 
#' installation.
#' 2) Using deconstructed board data
#' estimated by taking the ratio of intact board ER at time T = 0 
#' and the deconstructed board ER at time T = 0 with the assumption that 
#' the ratio would hold at any time T

library(data.table)
library(lubridate)
library(ggplot2)
library(car)
library(psych)
library(dunn.test)

setwd("/Users/asinghal/Desktop/LL/DataFiles/")

# Decay Model ------------------------------------------------------------------
#' Use the decay model from Sheehan et al. (2018) to estimate ER at time T = 0, 
#' i.e., immediately after installation 
#' Parameter values below are from the model published. 
A <- 0.7079399664092242
B <- 0.01637946082640082
C <- 0.00014477698872096014

#' Bi-exponential model -- fast decay followed by a slow gradual decay
decayModel <- function(t) {A * exp(-B * t) +  (1 - A) * exp(-C * t)}


# Data files -------------------------------------------------------------------

# Phase I intact board ER data
p1 <- fread("Phase I with ProdDate 01232017.csv")
str(p1)

# Phase II intact board ER data from boards removed from homes
p2 <- fread("LL Phase II data 2018-03-23.csv")
str(p2)

# Deconstructed data from Phase I LL testing
decon <- fread("Deconstructed.csv")
str(decon)


# Data clean-up and formatting -------------------------------------------------

# Phase I #

# Convert date into R-compatible format
p1[, ProdDate:=mdy(ProdDate)]

# Change column names for consistency between Phase I and Phase II
p1[, c("MfgDate", "Data", "ProdDate", "ER0", "ER") := list(ProdDate, 
                                                           "P1", NULL, 
                                                           ER, NULL)]

# Phase II #

# Convert date into R-compatible format
p2[, c("TestDate", "MfgDate", "OrderDate") := list(ymd(TestDate), 
                                                   ymd(MfgDate), 
                                                   ymd(OrderDate))]


head(sort(p2[, MfgDate]))
p2[MfgDate == "1014-11-01", MfgDate := ymd("2014-11-01")]

# Calculate duration of installation (NOTE: units are days not years!!!)
p2[, TimeInst := as.numeric(TestDate - OrderDate)]

# Calculate ER using BMI AER = 1/h and LR = 0.86 m2/m3 (from lab reports)
p2[, ER := P2Conc*1230 * 1/0.86]

# ER at time T = 0 -- back-calculated from ER at time T using the decay model 
p2[, ER0 := ER/decayModel(TimeInst)]

p2[, Data := "P2"]

# Merge P1 and P2 data
alldat <- rbind(p1[, .(MfgDate, ER0, Data)], p2[, .(MfgDate, ER0, Data)])

# Decon data #
decon[, ':=' (MfgDate = mdy(MfgDate), Vendor = as.character(Vendor))]


# Approach 1: using intact-board ER from Phase I data collected previously -----

#' Comparison between P1 & P2 
#' There should be no difference in ER at time T = 0 days, i.e., ER0, b/c
#' P2 ER0 is back-calculated from decay model, contingent upon P1 ER0 data

#' Only thing to check in the combined data is whether any significant diffs
#' by mfg month or by mfg year; mfg date too granular

#' First, a check to see whether the pooled data satisfy parametric assumptions
plot(density(alldat[, ER0], na.rm = TRUE))
# Looks roughly lognormally distributed
lines(density(rlnorm(1000, meanlog = 3.91, sdlog = 0.69)), col = "blue")

describe(log(alldat[, ER0]))
ks.test(p1[, ER0], "plnorm", meanlog = 3.91, sdlog = 0.69)
# Kolmogorov-Smirnov test disagrees. Not lognormally distributed (p = 0.0384)

# Kruskal-Wallis test to assess differences by manufacturing date/year
kruskal.test(ER0 ~ month(MfgDate), alldat)
# Interestingly, Mfg month significant (p = 0.01026)
with(alldat, dunn.test(ER0, factor(month(MfgDate))))
with(alldat[month(MfgDate) != 1, ], dunn.test(ER0, factor(month(MfgDate))))
#' Difference appears to be due to boards manufactured in January 
#' p = 0.05312 after removing January boards

alldat[, list(.N, mean(ER0, na.rm = TRUE)), by=month(MfgDate)]


kruskal.test(ER0 ~ factor(year(MfgDate)), alldat)
# Mfg year not significant (p = 0.6228)

alldat[, list(.N, mean(ER0, na.rm = TRUE)), by=year(MfgDate)]

# Plots #
ggplot(alldat[!is.na(MfgDate),], aes(x=factor(month(MfgDate)), y=ER0, color=Data)) + 
  geom_boxplot()+
  theme_minimal()+
  labs(x="Month  manufactured", y = "ER immediately after installation (ug/m2/h)")

ggplot(alldat[!is.na(MfgDate),], aes(x=factor(year(MfgDate)), y=ER0, color=Data)) + 
  geom_boxplot()+
  theme_minimal()+
  labs(x="Year  manufactured", y = "ER immediately after installation (ug/m2/h)")


# Same as above but only using P1 data since more appropriate and relevant
plot(density(p1[, ER0]))
# looks approximately lognormal
describe(log(p1[, ER0]))
lines(density(rlnorm(1000, meanlog = 3.86, sdlog = 0.71)), col = "blue")
ks.test(p1[, ER0], "plnorm", meanlog = 3.86, sdlog = 0.71)
# lognormally distributed (p = 0.3348)

# P1-specific comparison for mfg month and years
mod1 <- lm(log(ER0) ~ factor(month(MfgDate)), p1)
summary(mod1) 
# Appears that summer months have higher ER0s

p1[, list(.N, mean(ER0, na.rm = TRUE)), by=month(MfgDate)]
# Not due to a small sample size in any given month 
# The difference in means by month is rather large

ggplot(p1[!is.na(MfgDate),], aes(x=factor(month(MfgDate)), y=ER0)) + 
  geom_boxplot()+
  theme_minimal()+
  labs(x="Month  manufactured", y = "ER immediately after installation (ug/m2/h)")


mod2 <- lm(log(ER0) ~ factor(year(MfgDate)), p1)
summary(mod2)
# Appears that boards manufactured between 2013 and 2015 have higher ER0s

p1[, list(.N, mean(ER0, na.rm = TRUE)), by=year(MfgDate)]
# However, only 6 boards manufactured in 2011 and 2012

mod3 <- lm(log(ER0) ~ factor(year(MfgDate)), p1[!year(MfgDate) %in% c(2011, 2012), ])
summary(mod3)
#' Boards manufactured in 2015 have lower ER0s than 2013 and 2014
#' However, the difference is small ~16% (see the means)
#' Interestingly, this effect goes away when the P1 and P2 data 
#' are pooled (see above)

ggplot(p1[!is.na(MfgDate),], aes(x=factor(year(MfgDate)), y=ER0)) + 
  geom_boxplot()+
  theme_minimal()+
  labs(x="Month  manufactured", y = "ER immediately after installation (ug/m2/h)")


# Stats #
#' Mean, median and 95% ER0 from new board emission data 
#' For now, ignore differences by year and by month, because mfg month and year 
#' info are not available at the moment for the Plaintiffs 
(ER0 <- c(Mean = mean(p1$ER), quantile(p1$ER0, c(0.50, 0.95))))


# Approach 2: using Phase I ER and deconstructed data ratio --------------------

plot(density(decon[, ER]), ylim=c(0, 6e-4))
describe(log(decon[, ER]))
lines(density(decondistr), col = "blue")

shapiro.test(log(decon[, ER])) 
# approx. lognormally distributed (p = 0.1116)
ks.test(decon[, ER], "plnorm", meanlog = 7.22, sdlog = 0.58) 
# approx. lognormally distributed (p = 0.492)

#' Compare to see if any differences by vendor or by manufacturing date 
mod1 <- lm(log(ER) ~ factor(Vendor) + factor(month(MfgDate)), decon)
summary(mod1)

decon[, list(.N, mean(ER, na.rm = TRUE)), by=.(Vendor, month(MfgDate))]
# N is way too small for this comparison
# Assume no differences by either 

# Stats #
er.decon <- c(Mean = mean(decon$ER, na.rm = TRUE), 
             quantile(decon$ER, c(0.50, 0.95), na.rm = TRUE))


# Ratio of intact to decon ER estimate -----------------------------------------
# Using Monte Carlo simulation
describe(log(p1[, ER0]))
intactdistr <- rlnorm(10000, meanlog = 3.86, sdlog = 0.71)

describe(log(decon[, ER]))
decondistr <- rlnorm(10000, meanlog = 7.22, sdlog = 0.58)

ratio <-  intactdistr/decondistr
ratio.stat <- c(Mean = mean(ratio), quantile(ratio, c(0.50, 0.95)))

# In-home FA concentration estimates due specifically to LL flooring ----------- 

# Data table of ERs and different parameters------------------------------------
#' AER estimates from Ken's analysis based on Persily et al. contingent upon 
#' age of home, type of home, geographic region
AER <- c(Thomas = 0.595, Choe = 0.316, Stein = 0.287, Jensen = 0.287)

#' Fraction of flooring from Sawyer's expert report / fact sheets
Ffloor <- c(Thomas = 0.64, Choe = 0.63, Stein = 0.70, Jensen = 811/1760)

#' From Sawyer's expert report and exhibits
Tinst <- c(Thomas = round(7*365.25), 
           Choe = round(3.5*365.25), 
           Stein = round(8.33*365.25), 
           Jensen = round(3.5*365.25))

#' Create a data table of parameters
(params <- data.table(Plaintiff = names(AER), AER, Ffloor, Tinst))
params[, H := 2.59]

#' Use the decay model provided in Sheehan et al. (2018) to estimate 
#' decay factor using Tinst
params[, DecayFactor := sapply(Tinst, function (x) {
  integrate(decayModel, lower=0, upper=x)$value/x
  })]

#' Create a data table of different ERs
est <- data.table(Plaintiff = names(AER),t(ER0))
est <- melt(est, variable.name = "Stat", 
            value.name = "ER0", id.vars = "Plaintiff")

est$ERdecon <-  er.decon[match(est$Stat, names(er.decon))]
est$Ratio <- ratio.stat[match(est$Stat, names(ratio.stat))]
est <- est[order(est[, Plaintiff])]

est <- est[params, on = "Plaintiff"]


# Concentration estimates using NEW board emissions data -----------------------

# Concentration at the time of installation from P1 ER0 data
est[, C0 := ER0 * Ffloor / (AER * H)]

# ER at the time of removal, if applicable, using decay model 
est[, ERatT := ER0 * decayModel(Tinst)]

# Concentration at the time of removal, if applicable, using ERatT
est[, CatT := ERatT * Ffloor / (AER * H)]

# TWA concentration using decayfactor (see above)
est[, CTWA := C0 * DecayFactor]


# Concentration estimates using DECONSTRUCTED board emissions data -------------

#' Theoretical deconstructed ER at the time of removal, if applicable, 
#' based on P1 decon ER. Contingent on the decay model 
est[, TheorERdeconatT := ERdecon * decayModel(Tinst)]

#' From Plaintiff deconstruted results from Sawyer expert report
#' Thomas decon chamber conc = 0.074 ppm (average of 0.109 and 0.039 ppm)
#' Q/A = 1.17 m/h (not provided but assumed given that it's MDF)
est[Plaintiff == "Thomas", PlERdecon := 0.074 * 1230 * 1.17]

#' Choe decon chamber conc = 0.268 ppm (average of 8 boards tested) from Sawyer 
#' Q/A = 1.17 m/h (not provided but assumed given that it's MDF)
est[Plaintiff == "Choe", PlERdecon := 0.268 * 1230 * 1.17]

# Jensen concentration = 0.198 ppm and Q/A = 1.92 m/h
est[Plaintiff == "Jensen", PlERdecon := 0.198 * 1230 * 1.92]

#' For Stein there are no deconstructed data.  The sample collected at 
#' LA Testing is actually a bulk head space analysis, which is difficult to 
#' interpret.  So, not included in the deconstruction analysis. 

#' Estimate ER at time T from Plaintiff-provided deconstructed data and 
#' the ER0/ERdecon ratio from P1 data and assuming that the ratio remains
#' the same over time. 
est[, PlERatT := PlERdecon * Ratio]

#' Estimate PlER0 from PlERatT using the decay model from ER, i.e., ER0
#' back-calculated from estimated Plaintiff ER at T from decon data
est[, PlER0 := PlERatT / decayModel(Tinst)]

# Concentration at the time of installation
est[, PlC0 := PlER0 * Ffloor/(AER * H)]

# Concentration at the time of removal, if applicable, using ERatT
est[, PlCatT := PlERatT * Ffloor / (AER * H)]

# TWA concentration
est[, PlCTWA := PlC0 * DecayFactor]

#' In some cases, it's possible that FA concentration at time of installation  
#' from laminate flooring was greater than 100 ug/m3
#' Another way to look at it is to estimate the probablility of laminate
#' flooring contributing >= 100 ug/m3 FA 
params[, Param := Ffloor/(AER * H)]
params[, C0distr := lapply(Param, function(x) x * p1[, ER0])]
params[, P100ugperm3 := lapply(C0distr, function(x) ecdf(unlist(x))(100) * 100), by=Plaintiff]

params[, .(Plaintiff, P100ugperm3)]
##' Plaintiff P100ugperm3
##'    Choe    91.22807
##'  Jensen    95.48872
##'   Stein    86.21554
##'  Thomas    98.74687


# Write to file-----------------------------------------------------------------
write.csv(est, "Plaintiff concentration estimates AS.csv", row.names = FALSE)


################ATTORNEY-CLIENT PRIVILEGED######################################
####WARNING: DRAFT MATERIAL.  SUBJECT TO CHANGES AND REVISIONS##################
################################################################################
