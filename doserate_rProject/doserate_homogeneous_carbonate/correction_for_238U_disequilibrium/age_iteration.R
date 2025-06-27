# This is the R script used for speleothem age iteration due to the disequilibrium of the 238U decay chain.
# Reference: Zhang et al (2024). Isothermal thermoluminescence dating of speleothem growth - A case study from Ble?berg cave 2, Germany. Quaternary Geochronology 85, 101628.

### For sample LUM4347 from the reference paper, use the iteration method to get the age, modified from from Ikeya and Ohmura (1983).

rm(list = ls())

De <- 257.5 # unit Gys
De.err <- 15 # unit Gy

##################
U <- 0.67 # total U (U238+U235), ppm
U.err <- 0.006 # unit ppm

Th232 <- 0.0123 # Th232 ppm
Th232.err <- 0.001

##################
r0 <- 2.75 # initial [U234]/[U238] activity ratio, from U series dating
##################
Sa <- 20.8 # Measured Sa value, unit uGy/(1000*alpha*cm-2)
errSa <- 0.5 #
##################

gamma.attenuation <- 0.5 # Because the studied sample (BB2-1) is not infinitely large considering the range of gamma ray.
                         # It should be 1, when the sample is sufficiently large.

############### no need to change the code below ###############################
################################################################################
cosmicray <- 0.051 # Gy/ka, calculated cosmic ray

U238 <- 0.9929 * U # unit ppm
U235 <- 0.0071 * U # unit ppm

# U238 decay system
halflife.U238 <- 1.41 * 10^17 # unit s, From Guerin2011
halflife.U234 <- 7.75 * 10^12
halflife.Th230 <- 2.38 * 10^12

lamda.U238 <- log(2) / halflife.U238 # unit s-1
lamda.U234 <- log(2) / halflife.U234
lamda.Th230 <- log(2) / halflife.Th230

lamda.Th230 <- lamda.Th230 * 1000 * 365 * 24 * 3600 # unit ka-1
lamda.U234 <- lamda.U234 * 1000 * 365 * 24 * 3600

# conversion factor, unit Gy/ka

conversion.beta.U238 <- 0.0548
conversion.beta.U234 <- 0.0007
conversion.beta.Th230 <- 0.0868

conversion.gamma.U238 <- 0.0016
conversion.gamma.U234 <- 0.0001
conversion.gamma.Th230 <- 0.1079

conversion.beta.U235 <- 0.0037 ## 1 ppm to Gy/ka, from Guerin2011
conversion.gamma.U235 <- 0.0020

conversion.beta.Th232 <- 0.0277
conversion.gamma.Th232 <- 0.0479

# alpha flux for 1ppm natural U, when in equilibrium, Norbert Mercier table #######
flux.U238 <- 1529 ## unit /cm2/year
flux.U234 <- 1829
flux.Th230 <- 14290
totalflux.U238 <- 17647

totalflux.U235 <- 821
totalflux.Th232 <- 5166

# dose rates from U235 and Th232. constant
Dr.U235 <- U * (totalflux.U235 * Sa * 0.92 * 10^-6 + conversion.beta.U235 + conversion.gamma.U235 * gamma.attenuation)
Dr.Th232 <- Th232 * (totalflux.Th232 * Sa * 0.96 * 10^-6 + conversion.beta.Th232 + Th232 * conversion.gamma.Th232 * gamma.attenuation)
# note that 0.92 and 0.96 are correction factors for U and Th.

Dr.U238 <- U * (flux.U238 * Sa * 0.92 * 10^-6 + conversion.beta.U238 + conversion.gamma.U238 * gamma.attenuation)
Dr.U234 <- U * (flux.U234 * Sa * 0.92 * 10^-6 + conversion.beta.U234 + conversion.gamma.U234 * gamma.attenuation)
Dr.Th230 <- U * (flux.Th230 * Sa * 0.92 * 10^-6 + conversion.beta.Th230 + conversion.gamma.Th230 * gamma.attenuation)
Dr.U238.total <- Dr.U238 + Dr.U234 + Dr.Th230

Dr <- Dr.Th232 + Dr.U235 + Dr.U238.total + cosmicray ## total dose rate

Dr.alpha <- U * (totalflux.U235 + totalflux.U238) * Sa * 0.92 * 10^-6 + Th232 * totalflux.Th232 * Sa * 0.96 * 10^-6
Dr.beta <- U * 0.1457 + Th232 * 0.0277
Dr.gamma <- (U * 0.1116 + Th232 * 0.0479) * gamma.attenuation

age <- De / Dr # apparent age

for (i in 0:7) # 7 cycles is already enough for it to converge
{
  if (i == 0) {
    cat("\n********************************************************\n")
    cat("Equivalent dose (De)=", De, "Gy", "\n")
    cat("Dose rate assuming equilibrium=", Dr, "Gy/ka", "\n")
    cat("Alpha Dose rate assuming equilibrium=", Dr.alpha, "Gy/ka", "\n")
    cat("Beta Dose rate assuming equilibrium=", Dr.beta, "Gy/ka", "\n")
    cat("Gamma rate assuming equilibrium=", Dr.gamma, "Gy/ka", "\n")
    cat("********************************************************\n")
    cat("iteration cycle=", i, " ")
    cat("De_measure=", De)
    cat(" age=", round(age, 2), "\n")
  }

  cat("iteration cycle=", i + 1, " ")
  Dr * age
  De_cal <- Dr * age + Dr.U234 * (r0 - 1) * (1 - exp(-lamda.U234 * age)) / lamda.U234 - Dr.Th230 * ((1 - exp(-lamda.Th230 * age)) / lamda.Th230 - 1 / lamda.U234 * (r0 - 1) * (1 - (lamda.Th230 * exp(-lamda.U234 * age) - lamda.U234 * exp(-lamda.Th230 * age)) / (lamda.Th230 - lamda.U234)))
  age <- age * De / De_cal
  cat("De_cal=", round(De_cal, 2), "   ")
  cat("age=", round(age, 2), "\n")

  if (i == 7) {
    cat("********************************************************\n")
  }
}
