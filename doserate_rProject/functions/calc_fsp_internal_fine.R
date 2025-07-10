# function for internal dose rate of fine grains, suitable for 1-20 um

calc_fsp_internal_fine <- function() {
  # internal K and Rb inside K-feldspar, can be modified
  K.internal <- 12.5 ## Internal K of KF is set as 12.5 +- 0.5%, from Huntley and Baril (1997) and Zhao and Li (2005).
  K.internal.err <- 0.5 ## unit percent
  Rb.internal <- 400 ## Internal Rb of KF is 400 +- 100 ppm, from Huntley and Hancock (2001).
  Rb.internal.err <- 100 ## unit ppm

  # conversion parameters for K
  conversion.beta.K <- conversion[1, "K_beta"]
  conversion.beta.K.err <- conversion[2, "K_beta"]


  # beta absorption factor in 1-20 um, fitted by origin2020, with data from Brennan(2003)
  grain <- (grain.min + grain.max) / 2
  absorb.K <- 0.000364 * grain
  absorb.K.err <- 0.000364 * (grain.max - grain.min) / 2

  # absorbed dose of 87Rb in 1-20 um, fit data from Readhead(2002). It include the conversion factor. The unit is uGy/a/(ppm Rb)
  absorb.Rb <- 0.1782 * (1 - exp(- grain / 56.83))
  absorb.Rb.err <- (0.1782 * (1 - exp(- grain.max / 56.83)) - 0.1782 * (1 - exp(- grain.min / 56.83))) / 2

  # to calculate internal dose rate
  internal <- absorb.K * conversion.beta.K * K.internal + absorb.Rb * Rb.internal / 1000
  internal.error <- sqrt((absorb.K * conversion.beta.K * K.internal)^2 * ((absorb.K.err / absorb.K)^2 + (conversion.beta.K.err / conversion.beta.K)^2 + (K.internal.err / K.internal)^2) + (absorb.Rb * Rb.internal / 1000)^2 * ((absorb.Rb.err / absorb.Rb)^2 + (Rb.internal.err / Rb.internal)^2))

  c(internal, internal.error) # The function will return these two values
}
