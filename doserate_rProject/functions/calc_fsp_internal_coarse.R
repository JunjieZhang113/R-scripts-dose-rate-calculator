# function for internal dose rate of coarse grains, suitable for 20-1000 um

source("functions/calc_beta_absorb.R")

calc_fsp_internal_coarse <- function() {
  # internal K and Rb inside K-feldspar
  K.internal <- 12.5 ## Internal K of KF is set as 12.5 +- 0.5%, from Huntley and Baril (1997) and Zhao and Li (2005).
  K.internal.err <- 0.5 ## unit percent
  Rb.internal <- 400 ## Internal Rb of KF is 400 +- 100 ppm, from Huntley and Hancock (2001).
  Rb.internal.err <- 100 ## unit ppm

  # conversion parameters for K
  conversion.beta.K <- conversion[1, "K_beta"]
  conversion.beta.K.err <- conversion[2, "K_beta"]

  # to get the absorption parameters
  beta_absorb <- calc_beta_absorb()
  absorb.K <- beta_absorb[5]
  absorb.K.err <- beta_absorb[6]
  absorb.Rb <- beta_absorb[7] # absorb.Rb is fitted from Readhead(2002),unit µGy/a/(ppm Rb). It is only for the internal dose, and it already includes the conversion factor
  absorb.Rb.err <- beta_absorb[8]

  if (etch.depth != 0) {
    # combined etching factor: the degree of increase of the beta absorption factor after etching, etch depth from 0-40 um, raw data from Brennan (2003), cited in DARC
    etch.factor <- 1.00 + 0.00954 * etch.depth - 1.618E-4 * etch.depth^2 + 1.425E-6 * etch.depth^3
    absorb.K <- absorb.K * etch.factor
  }

  # to calculate dose rate
  internal <- absorb.K * conversion.beta.K * K.internal + absorb.Rb * Rb.internal / 1000
  internal.error <- sqrt((absorb.K * conversion.beta.K * K.internal)^2 * ((absorb.K.err / absorb.K)^2 + (conversion.beta.K.err / conversion.beta.K)^2 + (K.internal.err / K.internal)^2) + (absorb.Rb * Rb.internal / 1000)^2 * ((absorb.Rb.err / absorb.Rb)^2 + (Rb.internal.err / Rb.internal)^2))

  c(internal, internal.error) # The function will return these two values
}
