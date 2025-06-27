# function for beta dose rate for fine grains, between 1 um and 20 um

calc_dr_beta_fine <- function() {
  # conversion factors of beta
  conversion.beta.U.preRadon <- conversion[1, "U_beta"]
  conversion.beta.U.preRadon.err <- conversion[2, "U_beta"]
  conversion.beta.U.afterRadon <- conversion[3, "U_beta"]
  conversion.beta.U.afterRadon.err <- conversion[4, "U_beta"]
  conversion.beta.Th <- conversion[1, "Th_beta"]
  conversion.beta.Th.err <- conversion[2, "Th_beta"]
  conversion.beta.K <- conversion[1, "K_beta"]
  conversion.beta.K.err <- conversion[2, "K_beta"]

  # beta absorption factor in 1-20 um, fitted by origin2020, with data from Brennan(2003)
  grain <- (grain.min + grain.max) / 2
  absorb.K <- 0.000364 * grain
  absorb.K.err <- 0.000364 * (grain.max - grain.min) / 2
  absorb.U <- 0.06812 * (1 - exp(- grain / 27.81))
  absorb.U.err <- (0.06812 * (1 - exp(- grain.max / 27.81)) - 0.06812 * (1 - exp(- grain.min / 27.81))) / 2
  absorb.Th <- 0.1017 * (1 - exp(- grain / 33.39))
  absorb.Th.err <- (0.1017 * (1 - exp(- grain.max / 33.39)) - 0.1017 * (1 - exp(- grain.min / 33.39))) / 2
  
  # calculate beta dose rate
  beta.U <- U * (1 - absorb.U) * (conversion.beta.U.preRadon + (1 - Radon.loss) * conversion.beta.U.afterRadon)
  beta.U.error <- beta.U * sqrt((U.error / U)^2 + absorb.U.err^2 / (1 - absorb.U)^2 + (conversion.beta.U.preRadon.err^2 + (1 - Radon.loss)^2 * conversion.beta.U.afterRadon.err^2) / (conversion.beta.U.preRadon + (1 - Radon.loss) * conversion.beta.U.afterRadon)^2)

  beta.Th <- Th * (1 - absorb.Th) * conversion.beta.Th
  beta.Th.error <- beta.Th * sqrt((Th.error / Th)^2 + absorb.Th.err^2 / (1 - absorb.Th)^2 + (conversion.beta.Th.err / conversion.beta.Th)^2)

  beta.K <- K * (1 - absorb.K) * conversion.beta.K
  beta.K.error <- beta.K * sqrt((K.error / K)^2 + absorb.K.err^2 / (1 - absorb.K)^2 + (conversion.beta.K.err / conversion.beta.K)^2)

  beta <- (beta.K + beta.U + beta.Th) / (1 + 1.25 * water)
  beta.error <- sqrt((1 / (1 + 1.25 * water))^2 * beta.K.error^2 + (1 / (1 + 1.25 * water))^2 * beta.U.error^2 + (1 / (1 + 1.25 * water))^2 * beta.Th.error^2 + ((beta.K + beta.U + beta.Th) / (1 + 1.25 * water)^2)^2 * 1.25^2 * water.error^2)

  c(beta, beta.error)
}
