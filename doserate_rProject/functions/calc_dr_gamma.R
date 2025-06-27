# function for gamma dose rate

calc_dr_gamma <- function() {
  # gamma conversion factors
  conversion.gamma.U.preRadon <- conversion[1, "U_gamma"]
  conversion.gamma.U.preRadon.err <- conversion[2, "U_gamma"]
  conversion.gamma.U.afterRadon <- conversion[3, "U_gamma"]
  conversion.gamma.U.afterRadon.err <- conversion[4, "U_gamma"]
  conversion.gamma.Th <- conversion[1, "Th_gamma"]
  conversion.gamma.Th.err <- conversion[2, "Th_gamma"]
  conversion.gamma.K <- conversion[1, "K_gamma"]
  conversion.gamma.K.err <- conversion[2, "K_gamma"]

  # calculation
  gamma.U <- U * (conversion.gamma.U.preRadon + (1 - Radon.loss) * conversion.gamma.U.afterRadon)
  gamma.U.error <- gamma.U * sqrt((U.error / U)^2 + (conversion.gamma.U.preRadon.err^2 + (1 - Radon.loss)^2 * conversion.gamma.U.afterRadon.err^2) / (conversion.gamma.U.preRadon + (1 - Radon.loss) * conversion.gamma.U.afterRadon)^2)

  gamma.Th <- Th * conversion.gamma.Th
  gamma.Th.error <- gamma.Th * sqrt((Th.error / Th)^2 + (conversion.gamma.Th.err / conversion.gamma.Th)^2)

  gamma.K <- K * conversion.gamma.K
  gamma.K.error <- gamma.K * sqrt((K.error / K)^2 + (conversion.gamma.K.err / conversion.gamma.K)^2)

  gamma <- (gamma.K + gamma.U + gamma.Th) / (1 + 1.14 * water)
  gamma.error <- sqrt((1 / (1 + 1.14 * water))^2 * gamma.K.error^2 + (1 / (1 + 1.14 * water))^2 * gamma.U.error^2 + (1 / (1 + 1.14 * water))^2 * gamma.Th.error^2 + ((gamma.K + gamma.U + gamma.Th) / (1 + 1.14 * water)^2)^2 * 1.14^2 * water.error^2)

  c(gamma, gamma.error)
}
