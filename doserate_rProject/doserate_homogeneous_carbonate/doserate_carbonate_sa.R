## Dose rate calculation for carbonates, with the Sa system, assuming an infinite homogeneous medium

rm(list = ls())

Dr <- read.csv("doserate_homogeneous_carbonate/Template_carbonate_sa.csv", header = TRUE)
# head(Dr)
output_filename <- "doserate_homogeneous_carbonate/Doserate_results_carbonate_sa.csv"

Radon.loss <- 0.0 # for radon loss. 0.2 means 20% of the radon is escaped.
conversion <- read.csv("conversion_data/Liritzis_2013_conversion.csv", header = TRUE)


#################################################################################
### The code below are for the calculations. Select all (ctrl+A), and run.#######
#################################################################################

# load the functions
source("functions/calc_cosmicray.R")
# source('functions/calc_dr_alpha_infinite.R') #not needed for the sa-system
source("functions/calc_dr_beta_infinite.R")
source("functions/calc_dr_gamma.R")

# calculation
for (i in 1:nrow(Dr)) {
  # mineral
  mineral <- Dr[i, "mineral"]

  # alpha efficiency
  Sa <- Dr[i, "Sa"] # unit
  Sa.err <- Dr[i, "Sa_err"]

  # geography information
  geo.lat <- Dr[i, "latitude"] # latitude
  geo.lon <- Dr[i, "longitude"] # longitude
  altitude <- Dr[i, "altitude_m"] # altitude，unit m
  depth <- Dr[i, "depth_m"] # the sampling depth from the surface，unit m

  # water
  water <- Dr[i, "water_percent"] / 100
  water.error <- Dr[i, "water_err"] / 100

  # U, Th, K concentrations
  U <- Dr[i, "U_ppm"]
  U.error <- Dr[i, "U_err"]

  Th <- Dr[i, "Th_ppm"]
  Th.error <- Dr[i, "Th_err"]

  K <- Dr[i, "K_percent"]
  K.error <- Dr[i, "K_err"]


  if (!(mineral %in% c("calcite", "dolomite"))) {
    stop("The input of 'mineral' should be either 'calcite' or 'dolomite' in the template file.")
  }

  ##########################################################################
  ################################## calculations ##########################

  ############## cosmic ray dose rate ###########
  cosmic <- calc_cosmicray()
  Dr[i, "cosmicray"] <- cosmic[1]
  Dr[i, "cosmicray.err"] <- cosmic[2]

  ################### alpha dose rate ###########
  # alpha flux is calculated from ranges of alpha particles, with an excel table Norbert Mercier
  if (mineral == "calcite") {
    flux.alpha.U.preRadon <- 7793 # for calcite, alpha flux of 1ppm U is 18468 alpha particles (cm-2*year)
    flux.alpha.U.afterRadon <- 10675
    flux.alpha.Th <- 5166
  }
  if (mineral == "dolomite") {
    flux.alpha.U.preRadon <- 7585 # for dolomite, alpha flux of 1ppm U is 18013 alpha particles (cm-2*year)
    flux.alpha.U.afterRadon <- 10428
    flux.alpha.Th <- 5047
  }

  # Note that for Sa value measured in Bordeaux, a correction factor of 0.92 is used for U and 0.96 for Th.
  # Correction factors are from Norbert Mercier
  # Because the alpha source Bordeaux used has a mean energy of 3.3 MeV, which is different from the natural alpha energy spectrum of U and Th
  alpha.U <- 0.92 * U * Sa * (flux.alpha.U.preRadon + (1 - Radon.loss) * flux.alpha.U.afterRadon) * 10^-6 # unit Gy/ka
  alpha.U.error <- alpha.U * sqrt((U.error / U)^2 + (Sa.err / Sa)^2)
  alpha.Th <- 0.96 * Th * Sa * flux.alpha.Th * 10^-6
  alpha.Th.error <- alpha.Th * sqrt((Th.error / Th)^2 + (Sa.err / Sa)^2)
  alpha <- (alpha.U + alpha.Th) / (1 + 1.50 * water)
  alpha.err <- sqrt((1 / (1 + 1.5 * water))^2 * alpha.U.error^2 + (1 / (1 + 1.5 * water))^2 * alpha.Th.error^2 + ((alpha.U + alpha.Th) / (1 + 1.5 * water)^2)^2 * 1.5^2 * water.error^2)

  Dr[i, "alpha"] <- alpha
  Dr[i, "alpha_err"] <- alpha.err

  #################### beta dose rate  ##############
  beta <- calc_dr_beta_infinite()
  Dr[i, "beta"] <- beta[1]
  Dr[i, "beta_err"] <- beta[2]

  #################### gamma dose rate ##############
  gamma <- calc_dr_gamma()
  Dr[i, "gamma"] <- gamma[1]
  Dr[i, "gamma_err"] <- gamma[2]

  ##################  total dose rate ###############
  doserate <- cosmic[1] + alpha + beta[1] + gamma[1]
  doserate.error <- sqrt(cosmic[2]^2 + alpha.err^2 + beta[2]^2 + gamma[2]^2)
  Dr[i, "doserate_Gyperka"] <- doserate
  Dr[i, "doserate_err"] <- doserate.error

  if (i == nrow(Dr)) {
    write.csv(Dr, output_filename, row.names = FALSE)
    cat(
      "\n####################### Well done! ##########################\n",
      "Dose rate calculation of", nrow(Dr), "samples is finished.\n",
      "The results are saved in the csv file:\n", output_filename, "\n",
      "The dose rates of the first five samples are shown below:\n\n "
    )
    print(Dr[1:5, c("Sample", "mineral", "Sa", "doserate_Gyperka", "doserate_err")])
  }
}
