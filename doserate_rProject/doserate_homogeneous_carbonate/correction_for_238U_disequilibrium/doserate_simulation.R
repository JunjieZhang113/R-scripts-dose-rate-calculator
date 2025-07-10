# This is the R script used for speleothem dose rate modelling due to the disequilibrium of the 238U decay chain.
# reference: Zhang et al (2024). Isothermal thermoluminescence dating of speleothem growth - A case study from Ble?berg cave 2, Germany. Quaternary Geochronology 85, 101628.

### For sample LUM4347 from the Zhang et al. (2024)
rm(list = ls())

output.file.name <- "doserate_homogeneous_carbonate/correction_for_238U_disequilibrium/doserate_simulated_238U_disequilibrium.csv" # define the output file name

###############  parameters below to change #####################################
################################################################################

gamma.attenuation <- 0.5 # Because the studied sample (BB2-1) is not infinitely large considering the range of gamma ray (~30 cm).
                         # It should be 1, when the sample is sufficiently large.

cosmicray <- 0.051 # Gy/ka, calculated cosmic ray for Blessberg cave 2.
##################
U <- 0.67 # total U (U238+U235), ppm
Th232 <- 0.0123 # Th232, ppm
##################
r0 <- 2.75 # initial [U234]/[U238] activity ratio, from U series dating
##################
Sa <- 20.8 # Measured Sa value, unit uGy/(1000*alpha*cm-2)
errSa <- 0.5 #



################################################################################
############### no need to change the code below ###############################
################################################################################

U238 <- 0.9929 * U # unit ppm
U235 <- 0.0071 * U # unit ppm

# U238 decay system
halflife.U238 <- 1.41 * 10^17 # unit s, from Guerin et al. (2011)
halflife.U234 <- 7.75 * 10^12
halflife.Th230 <- 2.38 * 10^12
halflife.Ra226 <- 5.05 * 10^10

lamda.U238 <- log(2) / halflife.U238 # unit s-1
lamda.U234 <- log(2) / halflife.U234
lamda.Th230 <- log(2) / halflife.Th230
lamda.Ra226 <- log(2) / halflife.Ra226

####################
a.U238 <- 12.44 * U238 #   The activity of 1 ppm of U238 is 12.44 Bg/kg, from Guerin et al. (2011)

a.U234.initial <- r0 * a.U238
N.U234.initial <- a.U234.initial / lamda.U234
N.Th230.initial <- 0
#N.Th230.initial<-(Th232*10^-3/232)*6.02*10^23*(4.4*10^-6)  ## atom ratio of 230Th/232Th is 4.4+-2.2*10^-6

# alpha flux for 1 ppm natural U, when in equilibrium,from Norbert Mercier's table
flux.U238 <- 1529 ## unit /cm2/year
flux.U234 <- 1829
flux.Th230 <- 14290
# totalflux.U238<-17647

totalflux.U235 <- 821
totalflux.Th232 <- 5166

conversion.beta.U235 <- 0.0037 ## 1 ppm to Gy/ka, from Guerin et al. (2011).
conversion.gamma.U235 <- 0.0020

conversion.beta.Th232 <- 0.0277
conversion.gamma.Th232 <- 0.0479

# dose rates from U235 and Th232. constant
Dr.U235 <- U * (totalflux.U235 * Sa * 0.92 * 10^-6) + U * conversion.beta.U235 + U * conversion.gamma.U235 * gamma.attenuation
Dr.Th232 <- Th232 * (totalflux.Th232 * Sa * 0.96 * 10^-6) + Th232 * conversion.beta.Th232 + Th232 * conversion.gamma.Th232 * gamma.attenuation
# note that 0.92 and 0.96 are correction factors for U and Th.

##### Below is the calculation of dose rate of U238 decay chain
## conversion from from 1 ppm to Gy/ka, Guerin2011  ###############
conversion.beta.U <- 0.1419 ## total chain
conversion.beta.U238 <- 0.0548
conversion.beta.U234 <- 0.0007
conversion.beta.Th230 <- 0.0864
conversion.gamma.U <- 0.1096 ## total chain
conversion.gamma.U238 <- 0.0017
conversion.gamma.U234 <- 0.0001
conversion.gamma.Th230 <- 0.1079

De.initial <- 0 # Gy/ka
time.sum <- integer(0)

# activity
a.U238.sum <- numeric(0)
a.U234.sum <- numeric(0)
a.Th230.sum <- numeric(0)
alphaflux.sum <- numeric(0)

# dose rate of 238U to 234U
Dr238.alpha.sum <- numeric(0)
Dr238.beta.sum <- numeric(0)
Dr238.gamma.sum <- numeric(0)
Dr238.sum <- numeric(0)

# dose rate of 234U to 230Th
Dr234.alpha.sum <- numeric(0)
Dr234.beta.sum <- numeric(0)
Dr234.gamma.sum <- numeric(0)
Dr234.sum <- numeric(0)

# dose rate of 230Th to end
Dr230.alpha.sum <- numeric(0)
Dr230.beta.sum <- numeric(0)
Dr230.gamma.sum <- numeric(0)
Dr230.sum <- numeric(0)

Dr.totalalpha.sum <- numeric(0) # alpha dose rate of U238 decay chain
Dr.totalbeta.sum <- numeric(0)
Dr.totalgamma.sum <- numeric(0)
Dr.total.sum <- numeric(0) # total dose rate of U238

Dr.U235.sum <- numeric(0) # total dose rate of U235
Dr.Th232.sum <- numeric(0) # total dose rate of Th232

Dr.cosmic.sum <- numeric()

Dr.sum <- numeric(0) # overall total dose rate
De.sum <- numeric(0)

for (i in 0:3000000) # unit year
{
  if (i == 0) {
    N.U234 <- N.U234.initial
    N.Th230 <- N.Th230.initial
    De <- De.initial
  }

  a.U234 <- lamda.U234 * N.U234
  a.Th230 <- lamda.Th230 * N.Th230

  ############################## calculate alpha dose rate Gy/ka ####
  ra.U234 <- a.U234 / a.U238 # relative activity of U234 to U238  (or the U234 in equilibrium state)
  ra.Th230 <- a.Th230 / a.U238 # relative activity of Th230 to U238 (or the Th230 in equilibrium state)

  totalflux <- flux.U238 + flux.U234 * ra.U234 + flux.Th230 * ra.Th230
  Dr.totalalpha <- U * totalflux * Sa * 0.92 * 10^-6 ## correction factor of 0.92 from Norbert

  Dr.alpha.U238 <- U * flux.U238 * Sa * 0.92 * 10^-6
  Dr.alpha.U234 <- U * flux.U234 * ra.U234 * Sa * 0.92 * 10^-6
  Dr.alpha.Th230 <- U * flux.Th230 * ra.Th230 * Sa * 0.92 * 10^-6

  ########################## calculate beta dose rate unit Gy/ka ##########################

  Dr.beta.U238 <- U * conversion.beta.U238
  Dr.beta.U234 <- U * conversion.beta.U234 * ra.U234
  Dr.beta.Th230 <- U * conversion.beta.Th230 * ra.Th230

  Dr.totalbeta <- Dr.beta.U238 + Dr.beta.U234 + Dr.beta.Th230

  ########################## calculate gamma dose rate unit Gy/ka ##########################

  ## here we assume only 20% of the total gamma dose rate, because the sample is not infinite large for the gamma ray
  Dr.gamma.U238 <- U * conversion.gamma.U238 * gamma.attenuation
  Dr.gamma.U234 <- U * conversion.gamma.U234 * ra.U234 * gamma.attenuation
  Dr.gamma.Th230 <- U * conversion.gamma.Th230 * ra.Th230 * gamma.attenuation

  Dr.totalgamma <- Dr.gamma.U238 + Dr.gamma.U234 + Dr.gamma.Th230

  ######################### total dose rate of U238 ######################################################

  Dr.total <- Dr.totalalpha + Dr.totalbeta + Dr.totalgamma # means total dose rate of U238

  Dr.U238 <- Dr.alpha.U238 + Dr.beta.U238 + Dr.gamma.U238 # total dose rate of the U238- to U234 segment
  Dr.U234 <- Dr.alpha.U234 + Dr.beta.U234 + Dr.gamma.U234 # total dose rate of the U234- to Th230 segment
  Dr.Th230 <- Dr.alpha.Th230 + Dr.beta.Th230 + Dr.gamma.Th230 # total dose rate of the Th230- to end of decay chain

  ####################### overall total dose rate  #########################################################

  Dr <- Dr.total + Dr.U235 + Dr.Th232 + cosmicray

  #################################################

  if (i %% 1000 == 0) {
    j <- i / 1000 + 1
    time.sum[j] <- i / 1000
    a.U238.sum[j] <- a.U238
    a.U234.sum[j] <- a.U234
    a.Th230.sum[j] <- a.Th230
    alphaflux.sum[j] <- totalflux


    Dr238.alpha.sum[j] <- Dr.alpha.U238
    Dr238.beta.sum[j] <- Dr.beta.U238
    Dr238.gamma.sum[j] <- Dr.gamma.U238
    Dr238.sum[j] <- Dr.U238

    Dr234.alpha.sum[j] <- Dr.alpha.U234
    Dr234.beta.sum[j] <- Dr.beta.U234
    Dr234.gamma.sum[j] <- Dr.gamma.U234
    Dr234.sum[j] <- Dr.U234

    Dr230.alpha.sum[j] <- Dr.alpha.Th230
    Dr230.beta.sum[j] <- Dr.beta.Th230
    Dr230.gamma.sum[j] <- Dr.gamma.Th230
    Dr230.sum[j] <- Dr.Th230

    Dr.totalalpha.sum[j] <- Dr.totalalpha
    Dr.totalbeta.sum[j] <- Dr.totalbeta
    Dr.totalgamma.sum[j] <- Dr.totalgamma
    Dr.total.sum[j] <- Dr.total # U238 total dose rate

    Dr.U235.sum[j] <- Dr.U235
    Dr.Th232.sum[j] <- Dr.Th232

    Dr.cosmic.sum[j] <- cosmicray

    Dr.sum[j] <- Dr
    De.sum[j] <- De
  }

  De.increase <- Dr / 1000 ## De increase Gy per year
  De <- De + De.increase

  N.U234 <- N.U234 + (a.U238 - a.U234) * 3600 * 24 * 365 # per year
  N.Th230 <- N.Th230 + (a.U234 - a.Th230) * 3600 * 24 * 365 # per year
}

dataset <- data.frame(
  "Time_ka" = time.sum, "alphaflux" = alphaflux.sum, "a_U238" = a.U238.sum, "a_U234" = a.U234.sum, "a_Th230" = a.Th230.sum,
  "r234_238" = a.U234.sum / a.U238.sum, "r230_234" = a.Th230.sum / a.U234.sum,
  "Dr238_alpha" = Dr238.alpha.sum, "Dr238_beta" = Dr238.beta.sum, "Dr238_gamma" = Dr238.gamma.sum, "Dr238" = Dr238.sum,
  "Dr234_alpha" = Dr234.alpha.sum, "Dr234_beta" = Dr234.beta.sum, "Dr234_gamma" = Dr234.gamma.sum, "Dr234" = Dr234.sum,
  "Dr230_alpha" = Dr230.alpha.sum, "Dr230_beta" = Dr230.beta.sum, "Dr230_gamma" = Dr230.gamma.sum, "Dr230" = Dr230.sum,
  "Dr_U238chain_alpha" = Dr.totalalpha.sum, "Dr_U238chain_beta" = Dr.totalbeta.sum, "Dr_U238chain_gamma" = Dr.totalgamma.sum,
  "Dr_U238chain_total" = Dr.total.sum, "Dr_U235" = Dr.U235.sum, "Dr_Th232" = Dr.Th232.sum, "Dr.cosmic" = Dr.cosmic.sum,
  "Dr_total" = Dr.sum, "De_Gy" = De.sum
)
head(dataset)

write.csv(dataset, output.file.name, row.names = FALSE)

## Graph ###
par(mfrow = c(2, 3))
par(mar = c(5, 5, 2, 2))
plot(dataset$Time_ka, dataset$a_U238, xlab = "Time (ka)", ylab = "Activity (Bq/kg)", type = "n", xlim = c(0, 2500), ylim = c(0, 30))
lines(dataset$Time_ka, dataset$a_U238, col = "black")
lines(dataset$Time_ka, dataset$a_U234, col = "red")
lines(dataset$Time_ka, dataset$a_Th230, col = "blue")
legend("topright",
  legend = c("[238U]", "[234U]", "[230Th]"),
  col = c("black", "red", "blue"), lty = 1, lwd = 1, cex = 1
)

plot(dataset$Time_ka, dataset$a_U238, xlab = "Time (ka)", ylab = "Dr (Gy/ka)", type = "n", xlim = c(0, 2500), ylim = c(0, 0.80))
lines(dataset$Time_ka, dataset$Dr238, col = "black")
lines(dataset$Time_ka, dataset$Dr234, col = "red")
lines(dataset$Time_ka, dataset$Dr230, col = "blue")
lines(dataset$Time_ka, dataset$Dr_U238chain_total, col = "cyan3", lwd = 1)
legend("topright",
  legend = c("238U-234U", "234U-230Th", "230Th-206Pb", "238U full chain"),
  col = c("black", "red", "blue", "cyan3"), lty = 1, lwd = 1, cex = 1
)

plot(dataset$Time_ka, dataset$De, xlab = "Time (ka)", ylab = "De (Gy)", type = "n", xlim = c(0, 2500), ylim = c(0, 1500))
lines(dataset$Time_ka, dataset$De, col = "brown")
points(x = 0, y = 258, col = "blue")
arrows(0, 258 - 12, 0, 258 + 12, lwd = 1, length = 0.05, angle = 90, code = 3, col = "blue")
lines(x = 1:421, y = rep(258, 421), lty = 2, col = "blue")
lines(y = -20:258, x = rep(421, 258 + 21), lty = 2, col = "blue")

################################################

plot(dataset$Time_ka, dataset$a_U238, xlab = "Time (ka)", ylab = "Activity (Bq/kg)", type = "n", xlim = c(0, 500), ylim = c(0, 30))
lines(dataset$Time_ka, dataset$a_U238, col = "black")
lines(dataset$Time_ka, dataset$a_U234, col = "red")
lines(dataset$Time_ka, dataset$a_Th230, col = "blue")
legend("topright",
  legend = c("[238U]", "[234U]", "[230Th]"),
  col = c("black", "red", "blue"), lty = 1, lwd = 1, cex = 1
)

plot(dataset$Time_ka, dataset$a_U238, xlab = "Time (ka)", ylab = "Dr (Gy/ka)", type = "n", xlim = c(0, 500), ylim = c(0, 0.8))
lines(dataset$Time_ka, dataset$Dr238, col = "black")
lines(dataset$Time_ka, dataset$Dr234, col = "red")
lines(dataset$Time_ka, dataset$Dr230, col = "blue")
lines(dataset$Time_ka, dataset$Dr_U238chain_total, col = "cyan3", lwd = 1)
legend("topright",
  legend = c("238U-234U", "234U-230Th", "230Th-206Pb", "238U full chain"),
  col = c("black", "red", "blue", "cyan3"), lty = 1, lwd = 1, cex = 1
)


plot(dataset$Time_ka, dataset$De, xlab = "Time (ka)", ylab = "De (Gy)", type = "n", xlim = c(0, 500), ylim = c(0, 400))
lines(dataset$Time_ka, dataset$De, col = "brown")
points(x = 0, y = 258, col = "blue")
arrows(0, 258 - 12, 0, 258 + 12, lwd = 1, length = 0.05, angle = 90, code = 3, col = "blue")
lines(x = 1:421, y = rep(258, 421), lty = 2, col = "blue")
lines(y = -20:258, x = rep(421, 258 + 21), lty = 2, col = "blue")
