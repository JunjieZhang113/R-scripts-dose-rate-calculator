# function for beta dose rate, suitable for 20-1000 um

source('functions/calc_beta_absorb.R')

calc_dr_beta_coarse <- function() {
  
  #conversion factors of beta
  conversion.beta.U.preRadon <- conversion[1,'U_beta']
  conversion.beta.U.preRadon.err <- conversion[2,'U_beta']
  conversion.beta.U.afterRadon <- conversion[3,'U_beta']
  conversion.beta.U.afterRadon.err <- conversion[4,'U_beta']
  conversion.beta.Th <- conversion[1,'Th_beta']
  conversion.beta.Th.err <- conversion[2,'Th_beta']
  conversion.beta.K <- conversion[1,'K_beta']
  conversion.beta.K.err <- conversion[2,'K_beta']
  
  #absorption factor
  beta_absorb <- calc_beta_absorb()
  absorb.U <- beta_absorb[1]
  absorb.U.err <- beta_absorb[2]
  absorb.Th <- beta_absorb[3]
  absorb.Th.err <- beta_absorb[4]
  absorb.K <- beta_absorb[5]
  absorb.K.err <- beta_absorb[6]
  
  if (etch.depth != 0) {
  #combined etching factor: the degree of increase of the beta absorption factor after etching, etch depth (decrease in diameter, not radius) from 0-40 um, raw data from Brennan (2003), in supplementary of DARC
  etch.factor <- 1.00+0.00954*etch.depth-1.618E-4*etch.depth^2+1.425E-6*etch.depth^3
  absorb.U <- absorb.U*etch.factor
  absorb.Th <- absorb.Th*etch.factor
  absorb.K <- absorb.K*etch.factor
  }

beta.U <- U*(1-absorb.U)*(conversion.beta.U.preRadon+(1-Radon.loss)*conversion.beta.U.afterRadon);
beta.U.error <- beta.U*sqrt((U.error/U)^2+absorb.U.err^2/(1-absorb.U)^2+(conversion.beta.U.preRadon.err^2+(1-Radon.loss)^2*conversion.beta.U.afterRadon.err^2)/(conversion.beta.U.preRadon+(1-Radon.loss)*conversion.beta.U.afterRadon)^2)

beta.Th <- Th*(1-absorb.Th)*conversion.beta.Th
beta.Th.error <- beta.Th*sqrt((Th.error/Th)^2+absorb.Th.err^2/(1-absorb.Th)^2+(conversion.beta.Th.err/conversion.beta.Th)^2)

beta.K <- K*(1-absorb.K)*conversion.beta.K;
beta.K.error <- beta.K*sqrt((K.error/K)^2+absorb.K.err^2/(1-absorb.K)^2+(conversion.beta.K.err/conversion.beta.K)^2)

beta <- (beta.K+beta.U+beta.Th)/(1+1.25*water);
beta.error <- sqrt((1/(1+1.25*water))^2*beta.K.error^2+(1/(1+1.25*water))^2*beta.U.error^2+(1/(1+1.25*water))^2*beta.Th.error^2+((beta.K+beta.U+beta.Th)/(1+1.25*water)^2)^2*1.25^2*water.error^2)

#cat(" absorb.K=",absorb.K,'+-',absorb.K.err,"\n",
#    "absorb.U=",absorb.U,'+-',absorb.U.err,"\n",
#    "absorb.Th=",absorb.Th, '+-',absorb.Th.err,"\n")

c(beta, beta.error)

}