# function for alpha dose rate of 4-11 um

calc_dr_alpha_fine <- function() {

a2k <- 0.90 # Correction factor to convert a-value (k3.7) into keff.

#conversion factors
conversion.alpha.U.preRadon <- conversion[1,'U_alpha']
conversion.alpha.U.preRadon.err <- conversion[2,'U_alpha']
conversion.alpha.U.afterRadon <- conversion[3,'U_alpha']
conversion.alpha.U.afterRadon.err <- conversion[4,'U_alpha']
conversion.alpha.Th <- conversion[1,'Th_alpha']
conversion.alpha.Th.err <- conversion[2,'U_alpha']

#attenuation factors fitted for 1-15 um, from Brenann (1991).
grain <- (grain.min+grain.max)/2
attenuation.alphaU <- 1-0.0214*grain
attenuation.alphaU.err <- 0.0214*(grain.max - grain.min)/2
attenuation.alphaTh <- 1-0.01685*grain
attenuation.alphaTh.err <- 0.01685*(grain.max - grain.min)/2

#calculation
alpha.U <- a2k*U*attenuation.alphaU*alpha.value*(conversion.alpha.U.preRadon+(1-Radon.loss)*conversion.alpha.U.afterRadon)
alpha.U.error <- alpha.U*sqrt((U.error/U)^2+(attenuation.alphaU.err/attenuation.alphaU)^2+(alpha.value.err/alpha.value)^2+(conversion.alpha.U.preRadon.err^2+(1-Radon.loss)^2*conversion.alpha.U.afterRadon.err^2)/(conversion.alpha.U.preRadon+(1-Radon.loss)*conversion.alpha.U.afterRadon)^2)

alpha.Th <- a2k*Th*attenuation.alphaTh*alpha.value*conversion.alpha.Th
alpha.Th.error <- alpha.Th*sqrt((Th.error/Th)^2+(attenuation.alphaTh.err/attenuation.alphaTh)^2+(alpha.value.err/alpha.value)^2+(conversion.alpha.Th.err/conversion.alpha.Th)^2)

alpha <- (alpha.U+alpha.Th)/(1+1.50*water)
#alpha.error <- sqrt((1/(1+1.5*water))^2*alpha.U.error^2+(1/(1+1.5*water))^2*alpha.Th.error^2+((alpha.U+alpha.Th)/(1+1.5*water)^2)^2*1.5^2*water.error^2)
alpha.error <- alpha*sqrt((alpha.U.error^2 + alpha.Th.error^2)/(alpha.U+alpha.Th)^2 + (1.5*water.error)^2/(1+1.5*water)^2)

#cat(' Alpha attenuation factor for U=',attenuation.alphaU,'+-',attenuation.alphaU.err,"\n",
#    'Alpha attenuation factor for Th=',attenuation.alphaTh,'+-',attenuation.alphaTh.err,"\n")

c(alpha,alpha.error)
}