# function for alpha dose rate of infinite homogeneous medium, for carbonate

calc_dr_alpha_infinite <- function() {

#there is no a2k correction factor, as the keff is entered in the Template_carbonate_keff csv template 
  
#conversion factors
conversion.alpha.U.preRadon <- conversion[1,'U_alpha']
conversion.alpha.U.preRadon.err <- conversion[2,'U_alpha']
conversion.alpha.U.afterRadon <- conversion[3,'U_alpha']
conversion.alpha.U.afterRadon.err <- conversion[4,'U_alpha']
conversion.alpha.Th <- conversion[1,'Th_alpha']
conversion.alpha.Th.err <- conversion[2,'U_alpha']

##alpha dose rate 
alpha.U <- U*keff*(conversion.alpha.U.preRadon+(1-Radon.loss)*conversion.alpha.U.afterRadon)
alpha.U.error <- alpha.U*sqrt((U.error/U)^2+(keff.err/keff)^2+(conversion.alpha.U.preRadon.err^2+(1-Radon.loss)^2*conversion.alpha.U.afterRadon.err^2)/(conversion.alpha.U.preRadon+(1-Radon.loss)*conversion.alpha.U.afterRadon)^2)
alpha.Th <- Th*keff*conversion.alpha.Th
alpha.Th.error <- alpha.Th*sqrt((Th.error/Th)^2+(keff.err/keff)^2+(conversion.alpha.Th.err/conversion.alpha.Th)^2)  
alpha <- (alpha.U+alpha.Th)/(1+1.50*water)
#alpha.error <- sqrt((1/(1+1.5*water))^2*alpha.U.error^2+(1/(1+1.5*water))^2*alpha.Th.error^2+((alpha.U+alpha.Th)/(1+1.5*water)^2)^2*1.5^2*water.error^2)
alpha.error <- alpha*sqrt((alpha.U.error^2 + alpha.Th.error^2)/(alpha.U+alpha.Th)^2 + (1.5*water.error)^2/(1+1.5*water)^2)

c(alpha, alpha.error)

}