##Dose rate calculation for calcite and dolomite, with effective k-value, assuming infinite homogeneous environment


#######  Prepare the CSV template, set the work directory, then select all the script and run  ################

rm(list = ls())
setwd("D:/R scripts for dose rate calculation/R scripts and templates/R scripts for calcite and dolomite")
#change the working directory to the folder address the user is using.

Dr<-read.csv("Template with keff.csv",header = TRUE)

filename<-"Dr results from keff.csv" 

###########################################
Radon.loss<-0.0     #for radon loss. 0.2 means 20% of the radon is escaped.
###########################################


###conversion factors are from Liritzis 2013 (Table 1-3).
###U unit 1 ppm
conversion.alpha.U.preRadon<-1.265
conversion.alpha.U.preRadon.err<-0.011
conversion.alpha.U.afterRadon<-2.793-1.265
conversion.alpha.U.afterRadon.err<-0.011

conversion.beta.U.preRadon<-0.0602
conversion.beta.U.preRadon.err<-0.0002
conversion.beta.U.afterRadon<-0.1459-0.0602
conversion.beta.U.afterRadon.err<-0.0004

conversion.gamma.U.preRadon<-0.0044
conversion.gamma.U.preRadon.err<-0.0000
conversion.gamma.U.afterRadon<-0.1118-0.0044
conversion.gamma.U.afterRadon.err<-0.0002

###Th unit 1ppm
conversion.alpha.Th<-0.7375
conversion.alpha.Th.err<-0.0026

conversion.beta.Th<-0.0275
conversion.beta.Th.err<-0.0009

conversion.gamma.Th<-0.0481
conversion.gamma.Th.err<-0.0002

### K unit 1%
conversion.beta.K<-0.8011
conversion.beta.K.err<-0.0073

conversion.gamma.K<-0.2498
conversion.gamma.K.err<-0.0048


##Rb  for 50 ppm Rb (may be used)
conversion.beta.Rb<-0.0185
conversion.beta.Rb.err<-0.0004
  


for(i in 1:nrow(Dr))  ##recycling from the first row to the last row.
    {
  geo.lat<-Dr[i,2]        #latitude
  geo.lon<-Dr[i,3]        #longitude
  altitude<-Dr[i,4]        #altitude，unit km
  depth<-Dr[i,5]           #the sampling depth from the surface，unit m

  water=Dr[i,6]/100     
  water.error=Dr[i,7]/100
  U=Dr[i,8]
  U.error=Dr[i,9]
  Th=Dr[i,10]
  Th.error=Dr[i,11]
  K=Dr[i,12] 
  K.error=Dr[i,13]
  alpha.factor<-Dr[i,14]   
  alpha.factor.err<-Dr[i,15]  

######################################### cosmic ray dose rate calculation #############################
  
  geomag.lat=180/pi*asin(0.203*cos(geo.lat*pi/180)*cos((geo.lon-291)*pi/180)+0.979*sin(geo.lat*pi/180))
  geomag.lat.abs<-abs(geomag.lat)
  #To set the three parameters Fi，J，H in the function D=D0*[Fi+J*exp(-h/H)]
  #from Prescott and Hutton，1994 (Fig.2).Do not change these values.
  if(geomag.lat.abs<2.5)                     {Fi=0.40; J=0.52; H=4.40}
  if(2.5<=geomag.lat.abs & geomag.lat.abs<7.5)   {Fi=0.39; J=0.54; H=4.38}
  if(7.5<=geomag.lat.abs & geomag.lat.abs<12.5)  {Fi=0.38; J=0.55; H=4.36}
  if(12.5<=geomag.lat.abs & geomag.lat.abs<17.5) {Fi=0.37; J=0.57; H=4.33}
  if(17.5<=geomag.lat.abs & geomag.lat.abs<22.5) {Fi=0.35; J=0.60; H=4.30}
  if(22.5<=geomag.lat.abs & geomag.lat.abs<27.5) {Fi=0.32; J=0.63; H=4.25}
  if(27.5<=geomag.lat.abs & geomag.lat.abs<32.5) {Fi=0.28; J=0.69; H=4.18}
  if(32.5<=geomag.lat.abs & geomag.lat.abs<37.5) {Fi=0.25; J=0.75; H=4.12}
  if(37.5<=geomag.lat.abs)                       {Fi=0.24; J=0.76; H=4.10}
  
  #calculate the cosmic ray
  x=2.0*depth   ##Density is set as 2.0 g/cm3, unit 100 g/cm-2
  if(x<1.7) {D0=0.115*exp(-x*100/77.86)+0.172}      #fitted from data of Prescott and Hutton (1988), containing soft and hard components
  else {D0=exp(-0.00055*x)*6072/(((x+11.6)^1.68+75)*(x+212))} #Prescott and Hutton (1994), only  hard component is significant
  
  cosmic.ray=D0*(Fi+J*exp(altitude/1000/H)) 
  Dr[i,16]=cosmic.ray
  cosmic.error=0.1*cosmic.ray #usually it is taken as 10 %
  Dr[i,17]=cosmic.error


  sampleID = as.character(Dr[i,1])
  
  cat('*******************************   Sample ID:  ',sampleID, '*******************************************************************','\n')
 
    
  #################################################### alpha dose rate #########################################
  
  alpha.U=U*alpha.factor*(conversion.alpha.U.preRadon+(1-Radon.loss)*conversion.alpha.U.afterRadon)
  alpha.U.error=alpha.U*sqrt((U.error/U)^2+(alpha.factor.err/alpha.factor)^2+(conversion.alpha.U.preRadon.err^2+(1-Radon.loss)^2*conversion.alpha.U.afterRadon.err^2)/(conversion.alpha.U.preRadon+(1-Radon.loss)*conversion.alpha.U.afterRadon)^2)
  cat('Individual alpha dose rate before water:','\n',' alpha.U=',alpha.U,'+-',alpha.U.error,'Gy','\n')
  
  alpha.Th=Th*alpha.factor*conversion.alpha.Th
  alpha.Th.error=alpha.Th*sqrt((Th.error/Th)^2+(alpha.factor.err/alpha.factor)^2+(conversion.alpha.Th.err/conversion.alpha.Th)^2)  
  cat('  alpha.Th=',alpha.Th,'+-',alpha.Th.error,'Gy','\n')
  
  alpha=(alpha.U+alpha.Th)/(1+1.50*water)
  alpha.error=sqrt((1/(1+1.5*water))^2*alpha.U.error^2+(1/(1+1.5*water))^2*alpha.Th.error^2+((alpha.U+alpha.Th)/(1+1.5*water)^2)^2*1.5^2*water.error^2)
  cat('alpha dose rate after water=',alpha,'+-',alpha.error,'Gy','\n\n')
  
  Dr[i,18]=alpha
  Dr[i,19]=alpha.error 
  
########################################################## beta dose rate calculation ################################# 

  beta.K=K*conversion.beta.K;
  beta.K.error=beta.K*sqrt((K.error/K)^2+(conversion.beta.K.err/conversion.beta.K)^2)
  cat('Individual beta dose rate before water:','\n',' beta.K=',beta.K,'+-',beta.K.error,'Gy','\n')
  
  beta.U=U*(conversion.beta.U.preRadon+(1-Radon.loss)*conversion.beta.U.afterRadon);
  beta.U.error=beta.U*sqrt((U.error/U)^2+(conversion.beta.U.preRadon.err^2+(1-Radon.loss)^2*conversion.beta.U.afterRadon.err^2)/(conversion.beta.U.preRadon+(1-Radon.loss)*conversion.beta.U.afterRadon)^2)
  cat('  beta.U=',beta.U,'+-',beta.U.error,'Gy','\n')
  
  beta.Th=Th*conversion.beta.Th;
  beta.Th.error=beta.Th*sqrt((Th.error/Th)^2+(conversion.beta.Th.err/conversion.beta.Th)^2)
  cat('  beta.Th=',beta.Th,'+-',beta.Th.error,'Gy','\n')
  
  beta=(beta.K+beta.U+beta.Th)/(1+1.25*water);
  beta.error=sqrt((1/(1+1.25*water))^2*beta.K.error^2+(1/(1+1.25*water))^2*beta.U.error^2+(1/(1+1.25*water))^2*beta.Th.error^2+((beta.K+beta.U+beta.Th)/(1+1.25*water)^2)^2*1.25^2*water.error^2)
  cat('beta dose rate after water=',beta,'+-',beta.error,'Gy','\n\n')
  
    Dr[i,20]=beta
    Dr[i,21]=beta.error
 
########################################################## gama dose rate calculation ################################# 
    gamma.K=K*conversion.gamma.K
    gamma.K.error=gamma.K*sqrt((K.error/K)^2+(conversion.gamma.K.err/conversion.gamma.K)^2)
    cat('Individual gamma dose rate before water:','\n',
        ' gamma.K=',gamma.K,'+-',gamma.K.error,'Gy','\n')
    
    gamma.U=U*(conversion.gamma.U.preRadon+(1-Radon.loss)*conversion.gamma.U.afterRadon)
    gamma.U.error=gamma.U*sqrt((U.error/U)^2+(conversion.gamma.U.preRadon.err^2+(1-Radon.loss)^2*conversion.gamma.U.afterRadon.err^2)/(conversion.gamma.U.preRadon+(1-Radon.loss)*conversion.gamma.U.afterRadon)^2)
    cat('  gamma.U=',gamma.U,'+-',gamma.U.error,'Gy','\n')
    
    gamma.Th=Th*conversion.gamma.Th
    gamma.Th.error=gamma.Th*sqrt((Th.error/Th)^2+(conversion.gamma.Th.err/conversion.gamma.Th)^2)
    cat('  gamma.Th=',gamma.Th,'+-',gamma.Th.error,'Gy','\n')
    
    gamma=(gamma.K+gamma.U+gamma.Th)/(1+1.14*water);
    gamma.error=sqrt((1/(1+1.14*water))^2*gamma.K.error^2+(1/(1+1.14*water))^2*gamma.U.error^2+(1/(1+1.14*water))^2*gamma.Th.error^2+((gamma.K+gamma.U+gamma.Th)/(1+1.14*water)^2)^2*1.14^2*water.error^2)
    
    Dr[i,22]=gamma
    Dr[i,23]=gamma.error
    cat('gamma dose rate after water=',gamma,'+-',gamma.error,'Gy','\n\n')
    
   
###########################################  total dose rate ################################################################
   Doserate=cosmic.ray+alpha+beta+gamma
   Doserate.error=sqrt(cosmic.error^2+alpha.error^2+beta.error^2+gamma.error^2)
   Dr[i,24]=Doserate
   Dr[i,25]=Doserate.error
   

    if(i==nrow(Dr)) 
      {
      
      write.csv(Dr,filename,row.names=FALSE)
      cat('\n',"Well done!","\n",'The results are written in the file:',filename,"\n",
          "The dose rates and errors are also displayed below.","\n\n")
      print(Dr[,c(1,16,17,18,19,20,21,22,23,24,25)])
      }
}
