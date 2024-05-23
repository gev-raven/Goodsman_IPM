###This code was originally authored by Devin W. Goodsman, postdoctoral researcher at LANL#####

##############################################################################################
# This software has been authored by an employee or employees of Los Alamos National Security,
# LLC, operator of the Los Alamos National Laboratory (LANL) under contract No. 
# DE-AC52-06NA25396 with the U.S. Department of Energy. The U.S. government has rights to use,
# reproduce, and distribute this software. The public may copy, distribute, prepare derivative
# works and publicly display this software without charge, provided that this Notice and
# any statement of authorship are reproduced in all copies. Neither the government nor LANS
# makes any warranty, express or implied, or assumes any liability or responsibility for the
# use of this software. If software is modified to produce derivative works, such modified
# software should be clearly marked, so as not to confuse it with the version available from
# LANL
##############################################################################################

# Here is the the algorithm for the individual-based model 
# of mountain pine beetle phenology.

# This algorithm runs only a single stochastic run of the simulation.

# This version uses real temperature data from the Jasper
# warden site in Alberta.

######################################################################
# The rate equations are given in Regniere 2012
# I assume (as they did) that there is log-normal variability
# in development rates. The output is saved to a csv file 
######################################################################

# Reading in the temperature data (in degrees C)
Tmin = read.csv("JasperDailyMin2.csv", h = F)   # Min temp  
Tmax = read.csv("JasperDailyMax2.csv", h = F)   # Max temp

# Converting temperatures to under-bark temperatures
Temp = 0.5*(Tmax$V1 + Tmin$V1) + 0.9 + 6.6*(Tmax$V1 - Tmin$V1)/(2*24.4)
Tmax2 = Tmax$V1 + 6.6*(Tmax$V1 - Tmin$V1)/24.4
Tmin2 = Tmin$V1 + 1.8

#########################################################################################
# Components of the phenology model
#########################################################################################

# The Regniere function for temperature dependent beetle development
RegniereFunc = function(TC, TB, DeltaB, TM, DeltaM, omega, psi){
  
  Y = 0.0
  
  if(TC >= TB & TC <= TM){
    Y = psi*(exp(omega*(TC - TB)) - (TM - TC)/(TM - TB)*exp(-omega*(TC - TB)/DeltaB)
             - (TC - TB)/(TM - TB)*exp(omega*(TM - TB) - (TM - TC)/DeltaM))
  }
  
  if(Y <= 0.0) Y = 0.0
  
  return(Y)
}

# Defining all of the parameters

# Defining the time step
# We will assume that the reaction rate is measured per day
# like the time step
deltat = 1.0 # units are days

# The number of initial eggs layed in the brood
imax = 82

# parameters for oviposition rate
# (from Regniere et al 2012)
sigma0 = 0.2458         # controls rate variability
TB0 = 4.6341            # base temperature in degrees C
DeltaB0 = 0.1
TM0 = 27.7587           # max temperature in degree C
DeltaM0 = 3.0759
omega0 = 0.3684
psi0 = 0.005199

# parameters for development rate for eggs
# (from Regniere et al 2012)
sigma1 = 0.1799         # controls rate variability
TB1 = 7.0               # base temperature in degrees C
DeltaB1 = 0.019297569
TM1 = 30.0928           # max temperature in degree C
DeltaM1 = 4.4175
omega1 = 0.2563
psi1 = 0.02317

# parameters for development rate for L1 larvae
# (from Regniere et al 2012)
sigma2 = 0.2911
TB2 = 3.5559
DeltaB2 = 0.1
TM2 = 29.2647
DeltaM2 = 3.8227
omega2 = 0.2398
psi2 = 0.01082

# parameters for development rate for L2 larvae
# (from Regniere et al 2012)
sigma3 = 0.3799
TB3 = 6.9598
DeltaB3 = 0.097087379
TM3 = 28.9047
DeltaM3 = 3.0374
omega3 = 0.3714
psi3 = 0.01072

# parameters for development rate for L3 larvae
# (from Regniere et al 2012)
sigma4 = 0.3868
TB4 = 6.8462
DeltaB4 = 0.1
TM4 = 28.7013
DeltaM4 = 2.5359
omega4 = 0.4399
psi4 = 0.003892

# parameters for development rate for L4 larvae
# (from Regniere et al 2012)
sigma5 = 0.3932
TB5 = 16.2464
DeltaB5 = 0.039052889
TM5 = 28.0
DeltaM5 = 4.5504
omega5 = 0.2593
psi5 = 0.05034

# parameters for development rate for pupae
# (from Regniere et al 2012)
sigma6 = 0.2998
TB6 = 5.63
DeltaB6 = 0.10989011
TM6 = 28.55
DeltaM6 = 2.86
omega6 = 0.1532
psi6 = 0.02054

# parameters for development rate for teneral adults
# (from Regniere et al 2012)
sigma7 = 0.5284
TB7 = 4.24
DeltaB7 = 0.099967011
TM7 = 35.0
DeltaM7 = 7.1479
omega7 = 0.1463
psi7 = 0.01173

###################################################################################
# Initializing the integrated model
###################################################################################

# The start time
StartT = 15

# Containers for predictions
Fec = rep(0,443)
Fec[1:StartT] = rep(82, StartT)
Eggs = rep(0,443)
L1 = rep(0,443)
L2 = rep(0,443)
L3 = rep(0,443)
L4 = rep(0,443)
Pupae = rep(0,443)
TenAd = rep(0,443)
Adults = rep(0,443)

# Containers for level of development
FecVec = rep(0,82)
EggVec = rep(0,82)
L1Vec = rep(0,82)
L2Vec = rep(0,82)
L3Vec = rep(0,82)
L4Vec = rep(0,82)
PupaeVec = rep(0,82)
TeneralVec = rep(0,82)

# A container for survival (1 -> survived and 0 -> did not survive)
SurvVec = rep(1,82)

#################################################################################
# The stochastic individual-based model
#################################################################################
ptm <- proc.time()
# Doing the iteration
for(i in (StartT+1):443){
  
  # Making an adjustment to minimum air temperature to account for the buffering
  # effect of tree bark as described by Bolstad, Bentz, and Logan (1997)
  # This is used in many of the mortality calculations
  Tmn = Tmin2[i]
  
  # Computing the median development rate for each life stage in this time step
  median0 = RegniereFunc(TC = Temp[i], TB0, DeltaB0, TM0, DeltaM0, omega0, psi0)   # for oviposition
  median1 = RegniereFunc(TC = Temp[i], TB1, DeltaB1, TM1, DeltaM1, omega1, psi1)   # for eggs
  median2 = RegniereFunc(TC = Temp[i], TB2, DeltaB2, TM2, DeltaM2, omega2, psi2)   # for L1
  median3 = RegniereFunc(TC = Temp[i], TB3, DeltaB3, TM3, DeltaM3, omega3, psi3)   # for L2
  median4 = RegniereFunc(TC = Temp[i], TB4, DeltaB4, TM4, DeltaM4, omega4, psi4)   # for L3
  median5 = RegniereFunc(TC = Temp[i], TB5, DeltaB5, TM5, DeltaM5, omega5, psi5)   # for L4
  median6 = RegniereFunc(TC = Temp[i], TB6, DeltaB6, TM6, DeltaM6, omega6, psi6)   # for pupae
  median7 = RegniereFunc(TC = Temp[i], TB7, DeltaB7, TM7, DeltaM7, omega7, psi7)   # for teneral adults
  
  # The mu parameter of the lognormal distribution is given by ln(R[T]*deltat)
  mu0 = log(median0*deltat) # for oviposition
  mu1 = log(median1*deltat) # for eggs
  mu2 = log(median2*deltat) # for L1
  mu3 = log(median3*deltat) # for L2
  mu4 = log(median4*deltat) # for L3
  mu5 = log(median5*deltat) # for L4
  mu6 = log(median6*deltat) # for pupae
  mu7 = log(median7*deltat) # for teneral adults
  
  # In each time step we set the counters back to zero
  EggCount = 0.0
  L1Count = 0.0
  L2Count = 0.0
  L3Count = 0.0
  L4Count = 0.0
  PupaeCount = 0.0
  TeneralCount = 0.0
  AdultCount = 0.0
  
  for(j in 1:imax){ # imax is set to 82 (82 individuals)
    #R0 = rlnorm(1, meanlog = mu0, sdlog = sigma0)
    R1 = rlnorm(1, meanlog = mu1, sdlog = sigma1)
    R2 = rlnorm(1, meanlog = mu2, sdlog = sigma2)
    R3 = rlnorm(1, meanlog = mu3, sdlog = sigma3)
    R4 = rlnorm(1, meanlog = mu4, sdlog = sigma4)
    R5 = rlnorm(1, meanlog = mu5, sdlog = sigma5)
    R6 = rlnorm(1, meanlog = mu6, sdlog = sigma6)
    R7 = rlnorm(1, meanlog = mu7, sdlog = sigma7)
    
    # Now applying mortality to every stage except the larval stage (and oviposition).
    # Initially, we follow Regniere et al 2015 and assume that if the temperature
    # drops below -18.0 degrees C, there is complete mortality of
    # eggs, pupae, teneral adults and adults.
    if(EggVec[j] < 1.0 & Tmn <= -18.0){
      # Mortality for the egg stage
      SurvVec[j] = 0
    }
    if(L4Vec[j] >= 1.0 & PupaeVec[j] < 1.0 & Tmn <= -18.0){
      # Mortality for pupae that were already pupae
      SurvVec[j] = 0
    }
    if(PupaeVec[j] >= 1.0 & TeneralVec[j] < 1.0 & Tmn <= -18.0){
      # Mortality for teneral adults
      SurvVec[j] = 0
    }
    if(TeneralVec[j] >= 1.0 & Tmn <= -18.0){
      # Mortality for emerged adults
      SurvVec[j] = 0
    }
    
    # Advancing development levels
    # Note that we need to make sure that individuals don't immediately 
    # develop within a stage after just having developed into it.
    # Thus we record all the previous levels of development before updating.
    EggVecm1 = EggVec[j]
    L1Vecm1 = L1Vec[j]
    L2Vecm1 = L2Vec[j]
    L3Vecm1 = L3Vec[j]
    L4Vecm1 = L4Vec[j]
    PupaeVecm1 = PupaeVec[j]
    
    # We do a probablistic calculation for oviposition
    # in each step.
    Draw1 = runif(1)
    ProbOvip = (1 - exp(-median0*deltat))
    
    # Here is the oviposition step. We only allow oviposition if 
    # the ovipositing adult has not been killed by temperatures
    # below -18.0 degrees centigrade.
    if(EggVecm1 == 0.0 & median0 > 0.0 & Draw1 < ProbOvip & min(Tmin2[StartT:i]) > -18.0){
      EggVec[j] = 1e-200
    }
    
    # All of the following steps use similar logic
    # to one another (unlike oviposition)
    if(EggVecm1 > 0.0 & median1 > 0.0){
      EggVec[j] = EggVec[j] + R1
    }
    if(EggVecm1 >= 1.0 & median2 > 0.0){
      L1Vec[j] = L1Vec[j] + R2
    }
    if(L1Vecm1 >= 1.0 & median3 > 0.0){
      L2Vec[j] = L2Vec[j] + R3
    }
    if(L2Vecm1 >= 1.0 & median4 > 0.0){
      L3Vec[j] = L3Vec[j] + R4
    }
    if(L3Vecm1 >= 1.0 & median5 > 0.0){
      L4Vec[j] = L4Vec[j] + R5
    }
    if(L4Vecm1 >= 1.0 & median6 > 0.0){
      PupaeVec[j] = PupaeVec[j] + R6
    }
    if(PupaeVecm1 >= 1.0 & median7 > 0.0){
      TeneralVec[j] = TeneralVec[j] + R7
    }
    
    # Counting up the number of individuals in each stage
    if(EggVec[j] > 0.0 & EggVec[j] < 1.0 & SurvVec[j] == 1){
      # First counting the number of eggs
      EggCount = EggCount + 1
    }
    if(EggVec[j] >= 1.0 & L1Vec[j] < 1.0 & SurvVec[j] == 1){
      # Counting up the number of L1 individuals
      L1Count = L1Count + 1
    }
    if(L1Vec[j] >= 1.0 & L2Vec[j] < 1.0 & SurvVec[j] == 1){
      # Counting up the number of L2 individuals
      L2Count = L2Count + 1
    }
    if(L2Vec[j] >= 1.0 & L3Vec[j] < 1.0 & SurvVec[j] == 1){
      # Counting up the number of L3 individuals
      L3Count = L3Count + 1
    }
    if(L3Vec[j] >= 1.0 & L4Vec[j] < 1.0 & SurvVec[j] == 1){
      # Counting up the number of L4 individuals
      L4Count = L4Count + 1
    }
    if(L4Vec[j] >= 1.0 & PupaeVec[j] < 1.0 & SurvVec[j] == 1){
      # Counting up the number of Pupae
      PupaeCount = PupaeCount + 1
    }
    if(PupaeVec[j] >= 1.0 & TeneralVec[j] < 1.0 & SurvVec[j] == 1){
      # Counting up the number of teneral adults
      TeneralCount = TeneralCount + 1
    }
    if(TeneralVec[j] >= 1.0 & SurvVec[j] == 1){
      # Counting up the number of adults
      AdultCount = AdultCount + 1
    }
    
  }   # end of j-loop
  
  # Now we slot in the observed number of individuals in each life stage
  # at this time step (corresponding to i)
  Eggs[i] = EggCount
  L1[i] = L1Count
  L2[i] = L2Count
  L3[i] = L3Count
  L4[i] = L4Count
  Pupae[i] = PupaeCount
  TenAd[i] = TeneralCount
  Adults[i] = AdultCount
  
}   # end of i-loop
proc.time() - ptm

# user  system elapsed 
# 1.78    0.06   2.4
# So this takes about 2 seconds.

times = 1:443

# Now I write the data to csv (20 times)
#IBM = data.frame(times, Eggs, L1, L2, L3, L4, Pupae, TenAd, Adults)
#write.csv(IBM, "MPBPhenMortSim520R.csv")

graphics.off()
plot(Eggs ~ times, type = 'l', ylim = c(0,82))
lines(L1 ~ times, lty = 2)
lines(L2 ~ times, lty = 3)
lines(L3 ~ times, lty = 4)
lines(L4 ~ times, lty = 1)
lines(Pupae ~ times, lty = 1)
lines(TenAd ~ times, lty = 2)
lines(Adults ~ times, lty = 3)
# This looks approximately right.

# Let's compare it to the analytical solution
source("IMAPfunctionNewJ2.R")

# Reading in some temperature data (in degrees C)
Tmin3 = read.csv("JasperDailyMin2.csv", h = F)   # Min temp  
Tmax3 = read.csv("JasperDailyMax2.csv", h = F)   # Max temp

Tmin3 = as.numeric(Tmin3$V1)
Tmax3 = as.numeric(Tmax3$V1)

ptm <- proc.time()
Out1 = ImapFuncNewJ(Tmin = Tmin3, Tmax = Tmax3, StartT = 15)
proc.time() - ptm

plot(Eggs ~ times, type = 'l', ylim = c(0,82))
lines(Out1[2,] ~ times, lty = 1, col = 'blue')
plot(L1 ~ times, type = 'l', lty = 2)
lines(Out1[3,] ~ times, lty = 2, col = 'orange')
plot(L2 ~ times, type = 'l', lty = 1)
lines(Out1[4,] ~ times, lty = 2, col = 'green')
plot(L3 ~ times, type = 'l', lty = 1, ylim = c(0,82))
lines(Out1[5,] ~ times, lty = 2, col = 'yellow')
lines(L4 ~ times, type = 'l', lty = 1)
lines(Out1[6,] ~ times, lty = 2, col = 'grey')
plot(Pupae ~ times, type = 'l', lty = 1, ylim = c(0,82))
lines(Out1[7,] ~ times, lty = 2, col = 'purple')
plot(TenAd ~ times, type = 'l', lty = 1, ylim = c(0,82))
lines(Out1[8,] ~ times, lty = 2, col = 'purple')
plot(Adults ~ times, type = 'l', lty = 1, ylim = c(0,82))
lines(Out1[9,] ~ times, lty = 2, col = 'red')
