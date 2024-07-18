##########This code was authored by Devin W. Goodsman, postdoctoral researcher at LANL########

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

# Here I run the IMAP function to determine the flight phenology as a function of
# the distribution of previous flight times.

# To do this more efficiently we use the parallel package
library(parallel)

source("C:/Users/Grace.Veenstra/Downloads/ece33590-sup-0004-supinfo/IMAPfunctionfixedLCPredFlight10.R")

MaxT2 = read.csv("C:/Users/Grace.Veenstra/Downloads/ece33590-sup-0004-supinfo/MaxDailyT08to16JunOctVal.csv", h = T)
MinT2 = read.csv("C:/Users/Grace.Veenstra/Downloads/ece33590-sup-0004-supinfo/MinDailyT08to16JunOctVal.csv", h = T)
PD2 = read.csv("C:/Users/Grace.Veenstra/Downloads/ece33590-sup-0004-supinfo/PhenologyTemp08to16JunOctVal2.csv", h = T)

colnames(MaxT2)
MaxT2 = MaxT2[,-1]
MinT2 = MinT2[,-1]
PD2 = PD2[,-1]

#################################################################################################
# Generating distributed start times for the 2008 data

PD2008 = subset(PD2, Year == "2008")
PD2008S = subset(PD2008, Project == "Smoky")
PD2008P = subset(PD2008, Project == "Peace")

# First I import the trap data for 2007  [Can't run because MPB doesn't exist in ece33590 folder]
Trap07PS = read.csv("MPBFunnelTraps2007PS.csv", h = T)

# Setting up the forcing for Peace
ForcingP08 = rep(0,504)
dim(Trap07PS)
for(i in 1:dim(Trap07PS)[1]) ForcingP08[Trap07PS$TimeTwo[i]] = Trap07PS$Peace[i]

ForcingP08Mat = matrix(rep(ForcingP08, dim(PD2008P)[1]), 
                       nrow = dim(PD2008P)[1], ncol = length(ForcingP08),
                       byrow = TRUE)

PD2008P2 = cbind(PD2008P, ForcingP08Mat)
rm(ForcingP08Mat); gc()

# Setting up the forcing for Smoky
ForcingS08 = rep(0,504)
dim(Trap07PS)
for(i in 1:dim(Trap07PS)[1]) ForcingS08[Trap07PS$TimeTwo[i]] = Trap07PS$Smoky[i]

ForcingS08Mat = matrix(rep(ForcingS08, dim(PD2008S)[1]),
                       nrow = dim(PD2008S)[1], ncol = length(ForcingP08),
                       byrow = TRUE)

PD2008S2 = cbind(PD2008S, ForcingS08Mat)
rm(ForcingS08Mat); gc()

#################################################################################################
# Generating distributed start times for the 2009 data

PD2009 = subset(PD2, Year == "2009")
PD2009S = subset(PD2009, Project == "Smoky")
PD2009P = subset(PD2009, Project == "Peace")

# First I import the trap data for 2008
Trap08PS = read.csv("MPBFunnelTraps2008PS2.csv", h = T)
Trap08PS$ProbsSmoky = rowSums(Trap08PS[,3:5])
Trap08PS$ProbsPeace = rowSums(Trap08PS[,6:7])

# Setting up the forcing for Peace
ForcingP09 = rep(0,504)
dim(Trap08PS)
for(i in 1:dim(Trap08PS)[1]) ForcingP09[Trap08PS$TimeTwo[i]] = Trap08PS$ProbsPeace[i]

ForcingP09Mat = matrix(rep(ForcingP09, dim(PD2009P)[1]), 
                       nrow = dim(PD2009P)[1], ncol = length(ForcingP09),
                       byrow = TRUE)

PD2009P2 = cbind(PD2009P, ForcingP09Mat)
rm(ForcingP09Mat); gc()

# Setting up the forcing for Smoky
ForcingS09 = rep(0,504)
dim(Trap08PS)
for(i in 1:dim(Trap08PS)[1]) ForcingS09[Trap08PS$TimeTwo[i]] = Trap08PS$ProbsSmoky[i]

ForcingS09Mat = matrix(rep(ForcingS09, dim(PD2009S)[1]),
                       nrow = dim(PD2009S)[1], ncol = length(ForcingS09),
                       byrow = TRUE)

PD2009S2 = cbind(PD2009S, ForcingS09Mat)
rm(ForcingS09Mat); gc()

#################################################################################################
# Generating distributed start times for the 2010 data

PD2010 = subset(PD2, Year == "2010")
PD2010S = subset(PD2010, Project == "Smoky")
PD2010P = subset(PD2010, Project == "Peace")

# First I import the trap data for 2009
Trap09PS = read.csv("MPBFunnelTraps2009PS2.csv", h = T)
Trap09PS$ProbsSmoky = rowSums(Trap09PS[,3:5])
Trap09PS$ProbsPeace = rowSums(Trap09PS[,6:7])

# Setting up the forcing for Peace
ForcingP10 = rep(0,504)
dim(Trap09PS)
for(i in 1:dim(Trap09PS)[1]) ForcingP10[Trap09PS$TimeTwo[i]] = Trap09PS$ProbsPeace[i]

ForcingP10Mat = matrix(rep(ForcingP10, dim(PD2010P)[1]), 
                       nrow = dim(PD2010P)[1], ncol = length(ForcingP10),
                       byrow = TRUE)

PD2010P2 = cbind(PD2010P, ForcingP10Mat)
rm(ForcingP10Mat); gc()

# Setting up the forcing for Smoky
ForcingS10 = rep(0,504)
dim(Trap09PS)
for(i in 1:dim(Trap09PS)[1]) ForcingS10[Trap09PS$TimeTwo[i]] = Trap09PS$ProbsSmoky[i]

ForcingS10Mat = matrix(rep(ForcingS10, dim(PD2010S)[1]),
                       nrow = dim(PD2010S)[1], ncol = length(ForcingS10),
                       byrow = TRUE)

PD2010S2 = cbind(PD2010S, ForcingS10Mat)
rm(ForcingS10Mat); gc()

#################################################################################################
# Generating distributed start times for the 2011 data

PD2011 = subset(PD2, Year == "2011")
PD2011S = subset(PD2011, Project == "Smoky")
PD2011P = subset(PD2011, Project == "Peace")

# First I import the trap data for 2010
Trap10PS = read.csv("MPBFunnelTraps2010PS2.csv", h = T)
Trap10PS$ProbsSmoky = rowSums(Trap10PS[,3:5])
Trap10PS$ProbsPeace = rowSums(Trap10PS[,6:7])

# Setting up the forcing for Peace
ForcingP11 = rep(0,504)
dim(Trap10PS)
for(i in 1:dim(Trap10PS)[1]) ForcingP11[Trap10PS$TimeTwo[i]] = Trap10PS$ProbsPeace[i]

ForcingP11Mat = matrix(rep(ForcingP11, dim(PD2011P)[1]), 
                       nrow = dim(PD2011P)[1], ncol = length(ForcingP11),
                       byrow = TRUE)

PD2011P2 = cbind(PD2011P, ForcingP11Mat)
rm(ForcingP11Mat); gc()

# Setting up the forcing for Smoky
ForcingS11 = rep(0,504)
dim(Trap10PS)
for(i in 1:dim(Trap10PS)[1]) ForcingS11[Trap10PS$TimeTwo[i]] = Trap10PS$ProbsSmoky[i]

ForcingS11Mat = matrix(rep(ForcingS11, dim(PD2011S)[1]),
                       nrow = dim(PD2011S)[1], ncol = length(ForcingS11),
                       byrow = TRUE)

PD2011S2 = cbind(PD2011S, ForcingS11Mat)
rm(ForcingS11Mat); gc()

#################################################################################################
# Generating distributed start times for the 2012 data

PD2012 = subset(PD2, Year == "2012")
PD2012W = subset(PD2012, Project == "Woodlands")

# First I import the trap data for 2011
Trap11W = read.csv("MPBFunnelTraps2011W.csv", h = T)

# Setting up the forcing for Woodlands
ForcingS12 = rep(0,504)
dim(Trap11W)
for(i in 1:dim(Trap11W)[1]) ForcingS12[Trap11W$TimeTwo[i]] = Trap11W$Woodlands[i]

ForcingW12Mat = matrix(rep(ForcingS12, dim(PD2012W)[1]),
                       nrow = dim(PD2012W)[1], ncol = length(ForcingS12),
                       byrow = TRUE)

PD2012W2 = cbind(PD2012W, ForcingW12Mat)
rm(ForcingW12Mat); gc()

####################################################################################################
# Generating random start times for the 2013 data

PD2013 = subset(PD2, Year == "2013")
PD2013F = subset(PD2013, Project == "Foothills")
PD2013W = subset(PD2013, Project == "Woodlands")

# First I import the trap data for 2012
Trap12FW = read.csv("MPBFunnelTraps2012FW.csv", h = T)

# Setting up the forcing for Foothills
ForcingF13 = rep(0,504)
dim(Trap12FW)
for(i in 1:dim(Trap12FW)[1]) ForcingF13[Trap12FW$TimeTwo[i]] = Trap12FW$Foothills[i]

ForcingF13Mat = matrix(rep(ForcingF13, dim(PD2013F)[1]), 
                       nrow = dim(PD2013F)[1], ncol = length(ForcingF13),
                       byrow = TRUE)

PD2013F2 = cbind(PD2013F, ForcingF13Mat)
rm(ForcingF13Mat); gc()

# Setting up the forcing for Woodlands
ForcingW13 = rep(0,504)
dim(Trap12FW)
for(i in 1:dim(Trap12FW)[1]) ForcingS11[Trap12FW$TimeTwo[i]] = Trap12FW$Woodlands[i]

ForcingW13Mat = matrix(rep(ForcingW13, dim(PD2013W)[1]),
                       nrow = dim(PD2013W)[1], ncol = length(ForcingW13),
                       byrow = TRUE)

PD2013W2 = cbind(PD2013W, ForcingW13Mat)
rm(ForcingW13Mat); gc()

####################################################################################################
# Generating distributed start times for the 2014 data

PD2014 = subset(PD2, Year == "2014")
PD2014P = subset(PD2014, Project == "Peace")

# First I import the trap data for 2013
Trap13P = read.csv("MPBFunnelTraps2013P.csv", h = T)

# Setting up the forcing for Peace
ForcingP14 = rep(0,504)
dim(Trap13P)
for(i in 1:dim(Trap13P)[1]) ForcingP14[Trap13P$TimeTwo[i]] = Trap13P$Peace[i]

ForcingP14Mat = matrix(rep(ForcingP14, dim(PD2014P)[1]),
                       nrow = dim(PD2014P)[1], ncol = length(ForcingP14),
                       byrow = TRUE)

PD2014P2 = cbind(PD2014P, ForcingP14Mat)
rm(ForcingP14Mat); gc()

####################################################################################################
# Generating distributed start times for the 2015 data

PD2015 = subset(PD2, Year == "2015")
PD2015P = subset(PD2015, Project == "Peace")

# First I import the trap data for 2014
Trap14P = read.csv("MPBFunnelTraps2014P.csv", h = T)
Trap14P$ProbsP = rowSums(Trap14P[,3:10])

# Setting up the forcing for Peace
ForcingP15 = rep(0,504)
dim(Trap14P)
for(i in 1:dim(Trap14P)[1]) ForcingP15[Trap14P$TimeTwo[i]] = Trap14P$ProbsP[i]

ForcingP15Mat = matrix(rep(ForcingP15, dim(PD2015P)[1]),
                       nrow = dim(PD2015P)[1], ncol = length(ForcingP15),
                       byrow = TRUE)

PD2015P2 = cbind(PD2015P, ForcingP15Mat)
rm(ForcingP15Mat); gc()

####################################################################################################
# Generating distributed start times for the 2016 data

PD2016 = subset(PD2, Year == "2016")
PD2016P = subset(PD2016, Project == "Peace")

# First I import the trap data for 2015
Trap15P = read.csv("MPBFunnelTraps2015P.csv", h = T)
Trap15P$ProbsP = rowSums(Trap15P[,3:4])

# Setting up the forcing for Peace
ForcingP16 = rep(0,504)
dim(Trap15P)
for(i in 1:dim(Trap15P)[1]) ForcingP16[Trap15P$TimeTwo[i]] = Trap15P$ProbsP[i]

ForcingP16Mat = matrix(rep(ForcingP16, dim(PD2016P)[1]),
                       nrow = dim(PD2016P)[1], ncol = length(ForcingP16),
                       byrow = TRUE)

PD2016P2 = cbind(PD2016P, ForcingP16Mat)
rm(ForcingP16Mat); gc()

##################################################################################################
## Now all of these data with vectors of 100 start times need to be bound back together.

PD3 = rbind(PD2008P2, PD2008S2, PD2009P2, PD2009S2, PD2010P2, PD2010S2, PD2011P2, PD2011S2,
            PD2012W2, PD2013F2, PD2013W2, PD2014P2,  PD2015P2,  PD2016P2)

SubIndex = row.names(PD3)

MaxT2 = MaxT2[SubIndex,]
MinT2 = MinT2[SubIndex,]

## Now I will test the IMAPfunction
Tmin1 = as.numeric(MinT2[150,7:510])
Tmax1 = as.numeric(MaxT2[150,7:510])

ptm <- proc.time()
Out1 = ImapFuncLC(Tmin = Tmin1, Tmax = Tmax1, StartT = as.numeric(PD3[150, 19:118]))
proc.time() - ptm

Out1
rm(Tmin1, Tmax1, Out1); gc()

# Making a function that can be parallelized
Index = 1:(dim(PD3)[1])

ParFunc1 = function(s){
  
  Forcing1 = as.numeric(PD3[s, 7:510])
  Tmin1 = as.numeric(MinT2[s,7:510])
  Tmax1 = as.numeric(MaxT2[s,7:510])
  
  Sim = ImapFuncLC(Tmin = Tmin1, Tmax = Tmax1, Forcing = Forcing1)
  
  return(Sim)
}

#################################################################################
#########################Here's where I run the code on clusters#################
# I need to make sure required variables and functions are available to 
# the workers (clusterExport function)

M = makeCluster(spec=4)

junk1 = clusterExport(M, varlist = c("PD3", "MinT2", "MaxT2", "ImapFuncLC"), 
                      envir = .GlobalEnv)

ptm <- proc.time()
Out1 = parSapply(M, X = Index, FUN = ParFunc1, simplify = TRUE)
proc.time() - ptm
# This takes about 1 minute on a laptop with 2 cores (four with hyperthreading).

stopCluster(M)

dim(Out1)
Out1 = t(Out1)
str(Out1)

colnames(PD3)
PD4 = cbind(PD3[,1:6], Out1)
str(PD4)
colnames(PD4)

# Now I will write the PD4 data frame to csv
write.csv(PD4, "C:/Users/Grace.Veenstra/Downloads/ece33590-sup-0004-supinfo/Grace_Outputs/PhenologyTemp08to16LCFlightPred10.csv")
