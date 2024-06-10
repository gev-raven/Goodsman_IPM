## Setting up model convolution testing


##############################################
###### Packages Set Up ----- 
##############################################
# install.packages("anytime")
# install.packages("tidyverse")
# install.packages("dbplyr")
# install.packages("dplyr")
# install.packages("dtplyr")
# install.packages("lubridate")

library(anytime)
library(tidyverse)
library(dbplyr)
library(dplyr)
library(dtplyr)
library(lubridate)


##############################################

## Reading in the temperature data (in degrees C)
temp_df <- read.csv("C:\\Users\\Grace.Veenstra\\Documents\\GitHub\\Goodsman_IPM\\Egg Modeling Test\\MarshCreek_Temp_2015-2016.csv",
                    header=TRUE) #takes csv file and creates dataframe
temp_df <- temp_df[c("Date","HUC_10","Mean")] 
  #removes 'min' 'max' and 'HUC_8' columns from frame
temp_df <- temp_df[which(temp_df$HUC_10 %in% "1706020503"),] 
  #filters to only get Marsh Creek (Creek Identifier: "1706020503")
temp_df$Date <- parse_date_time(temp_df$Date, orders=c('mdy','ymd')) 
  #reads different date formats
temp_df$Date <- format.Date(temp_df$Date, format="%Y-%m-%d")
  #sets the date format of Date list
temp_df$Ydate <- yday(temp_df$Date) 
  #convert date to 'day of year' (value 1 to 365)
temp_df$Jdate <- 1:nrow(temp_df) 
  #sets multi-year 'julian' values (1-730+)

temp <- temp_df$Mean #assigns variable to the mean temp in data
temp <- as.numeric(temp)


###### Functions -----
############################################

DevelopmentFunc = function(Tp, a, b, c, startspawn) {
  
  # develop_time is the development time (in days)
  # Tp is daily mean temperature
  
  develop_df <- temp_df
  develop_time <- exp(log(a)+log((Tp-c)^b))
  develop_rate <- (1/develop_time)
  develop_df$DailyDevelopment <- (1/develop_time) #daily development rate
  
  #define development period
 
  #develop_df <- subset(develop_df, develop_df$Jdate >= startspawn & develop_df$Jdate <= (startspawn + 366))
  ## restricts develop_df to the period of development
  
  #develop_df$TotalDevelopment <- cumsum(develop_df$DailyDevelopment)
    ## total development = sum of daily development
  
  #y <- develop_df[min(which((develop_df$TotalDevelopment) >= 1)), "Jdate"] 
    ## once development = 1, it is emergence time and we extract emergence date
  
  return(develop_rate)
  
}

EmergenceFunc = function(Tp, a, b, c, startspawn) {
  
  # develop_time is the development time (in days)
  # Tp is daily mean temperature
  
  develop_df <- temp_df
  develop_time <- exp(log(a)+log((Tp-c)^b))
  develop_df$DailyDevelopment <- (1/develop_time) #daily development rate
  
  if(startspawn %in% spawn.window) {
  
    develop_df <- subset(develop_df, 
                         develop_df$Jdate >= startspawn & develop_df$Jdate <= (startspawn + 366))
      ## restricts develop_df to the period of development
    
    develop_df$TotalDevelopment <- cumsum(develop_df$DailyDevelopment)
      ## total development = sum of daily development
    
    y <- develop_df[min(which((develop_df$TotalDevelopment) >= 1)), "Jdate"] 
      ## once development = 1, it is emergence time and we extract emergence date
    
    return(y)
    
  }else{
    
    y = 0
    
    return(y)
  }
  
  return(y)
  
}

LnormDist = function(x, mulog, sigmalog) {
  
  y_log = 1/(x*sigmalog*sqrt(2*pi))*exp(-1/(2*sigmalog^2)*(log(x)-mulog)^2)
  
  return(y_log)
}

ConvolveFunc = function(x1, y1, padsize) {
  
  #Padding
  x2 = c(x1, rep(0, padsize))
  y2 = c(y1, rep(0, padsize))
  
  #Fast Fourier
  fx2 = fft(x2)
  fy2 = fft(y2)

  #Convolution + Inverse Fourier Transform
  Convolution1 = Re(fft(fx2*fy2, inverse=T)/length(x2))
  
  #Clip padding
  Convolution2 = Convolution1[1:length(x1)]
  
  #Add in individuals that developd beyond next threshold
  Convolution2[length(x1)] = Convolution2[length(x1)] + 
    sum(Convolution1[(length(x1)+1):length(x2)])
  
  return(Convolution2)
  
}


# avec is just 'a vector' which functions to set a domain or physiological age
# by giving a vector against which to plot thresholds for life stages
  #avec = seq(1e-20, 2, length.out = 2^8) # domain for the larval stage
  avec = seq(1e-20, 2, 0.001)
  da = avec[3] - avec[2]
  avec1 <- min(which(avec >= 1)) #this allows you to vary the length of avec and still index the right cell below

# Figuring out where in the domain avec = 1 (upper breakpoint for egg stage)
# which.min(abs(avec - 1))
# avec[128] to avec[129]


##### Spawning Set-Up -----
######################################

## Spawning Distribution
spawn_mean <- 230 #aug 18; mean spawn start date [EXAMPLE]
spawn_sd <- 5 #spawning date standard deviation [EXAMPLE]
total.spawners <- 10000 #number of spawning fish [EXAMPLE]
spawn_dist <- rnorm(n = total.spawners, mean = spawn_mean, sd = spawn_sd) 
  #creates normal distribution of spawners
spawn_dist <- round(spawn_dist)
hist(spawn_dist)

## Spawn Window
spawndate_lower <- min(spawn_dist) #lower limit of spawn window
spawndate_upper <- max(spawn_dist) #upper limit of spawn window
spawn.window <- spawndate_lower:spawndate_upper

## Spawning Frequency
spawn_df <- data.frame(Date=temp_df$Date) #create 'spawn_df' data frame
spawn_df$Jdate <- temp_df$Jdate
spawn_df <- subset(spawn_df, spawndate_lower <= spawn_df$Jdate & spawn_df$Jdate <= spawndate_upper)
  #restricts spawn data frame only to the spawning window
SpawnFreq <- as.data.frame(table(spawn_dist, useNA="always")); names(SpawnFreq) <- c("Jdate","SpawnFreq")
  #creates frame with number of spawns per day from spawning distribution
spawn_df <- merge(spawn_df, SpawnFreq, by = "Jdate", all=T)
  #merges the 'SpawnFreq' dataframe with spawn_df based on the julian day
spawn_df[is.na(spawn_df$SpawnFreq),"SpawnFreq"] <- 0; spawn_df <- na.omit(spawn_df)
  #sets 'NA' values of 'number of spawners' to zero
spawn_df <- subset(spawn_df, min(which(spawn_df$SpawnFreq != 0)) | max(which(spawn_df$SpawnFreq != 0)))
  #Cuts any stray leading / tailing zeroes that can cause issues when calculating emergence


###### Parameters -----
############################################

# Defining the time step. Assume that the reaction rate is measured per day
deltat = 1.0 # units are days

## Egg Development Model
Tp <- temp
a <- exp(10.404)
b <- -2.043
c <- -7.575

## Egg Stage
sigma1 = 0.2


###### Initialization -----
############################################

# Prediction "Containers" / aka How Many Individuals in X life stage
Fec = rep(0,length(temp))
Fec[1:spawndate_lower] = rep(total.spawners, spawndate_lower)
Eggs = rep(0, length(temp))
Juv = rep(0, length(temp))

# Initializing the previous time step eggs
PrevEgg = rep(0,length(avec))
Egg.B = rep(0,length(avec))
NewEggsT1 = 0.0

# Initializing the previous time step juveniles
PrevJuv = rep(0,length(avec))
NewJuv = rep(0,length(avec))
NewJuvT1 = 0.0

# Assorted
develop_df <- temp_df
emergence <- NULL
DailyRate <- rep(0,length(temp))

###### Iteration -----
############################################

ptm <- proc.time()
for(i in 1:length(temp)) {
 
  ##### Spawning ------
  #######################
  
  ## Setting up 'fecundity'
  if(i %in% spawn.window){
    num.spawn <- spawn_df[which(spawn_df$Jdate == i), "SpawnFreq"] 
      #gets number of fish spawning at time step i
      #num.spawn not found outside this
  
    Fec[i] = 50
      # Fecundity is number of eggs each fish produces (say, 50)
    
  } else {
    num.spawn = 0
    Fec[i] = 0
  }
  
  ##### Egg Stage ------
  #######################
  
  # New eggs are individuals that developed into eggs in this time step
  # aka, eggs from spawning
  NewEggs = num.spawn*Fec[i]
  
  # Egg development rate
  
  if(i >= spawndate_lower){
    egg.rate = DevelopmentFunc(Tp = temp[i], a, b, c, startspawn = i) #daily development rate
    DailyRate[i] = egg.rate
    emergence[i] = EmergenceFunc(Tp, a, b, c, startspawn = i) #emergence day vector
  
  } else {
    egg.rate = 0
    DailyRate[i] = 0
    emergence[i] = 0
  }

  # To compute new juvenilles in next stage of time step, we collect those who were
  # already eggs in previous time step (Egg.B)
  PrevEgg = Egg.B
    # PrevEgg is the distribution of individuals of age 'f'
    # who advanced from age 'e' of the previous time step
  
  if(egg.rate > da) {
    
    #Aging Kernel - k_i(b-a)
    mu1 = log(egg.rate*deltat)
    egg.dist = LnormDist(x = avec, mulog = mu1, sigmalog = sigma1) #probability density distribution of the development rate
    egg.dist = egg.dist/sum(egg.dist) # normalizing
      # egg.dist is the 'aging kernal' of the convolution
      # It represents the probability of aging from age 'f' to age 'g'
      # Lisa: "adds noise"
    
    #Convolution = fft^-1[fft(Prev.Egg)*fft(egg.dist)]
    Egg.B = ConvolveFunc(x1 = PrevEgg, y1 = egg.dist, padsize = length(avec)) + NewEggsT1*egg.dist
      ## Egg.B is distribution of individuals of age 'g' in stage x at time i
    NewEggsT1 = NewEggs
      # Egg.B is the prob. density distribution of individuals of age 'g'
      # who advanced from age 'f' in the current time step

    #Computing Total Number of Ind. Currently in Egg Stage in this time step
    ### [Integral]
    Eggs[i] = sum(na.omit(Egg.B[1:(avec1-1)])) + NewEggs
  
  } else{
    Egg.B = PrevEgg
    NewEggsT1 = NewEggsT1 + NewEggs
    Eggs[i] = sum(na.omit(Egg.B[1:(avec1-1)])) + NewEggsT1
  }
  
  #### Juvenile Stage
  
  PrevJuv = sum(PrevEgg[avec1:length(avec)], na.rm = T)
  NewJuv = sum(Egg.B[avec1:length(avec)], na.rm = T) - PrevJuv
  Juv[i] = sum(na.omit(PrevJuv[1:(avec1-1)])) + NewJuv
  
}
proc.time() - ptm

Matrix1 = matrix(c(Fec, Eggs, Juv), nrow = 3, ncol=length(temp), byrow=TRUE)



## Emergence
if(i > spawndate_lower) {
  emergence[is.na(emergence)] <- 0
  emergence <- emergence[! emergence %in% c('0')]
  emerge.day <- rep(emergence, spawn_df$SpawnFreq) 
    #repeats a given emergence day by the number of spawns on the associated spawn day
  freq.emerge <- as.vector(table(emerge.day))
    #creates a vector that lists the number of fish emerging on a given day
  emerge.window <- min(emerge.day):max(emerge.day)
    #window of emergence days
  hist(emerge.day) #histogram displaying number of fish emerging on a given day
}


##### Graphing -----

times = temp_df$Jdate

#plot(DailyRate ~ times)
plot(Eggs ~ times)
plot(Eggs, xlim=range(emerge.window))


