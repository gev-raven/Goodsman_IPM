### Working on Perry + Juvenile


##### Packages Set Up -----

# install.packages("anytime")
# install.packages("tidyverse")
# install.packages("dbplyr")
# install.packages("dplyr")
# install.packages("dtplyr")
# install.packages("lubridate")
# install.packages("here")
# install.packages("zoo")

library(anytime)
library(tidyverse)
library(dbplyr)
library(dplyr)
library(dtplyr)
library(lubridate)
library(here)
library(zoo)


##### Temperature -----

## Reading in the temperature data (in degrees C)
temp_df <- read.csv("C:\\Users\\Grace.Veenstra\\Documents\\GitHub\\Goodsman_IPM\\Egg Modeling Test\\MarshCreek_Temp_2015-2018.csv", 
                    header=TRUE) 
  #takes csv file and creates dataframe
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
  #sets multi-year 'julian' values (1-1090+)

temp <- temp_df$Mean #assigns variable to the mean temp in data
temp <- as.numeric(temp)


##### Development Functions -----

DevelopmentFunc = function(Tp, a, b, c, startspawn) {
  
  # develop_time is the development time (in days)
  # Tp is daily mean temperature
  
  develop_df <- temp_df
  develop_time <- exp(log(a)+log((Tp-c)^b))
  develop_rate <- (1/develop_time)
  develop_df$DailyDevelopment <- (1/develop_time) #daily development rate
  
  #define development period
  
  #develop_df <- subset(develop_df, 
  #                     develop_df$Jdate >= startspawn & develop_df$Jdate <= (startspawn + 366))
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


##### Spawning Set-Up ------

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


###### Initializing -----

## Parameters
Tp <- temp
a <- exp(10.404)
b <- -2.043
c <- -7.575

## Initalizing
develop_df <- temp_df
emergence <- NULL
DailyRate <- rep(0,length(temp))


###### Egg Stage -----

for(i in 1:length(temp)) {
  
  if(i >= spawndate_lower){
    egg.rate = DevelopmentFunc(Tp = temp[i], a, b, c, startspawn = i) #daily development rate
    DailyRate[i] = egg.rate
    emergence[i] = EmergenceFunc(Tp, a, b, c, startspawn = i) #emergence day vector
    
  } else {
    egg.rate = 0
    DailyRate[i] = 0
    emergence[i] = 0
  }
}

times = temp_df$Jdate

### Emergence

if(i > spawndate_lower) {
  emergence[is.na(emergence)] <- 0
  emergence <- emergence[! emergence %in% c('0')]
    #cutting emergence vector to correct size
  emerge.day <- rep(emergence, spawn_df$SpawnFreq) 
    #repeats a given emergence day by the number of spawns on the associated spawn day
  freq.emerge <- as.vector(table(emerge.day))
    #creates a vector that lists the number of fish emerging on a given day
    #### NEEDS FIX
    #### Need to preserve the zeros from spawnfreq -> emerge freq
  
  spawn_df$EmergeDay <- emergence
  x <- cbind.data.frame(unique(emerge.day), freq.emerge); names(x)[1] <- c("EmergeDay")
  spawn_df <- merge(spawn_df, x, by = "EmergeDay", all=T) ; spawn_df[is.na(spawn_df$freq.emerge),"freq.emerge"] <- 0
  
  emerge.window <- min(emerge.day):max(emerge.day)
    #window of emergence days
  hist(emerge.day) #histogram displaying number of fish emerging on a given day
}

#plot(DailyRate ~ times)
#plot(cumsum(DailyRate) ~ times, ylim=c(0,1), xlim=c(200,500))


###### Juvenile Stage ------

SizeFunc = function(a1, b1, c1, d1, e1, j, date) {
  
  #restrict data frame to only dates after our specified emergence date
  size_df <- subset(temp_df, temp_df$Jdate >= j) 
  
  # restrict size_df to between emergence and a certain date
  size_df <- subset(size_df, size_df$Jdate <= date)
  
  emergence_weight = 35^2.953*0.00001432 #emergence mass (size 35 mm)
  omega <- (a1*(size_df$Mean-b1)*(1-exp(c1*(size_df$Mean-d1)))) #growth rate
  omega <- if_else(omega > 0, omega, 0) #growth rate can not be 'negative'
  
  juv_weight <- rep(0,(date-j))
  
  for(k in 1:(date - j + 1)) { #(k in j:date)
    if(k == 1) {
      juv_weight[k] <- (emergence_weight^(e1) + omega[k]*(1/100)*e1)^(1/e1) 
       
        #juv weight in grams
    }else{
      juv_weight[k] <- (juv_weight[k-1]^(e1) + omega[k]*(1/100)*e1)^(1/e1) 
        #juv weight in grams
    }
  }
  
  #juv_weight <- (emergence_weight^(e1) + omega*(1/100)*e1*[period])^(1/e1)
  juv_length <- ((juv_weight/0.00001432)^(1/2.953)) #juvenile length in mm
  
  #creates 'leading zeroes' until emergence day for matrix
  size_fix <- rep(0, (j - min(emerge.day)))
  juv_length1 <- c(size_fix, juv_length) 
  
  return(juv_length)
  
} 

# Parameters for Juv Size
a1 <- 0.415 #d
b1 <- 1.833 #T_L
c1 <- 0.315 #g
d1 <- 24.918 #T_U
e1 <- 0.338 #b

#initialize
juv_size <- rep(0,length(emerge.window))
j = min(emerge.day)


#Juv Size
for(j in emergence) {

  #extracts juvenile length on given data
  juv_length <- SizeFunc(a1, b1, c1, d1, e1, j, date = 562)
  
  #writes the size on the given date for a particular emergence day
  juv_size[j] <- tail(juv_length, n=1) 
  
}

#Modeled emergence is not continuous, so need to remove NAs
juv_size[is.na(juv_size)] <- 0 ; juv_size <- juv_size[! juv_size %in% c('0')]


### Graphing

if(i > min(emergence)) {
  size_dist <- rep(juv_size, spawn_df$freq.emerge) 
  #replicates size-by-emergence by the number emerged on that day
  hist(size_dist) 
  #histogram displaying number of fish emerging on a given day
}

hist(size_dist, xlim=c(50,58), xlab="Length", 
     main="  Distribution of Marsh Creek 
     Juvenile Lengths for July 16, 2016")

hist(emerge.day, xlim=c(420,525), xlab="Day of Emergence",
     main="   Distribution of Emergence Days
     for Marsh Creek Spawning of 2016")

hist(spawn_dist, xlim=c(210,250), xlab="Spawn Date", 
     main="Distribution of Spawning Date 
     for Marsh Creek, 2016")




#Dataframing
juv_size_df <- temp_df[c("Date", "Jdate")]
juv_size_df <- subset(juv_size_df, juv_size_df$Jdate >= min(emerge.day))
juv_size_df['Age'] <- row_number(juv_size_df)
juv_size_df$Length <- juv_length

##Make a Matrix
# Creates matrix with the rows as the date of emergence, columns as date
# and adds the juvenile size vector for a given emergence day on each loop
if(j == min(emerge.day)) {
  
  size_matrix <- rbind((min(emerge.day):length(temp)), juv_size)
  colnames(size_matrix) <- c(min(emerge.day):length(temp))
  rownames(size_matrix) <- c("Julian Date", j)
  
} else {
  
  size_matrix <- rbind(size_matrix, juv_size)
  size_matrix <- size_matrix[rownames(size_matrix) != 'Julian Date',]
  rownames(size_matrix)[(j-min(emerge.day)+1)] <- c(j)
  
}
size_matrix
