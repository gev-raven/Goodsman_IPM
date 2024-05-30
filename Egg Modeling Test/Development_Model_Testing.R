## Development_Model_Testing

# This code creates a distribution of egg emergence days for Chinook salmon
# based on the spawn date, with variance applied to the spawn date.
# The development rate equation comes from Abbie, and from the Beacham and Murray Model

# Temperature data is for Upper Middle Fork Salmon River at Marsh Creek.
# Model assumes no egg mortality.


##############################################
###### Packages Set Up ----- 
##############################################

install.packages("lubridate")
library(lubridate)
install.packages("anytime")
library(anytime)
install.packages("tidyverse")
library(tidyverse)
install.packages("dbplyr")
library(dbplyr)
install.packages("dplyr")
library(dplyr)
install.packages("dtplyr")
library(dtplyr)


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


############################################
###### Development model -----
############################################

development_func = function(Tp, a, b, c, startspawn) {
  
  # develop_time is the development time (in days)
  # Tp is daily mean temperature
  
  develop_df <- temp_df
  develop_time <- exp(log(a)+log((Tp-c)^b))
  develop_df$DailyDevelopment <- (1/develop_time) #daily development rate
  
  #define development period
  
  develop_df <- subset(develop_df, 
                       develop_df$Jdate >= startspawn & develop_df$Jdate <= (startspawn + 366))
  ## restricts develop_df to the period of development
  
  develop_df$TotalDevelopment <- cumsum(develop_df$DailyDevelopment)
  ## total development = sum of daily development
  
  y <- develop_df[min(which((develop_df$TotalDevelopment) >= 1)), "Jdate"] 
  ## once development = 1, it is emergence time and we extract emergence date
  
  return(y)
}

## Parameters
Tp <- temp
a <- exp(10.404)
b <- -2.043
c <- -7.575

## Spawning Distribution
spawn_mean <- 230 #aug 18; mean spawn start date
spawn_sd <- 5 #spawning date standard deviation
n.spawners <- 10000 #number of spawning fish
spawn_dist <- rnorm(n = n.spawners, mean = spawn_mean, sd = spawn_sd) 
            #creates normal distribution of spawners
spawn_dist <- round(spawn_dist)
hist(spawn_dist)

## Spawn Window
spawndate_lower <- min(spawn_dist) #lower limit of spawn window
spawndate_upper <- max(spawn_dist) #upper limit of spawn window
spawn.window <- spawndate_lower:spawndate_upper

## Initalizing
develop_df <- temp_df
emergence <- NULL
i = 1

## Iteration over Spawn Range

for(i in 1:length(spawn.window)) { #runs loop for each day in spawning window
  
  startspawn <- spawn.window[i] #set spawning start as day in spawn window
  
  emergence[i] <- development_func(Tp = Tp, a = a, b = b, c = c, startspawn = startspawn)
    #collects the emergence days into a vector
    
}


## Emergence
freq.spawn <- as.vector(table(spawn_dist)) 
    #turns spawning distribution into vector of number of spawns per day
freq.spawn

emerge.day <- rep(emergence, freq.spawn) 
    #repeats a given emergence day by the number of spawns on the associated spawn day
freq.emerge <- as.vector(table(emerge.day))
    #creates a vector that lists the number of fish emerging on a given day
emerge.window <- min(emerge.day):max(emerge.day)
    #window of emergence days
hist(emerge.day) #histogram displaying number of fish emerging on a given day

adults <- n.spawners - cumsum(freq.spawn) #looks at decline of spawners
plot(adults ~ spawn.window)
eggs <- 10000 - cumsum(freq.emerge) #tracks the decline of 'eggs' as they transition to 'emerged'
plot(eggs ~ emerge.window)
plot(eggs)


## qqnorm

qqnorm(emerge.day)
qqline(emerge.day)


## Graphing

spawn_df <- temp_df
spawn_df <- subset(temp_df, min(spawn_dist) <= temp_df$Jdate 
                   & temp_df$Jdate <= max(emergence)) #attempting a development window

times = 1:temp_df$Jdate

plot(develop_df$TotalDevelopment ~ develop_df$Jdate, ylim=c(0,1))
plot(develop_df$DailyDevelopment ~ develop_df$Mean)

plot(emergence ~ temp)
plot(emergence ~ spawn_df$Mean)

