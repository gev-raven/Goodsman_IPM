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
temp_df <- read.csv("Egg Modeling Test\\MarshCreek_Temp_2015-2016.csv",
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
  
  develop_df <- subset(develop_df, 
                       develop_df$Jdate >= startspawn & develop_df$Jdate <= (startspawn + 366))
  ## restricts develop_df to the period of development
  
  develop_df$TotalDevelopment <- cumsum(develop_df$DailyDevelopment)
  ## total development = sum of daily development
  
  y <- develop_df[min(which((develop_df$TotalDevelopment) >= 1)), "Jdate"] 
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
  emergence2 <- emergence[! emergence %in% c('0')]
    #cutting emergence vector to correct size
  emerge.day <- rep(emergence2, spawn_df$SpawnFreq) 
    #repeats a given emergence day by the number of spawns on the associated spawn day
  freq.emerge <- as.vector(table(emerge.day))
    #creates a vector that lists the number of fish emerging on a given day
  emerge.window <- min(emerge.day):max(emerge.day)
    #window of emergence days
  hist(emerge.day) #histogram displaying number of fish emerging on a given day
}

plot(DailyRate ~ times)
plot(cumsum(DailyRate) ~ times, ylim=c(0,1), xlim=c(200,500))


###### Juvenile Stage ------








##### irrelevant -----

###Clean temperature data

#ELK.temp$Site.name <- recode(ELK.temp$Site.name, "Bear Valley/Elk Creek" = "ELK")

#temp <- rbind(BVA.temp, CHO.temp, LAK.temp, MAR.temp, SFS.temp, ELK.temp) %>% select(-Water.depth)
#temp$Observe.date <- str_sub(temp$Observe.date, end = -10)

#VAL.temp$Observe.date <- str_sub(VAL.temp$Observe.date, end = -6)
#VAL.temp <- select(VAL.temp, -Water.depth)

#temp <- rbind(temp, VAL.temp)
#temp$Observe.date <- as.Date(temp$Observe.date, format = "%m/%d/%Y")
#temp$year <- str_sub(temp$Observe.date, end = -7)
#temp$doy <- strftime(temp$Observe.date, format = "%j")
#temp$doy <- as.numeric(temp$doy)
#temp <- rename_(temp, "stream" = "Site.name")
#temp$stream <- recode(temp$stream,
                      #"Bear Valley/Elk Creek" = "BVA",
                      #"Cape Horn Creek" = "CHO",
                      #"Lake Creek" = "LAK",
                      #"Marsh Creek" = "MAR",
                      #"South Fork Salmon" = "SFS",
                      #"Valley Creek" = "VAL")

###Calculate a daily avg temp by date and stream
#avg.temp <- temp %>%
#  group_by(Observe.date, stream) %>%
#  mutate(avg.daily.temp = mean(Temperature, na.rm = TRUE)) %>%
#  ungroup() %>%
#  distinct(Observe.date, .keep_all = TRUE) %>%
#  select(-Temperature) %>%
#  filter(doy >= 89 & doy <= 214)


#calculate size after once month, two months, and at August 1st using cumsum of daily development

#dailymodel<-avg.temp%>%
#  group_by(stream,year) %>%
#  mutate(totalgrowth214 = 
#             lag(rollapply(dailygrowth, (214-doy), sum, align = "left", fill=NA)))%>%
#  mutate(totalgrowth1month = 
#           lag(rollapply(dailygrowth, (30), sum, align = "left", fill=NA))) %>%
#  mutate(totalgrowth2month = 
#           lag(rollapply(dailygrowth, (60), sum, align = "left", fill=NA)))


#ggplot(data = filter(dailymodel, stream == "BVA", doy<150),  aes(x=doy)) + aes(color=year) +
#  geom_line(aes(y=totalgrowthmonth)) + 
#  geom_line(aes(y=totalgrowth2month)) +
#  geom_line(aes(y=totalgrowth214)) +
#  ggtitle("Size after one month for Different Emergence DOY and Years Bear Valley Creek") +
#  xlab("Day of Year") + ylab("Size")

#ggplot(data = filter(dailymodel, stream=="BVA", doy<150),  aes(x=doy)) + aes(color=year) +
#  geom_line(aes(y=totalgrowth300)) +
#  ggtitle("Size at DOY 300 for Different Emergence DOY and Years Bear Valley Creek") +
#  xlab("Day of Year") + ylab("Size")


#calculate average temperature over a time period of 30 days, 60 days, and at August 1st

#new <- avg.temp %>%
#  group_by(stream,year) %>%
#  mutate(averagetempfornext30days = 
#           lag(rollapply(avg.daily.temp, (30), mean, align = "left", fill=NA))) %>%
#  mutate(averagetempfornext60days = 
#           lag(rollapply(avg.daily.temp, (60), mean, align = "left", fill=NA))) %>%
#  mutate(averagetempuntilAugust1st = 
#           lag(rollapply(avg.daily.temp, (215-doy), mean, align = "left", fill=NA))) %>%
#  distinct(Observe.date, .keep_all = TRUE)


##### Temp stuff I want -------

SizeFunc = function(Tp, a1, b1, c1, d1, e1, emergence2) {
  
  tem_df <- subset(temp_df, temp_df$Jdate >= min(emerge.day)) 
  #restrict dataframe to only dates post-emergence
  size_df <- subset(tem_df, size_df$Jdate <= date)
  #size_df where we only want temp. dates between emergence and a given date
  X <- size_df$Mean
  
  size_at_emergence = 35^2.953*0.00001432
  size_at_date <- size_at_emergence + (a1*(X-b1)*(1-exp(c1*(X-d1))*1^(e1)))*(date-emergence2[j])
  
  juv_size[date] <- size_at_date
  
  return(juv_size)
  
} 

Parameters:
a1 <- 0.00415
b1 <- 1.833
c1 <- 0.315
d1 <- 24.918
e1 <- -0.338

size_matrix <- matrix(, nrow=length(emerge.window), ncol=c(1,lngth(temp_df$Jdate)-min(emerge.day))) %>%
      row.names(size_matrix) <- min(emerge.day):max(emerge.day) 
      col.names(size_matrix) <- 1:length(temp_df$Jdate)-min(emerge.day)) 
      # Size Matrix where rows are emergence day and columns are days since emergence
      # Value in a given cell is the size of fish on a given day

j = min(emerge.day)
date = min(emerge.day) 

if(j %in% emerge.window) {
  
  for(date in length(temp_df)) {

    juv_size[date] = juvSizeFunc(Tp, a1, b1, c1, d1, e1, emergence2)

    
    
  }
  
  size_matrix <- replace()
  
  #rbind
  #have trailing / leading NAs, how to keep and stick?
  #need to merge based on specific date...
  
}


#Turn average temperatures over a time period into size of fish, after 30 days, 60 days, and on August 1st 
temp_postemergence$sizeatemergence <- 35^2.953*0.00001432
temp_postemergence$size30dayslater = 
  temp_postemergence$sizeatemergence + (0.00415*(T-1.833)*(1-exp(0.315*(T-24.918))*1^(-0.338)))*30
temp_postemergence$size60dayslater = 
  temp_postemergence$sizeatemergence + (0.00415*(w-1.833)*(1-exp(0.315*(w-24.918))*1^(-0.338)))*60
temp_postemergence$sizeAugust1st = 
  temp_postemergence$sizeatemergence + (0.00415*(v-1.833)*(1-exp(0.315*(v-24.918))*1^(-0.338)))*(214-temp_postemergence$doy)

#ggplot(data=filter(temperaturespostemergence, stream=="BVA", doy<151),  aes(x=doy))+aes(color=year)+
#  geom_line(aes(y=sizeAugust1st))+
#  geom_line(aes(y=size30dayslater))+
#  ggtitle("Beacham Daily Development Bear Valley Creek") +
#  xlab("Day of emergence") + ylab("Size")

