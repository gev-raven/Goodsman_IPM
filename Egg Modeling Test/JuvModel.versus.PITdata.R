

###### Set your files, pick your regime ~ -----

#Stream_1 is what the fish start from (etc. spawn in Marsh, start from Marsh)
#Stream_2 is a supplemental stream for addition

TempSwitchDate = "Oct 1 2017" #date switch temp. data; will need to comment in/out relevant piece in for-loop

TempFile_AimeeRegional_1 = "C:\\Users\\Grace.Veenstra\\Documents\\GitHub\\Goodsman_IPM\\Data\\UpperMFSalmon_TempModel_EntireTimeSeries.csv"
HUC10_1 = "X1706020503"
StreamName_1 = "Marsh Creek"
StreamNameAbb_1 = "MAR"

TempFile_AimeeRegional_2 = "C:\\Users\\Grace.Veenstra\\Documents\\GitHub\\Goodsman_IPM\\Data\\LowerSalmon_TempModel_EntireTimeSeries.csv"
HUC10_2 = "X1706020906"
StreamName_2 = "White Bird Creek"
StreamNameAbb_2 = "MAR"

TempFile_MARSS_1 = "C:\\Users\\Grace.Veenstra\\Documents\\GitHub\\Goodsman_IPM\\Data\\Marsh_hourly_1993to2020.csv"
StreamName_1 = "Marsh Creek"
TempFile_MARSS_2 = "C:\\Users\\Grace.Veenstra\\Downloads\\whitebird_env2.csv"
StreamName_2 = "White Bird"

GrowthFile = "C:\\Users\\Grace.Veenstra\\Documents\\GitHub\\Goodsman_IPM\\Data\\Wild Chinook recaps_growth data all years.csv"
GrowthStream = "Marsh Creek"

StartYear = "2017"
#Choose ModelType and ModelStreams
TempModelType = "Aimee"
  TempModelType = "MARSS" 
TempModelStreams = "(Marsh)"
  TempModelStreams = "(Marsh/WhiteBird)"

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


##### TEMP CODING [USE ONLY RELEVANT SECTION] -----
###########################
  
#### TEMP - AIMEE (All Annual)
#temp_df <- read.csv("C:\\Users\\Grace.Veenstra\\Documents\\GitHub\\Goodsman_IPM\\Data\\TempModel_2003-2004.csv", header=TRUE) 
temp_df <- read.csv("~\\GitHub\\Goodsman_IPM\\Data\\TempModel_2011-2012.csv", 
                    header=TRUE)
temp_df <- temp_df[c("Date","HUC_10","Mean")] 
temp_df <- temp_df[which(temp_df$HUC_10 %in% "1706020503"),]
      ## Change Stream Identifier
temp_df$Date <- parse_date_time(temp_df$Date, orders=c('mdy','ymd')) 
temp_df$Date <- format.Date(temp_df$Date, format="%Y-%m-%d")
temp_df$Ydate <- yday(temp_df$Date)
temp_df$Jdate <- 1:nrow(temp_df)


#### TEMP 2 - AIMEE DATA (Regional - Entire Time Series)
#Temp Set #1
temp_df <- read.csv(paste(TempFile_AimeeRegional_1), header=TRUE)
temp_df <- temp_df[c("Date", paste(HUC10_1))] 
                            #Marsh Creek
colnames(temp_df) <- c("Date", "Mean")
temp_df$Date <- parse_date_time(temp_df$Date, orders=c('mdy','ymd')) 
temp_df$Date <- format.Date(temp_df$Date, format="%Y-%m-%d")
temp_df$Ydate <- yday(temp_df$Date)

#Temp Set #2
temp_df2 <- read.csv(paste(TempFile_AimeeRegional_2), header=TRUE)
temp_df2 <- temp_df2[c("Date", paste(HUC10_2))] 
                              #White Bird Creek
colnames(temp_df2) <- c("Date", "Mean")
temp_df2$Date <- parse_date_time(temp_df2$Date, orders=c('mdy','ymd')) 
temp_df2$Date <- format.Date(temp_df2$Date, format="%Y-%m-%d")
temp_df2$Ydate <- yday(temp_df2$Date)


#### TEMP 3 - MARSS DATA (Empirical)
temp_df <- read.csv(paste(TempFile_MARSS_1), header=TRUE) # SPECIFIC TO MARSH CREEK HORULY
temp_df <- temp_df[c("Observe.date", "Temperature")]
colnames(temp_df) <- c("Date", "Temp")
temp_df$Date <- parse_date_time(temp_df$Date, orders=c('mdy HM','ymd HM')) 
temp_df$Date <- format.Date(temp_df$Date, format="%Y-%m-%d")
temp_df <- temp_df %>% group_by(Date) %>% 
  mutate(Mean = (max(Temp)+min(Temp))/2) %>% 
  distinct(Date, .keep_all=TRUE)
temp_df$Ydate <- yday(temp_df$Date)

temp_df2 <- read.csv(paste(TempFile_MARSS_2), header=TRUE) #SPECIFIC TO WHITEBIRD CSV
temp_df2 <- temp_df2[c("Date", "TMEAN")]
colnames(temp_df2) <- c("Date", "Mean")
temp_df2$Date <- parse_date_time(temp_df2$Date, orders=c('mdy','ymd')) 
temp_df2$Date <- format.Date(temp_df2$Date, format="%Y-%m-%d")
#temp_df2 <- temp_df2 %>% group_by(Date) %>% 
#  mutate(Mean = (max(Temp)+min(Temp))/2) %>% 
#  distinct(Date, .keep_all=TRUE)
temp_df2$Ydate <- yday(temp_df2$Date)


#### GROWTH / RECAPTURE DATA -----
growth_df <- read.csv(paste(GrowthFile), header=TRUE) #reads csv files
growth_df <- growth_df[which(growth_df$Release.Site.Name %in% GrowthStream),] #selects stream
                            ## Change Stream Name
growth_df$Mark.Date <- parse_date_time(growth_df$Mark.Date, orders=c('mdy','ymd')) 
  growth_df$Mark.Date <- format.Date(growth_df$Mark.Date, format="%Y-%m-%d")
growth_df$Recap.Date <- parse_date_time(growth_df$Recap.Date, orders=c('mdy','ymd')) 
  growth_df$Recap.Date <- format.Date(growth_df$Recap.Date, format="%Y-%m-%d")
growth_df <- growth_df[grep(paste(StartYear), growth_df$Mark.Date),] #selects start year
                            ## Change Year


###################### -----

# Parameters for Juv Size
a1 <- 0.415 #d
b1 <- 1.833 #T_L
c1 <- 0.315 #g
d1 <- 24.918 #T_U
e1 <- 0.338 #b

model.recap.weight <- rep(0,length(growth_df$Tag.Code))
model.recap.length <- rep(0,length(growth_df$Tag.Code))

for(i in 1:nrow(growth_df)) {
  
  ind.fish <- growth_df[i,] #creates dataframe for just single fish we're running for
  j = ind.fish$Mark.DOY
  r = ind.fish$Recap.DOY + 365
  start_weight = ind.fish$Mark.Weight.g
  #date <- subset(temp_df$Jdate, temp_df$Jdate > j)
  
  #If DOING TEMP DIVIDE, USE BELOW
  divide.date <- parse_date_time(paste(TempSwitchDate), c("mdy", "ydm"))
  Tp1 <- subset(temp_df$Mean, temp_df$Date >= ind.fish$Mark.Date & temp_df$Date < divide.date)
  Tp2 <- subset(temp_df2$Mean, temp_df2$Date >= divide.date & temp_df2$Date <= ind.fish$Recap.Date)
  temp <- c(Tp1, Tp2)
  
  #IF USING ONE TEMP FILE, USE BELOW
  #temp <- subset(temp_df$Mean, temp_df$Date >= ind.fish$Mark.Date & temp_df$Date <= ind.fish$Recap.Date)
  
  juv_weight <- rep(0,(r-j))
  omega <- (a1*(temp-b1)*(1-exp(c1*(temp-d1)))) #growth rate
  omega <- if_else(omega > 0, omega, 0) #growth rate cannot be less than zero
  
  for(k in 1:(r - j + 1)) { #(k in j:date)
    if(k == 1) {
      juv_weight[k] <- (start_weight^(e1) + omega[k]*(1/100)*e1)^(1/e1) #weight (grams)
    } else{
      juv_weight[k] <- (juv_weight[k-1]^(e1) + omega[k]*(1/100)*e1)^(1/e1) #accumulates with each time step
    }
  }
  
  juv_length <- ((juv_weight/0.00001432)^(1/2.953)) #juvenile length (in mm)
  
  model.recap.weight[i] <- tail(juv_weight, n=1) #model-predicated recapture weight on recapture day
  model.recap.length[i] <- tail(juv_length, n=1) #model-predicated recapture length
  
}

  
##GRAPHING

#weight
plot(model.recap.weight, type='l', mgp=c(2,.6,0), las=1,
     ylim=c(min(model.recap.weight, growth_df$Recap.Weight.g), 
            max(model.recap.weight, growth_df$Recap.Weight.g)),
     xlab="Individual Fish", ylab="Weight (g)", 
     main="Comparison of Model Weight versus Actual Weight", 
     sub =paste(StreamName_1,StartYear,TempModelType,TempModelStreams))
legend("topleft", legend=c("Modeled","Recapture"),lty=1,lwd=2,col=c(1,2),bty='n')
lines(growth_df$Recap.Weight.g, col="red")
x.axis <- seq(1,length(growth_df$Tag.Code),1)
polygon(c(x.axis, rev(x.axis)), c(model.recap.weight, rev(growth_df$Recap.Weight.g)),
        density=10, angle=45)
lines(growth_df$Recap.Weight.g, lwd=2, col="red")
lines(model.recap.weight, lwd = 2)

#length
plot(model.recap.length, type='l', mgp=c(2,.6,0), las=1,
     ylim=c(min(model.recap.length, growth_df$Recap.Length.mm), 
            max(model.recap.length, growth_df$Recap.Length.mm)),
     xlab="Individual Fish", ylab="Length (mm)", 
     main="Comparison of Model Length versus Actual Length", 
     sub = paste(StreamName_1,StartYear,TempModelType,TempModelStreams))
legend("topleft", legend=c("Modeled","Recapture"),lty=1,lwd=2,col=c(1,2),bty='n')
lines(growth_df$Recap.Length.mm, col="red")
x.axis <- seq(1,length(growth_df$Tag.Code),1)
polygon(c(x.axis, rev(x.axis)), c(model.recap.length, rev(growth_df$Recap.Length.mm)),
        density=10, angle=45)
lines(growth_df$Recap.Length.mm, lwd=2, col="red")
lines(model.recap.length, lwd = 2)


#predicted v observed
par(mfrow=c(2,2))

plot(model.recap.weight, growth_df$Recap.Weight.g,
     mgp=c(2,.6,0), las=1,
     xlab="Predicted Weight (Model)", ylab="Observed Weight (Recapture)",
     main="Predicted vs. Observed Weight",
     sub=paste(StreamName_1,StartYear,TempModelType,TempModelStreams))
abline(a=0,b=1)

plot(sort(model.recap.weight), sort(growth_df$Recap.Weight.g),
     mgp=c(2,.6,0), las=1,
     xlab="Predicted Weight (Model)", ylab="Observed Weight (Recapture)",
     main="Quantile Predicted vs. Observed Weight",
     sub=paste(StreamName_1,StartYear,TempModelType,TempModelStreams))
abline(a=0,b=1)

plot(model.recap.length, growth_df$Recap.Length.mm,
     mgp=c(2,.6,0), las=1,
     xlab="Predicted Length (Model)", ylab="Observed Length (Recapture)",
     main="Predicted vs. Observed Length",
     sub=paste(StreamName_1,StartYear,TempModelType,TempModelStreams))
abline(a=0,b=1)

plot(sort(model.recap.length), sort(growth_df$Recap.Length.mm),
     mgp=c(2,.6,0), las=1,
     xlab="Predicted Length (Model)", ylab="Observed Length (Recapture)",
     main="Quantile Predicted vs. Observed Length",
     sub=paste(StreamName_1,StartYear,TempModelType,TempModelStreams))
abline(a=0,b=1)


