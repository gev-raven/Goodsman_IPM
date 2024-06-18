

#temp_df <- read.csv("C:\\Users\\Grace.Veenstra\\Documents\\GitHub\\Goodsman_IPM\\Data\\chinook data for Chittaro bioenergetics.csv", header=TRUE) 
temp_df <- read.csv("~\\GitHub\\Goodsman_IPM\\Data\\chinook data for Chittaro bioenergetics.csv", 
                    header=TRUE) 
temp_df <- temp_df[which(temp_df$stream %in% "MAR"),]
    ## Change Stream
temp_df <- temp_df[which(temp_df$year %in% "2003"),]
    ## Change Year
temp_df <- temp_df[grep(".9", temp_df$id),]
    ## Change Obs Month - .7 or .9

temp <- subset(temp_df$T, temp_df$id %in% temp_df$id[which.min(temp_df$doy)])
    ## sets temp to the time period of the earliest cohort


j = min(temp_df$doy)
  
date = max(temp_df$doy, temp_df$id %in% temp_df$id[which.max(temp_df$doy)])
    ## sets date to be the max date in the cohort range


emergence_length = 35  #(in mm)
emergence_weight = (emergence_length)^2.953*0.00001432 #(in grams)

juv_weight <- rep(0,(date-j))
omega <- (a1*(temp-b1)*(1-exp(c1*(temp-d1)))) #growth rate
omega <- if_else(omega > 0, omega, 0)

for(k in 1:(date - j + 1)) { #(k in j:date)
  if(k == 1) {
    juv_weight[k] <- (emergence_weight^(e1) + omega[k]*(1/100)*e1)^(1/e1) 
    #juv weight in grams
  }else{
    juv_weight[k] <- (juv_weight[k-1]^(e1) + omega[k]*(1/100)*e1)^(1/e1)
  }
}

juv_length <- ((juv_weight/0.00001432)^(1/2.953)) #juvenile length in mm




#2003.9 -----
plot(juv_weight, type='l', lty=1, lwd=2, ylim=c(0,6), xlab="Days since 'emergence'", ylab="Juvenile Weight", main="Marsh Creek, 2003 - Juvenile Weight")
legend("bottomright", legend=c("Modeled","Otolith"),lty=1,col=c(1,2),bty='n')
lines(temp_df$mass[which(temp_df$id %in% "MAR.2003.9.C1")], col="red")
lines(temp_df$mass[which(temp_df$id %in% "MAR.2003.9.C2")], col="red")
lines(temp_df$mass[which(temp_df$id %in% "MAR.2003.9.C3")], col="red")
lines(temp_df$mass[which(temp_df$id %in% "MAR.2003.9.C4")], col="red")
lines(temp_df$mass[which(temp_df$id %in% "MAR.2003.9.C5")], col="red")

plot(juv_length, type='l', lty=1, lwd=2, ylim=c(35,100), xlab="Days since 'emergence'", ylab="Juvenile Length", main="Marsh Creek, 2003 - Juvenile Length")
legend("bottomright", legend=c("Modeled","Otolith"),lty=1,col=c(1,2),bty='n')
lines(temp_df$fl[which(temp_df$id %in% "MAR.2003.9.C1")], col="red")
lines(temp_df$fl[which(temp_df$id %in% "MAR.2003.9.C2")], col="red")
lines(temp_df$fl[which(temp_df$id %in% "MAR.2003.9.C3")], col="red")
lines(temp_df$fl[which(temp_df$id %in% "MAR.2003.9.C4")], col="red")
lines(temp_df$fl[which(temp_df$id %in% "MAR.2003.9.C5")], col="red")

#2004 -----
plot(juv_weight, type='l', lty=1, lwd=2, ylim=c(0,6), xlab="Days since 'emergence'", ylab="Juvenile Weight", main="Marsh Creek, 2004 - Juvenile Weight")
legend("bottomright", legend=c("Modeled","Otolith"),lty=1,col=c(1,2),bty='n')
lines(temp_df$mass[which(temp_df$id %in% "MAR.2004.9.C1")], col="red")
lines(temp_df$mass[which(temp_df$id %in% "MAR.2004.9.C2")], col="red")
lines(temp_df$mass[which(temp_df$id %in% "MAR.2004.9.C3")], col="red")
lines(temp_df$mass[which(temp_df$id %in% "MAR.2004.9.C4")], col="red")
lines(temp_df$mass[which(temp_df$id %in% "MAR.2004.9.C5")], col="red")

plot(juv_length, type='l', lty=1, lwd=2, ylim=c(35,100), xlab="Days since 'emergence'", ylab="Juvenile Length", main="Marsh Creek, 2004 - Juvenile Length")
legend("bottomright", legend=c("Modeled","Otolith"),lty=1,col=c(1,2),bty='n')
lines(temp_df$fl[which(temp_df$id %in% "MAR.2004.9.C1")], col="red")
lines(temp_df$fl[which(temp_df$id %in% "MAR.2004.9.C2")], col="red")
lines(temp_df$fl[which(temp_df$id %in% "MAR.2004.9.C3")], col="red")
lines(temp_df$fl[which(temp_df$id %in% "MAR.2004.9.C4")], col="red")
lines(temp_df$fl[which(temp_df$id %in% "MAR.2004.9.C5")], col="red")


#2006.7 -----
plot(juv_weight, type='l', lty=1, lwd=2, ylim=c(0,6), xlab="Days since 'emergence'", ylab="Juvenile Weight", main="Marsh Creek, 2006 - Juvenile Weight")
legend("bottomright", legend=c("Modeled","Otolith"),lty=1,col=c(1,2),bty='n')
lines(temp_df$mass[which(temp_df$id %in% "MAR.2006.7.C1")], col="red")
lines(temp_df$mass[which(temp_df$id %in% "MAR.2006.7.C2")], col="red")
lines(temp_df$mass[which(temp_df$id %in% "MAR.2006.7.C3")], col="red")
lines(temp_df$mass[which(temp_df$id %in% "MAR.2006.7.C4")], col="red")
lines(temp_df$mass[which(temp_df$id %in% "MAR.2006.7.C5")], col="red")

#2006.9 ------
lines(juv_weight, )
lines(temp_df$mass[which(temp_df$id %in% "MAR.2006.9.C1")], col="red")
lines(temp_df$mass[which(temp_df$id %in% "MAR.2006.9.C2")], col="red")
lines(temp_df$mass[which(temp_df$id %in% "MAR.2006.9.C3")], col="red")






#2011 ------
plot(juv_weight, type='l', lty=1, lwd=2, ylim=c(0,6), xlab="Days since 'emergence'", ylab="Juvenile Weight", main="Marsh Creek, 2011 - Juvenile Weight")
legend("bottomright", legend=c("Modeled","Otolith"),lty=1,col=c(1,2),bty='n')
lines(temp_df$mass[which(temp_df$id %in% "MAR.2011.7.C1")], col="red")
lines(temp_df$mass[which(temp_df$id %in% "MAR.2011.7.C2")], col="red")
lines(temp_df$mass[which(temp_df$id %in% "MAR.2011.7.C3")], col="red")
lines(temp_df$mass[which(temp_df$id %in% "MAR.2011.7.C4")], col="red")
lines(temp_df$mass[which(temp_df$id %in% "MAR.2011.7.C5")], col="red")
lines(temp_df$mass[which(temp_df$id %in% "MAR.2011.7.C6")], col="red")
lines(temp_df$mass[which(temp_df$id %in% "MAR.2011.7.C7")], col="red")
lines(temp_df$mass[which(temp_df$id %in% "MAR.2011.7.C8")], col="red")
lines(temp_df$mass[which(temp_df$id %in% "MAR.2011.7.C9")], col="red")
lines(temp_df$mass[which(temp_df$id %in% "MAR.2011.7.C10")], col="red")
lines(temp_df$mass[which(temp_df$id %in% "MAR.2011.7.C11")], col="red")
lines(temp_df$mass[which(temp_df$id %in% "MAR.2011.7.C12")], col="red")
lines(temp_df$mass[which(temp_df$id %in% "MAR.2011.7.C13")], col="red")

temp_df$length <- ((temp_df$mass/0.00001432)^(1/2.953))
plot(juv_length, type='l', lty=1, lwd=2, ylim=c(35,70), xlab="Days since 'emergence'", ylab="Juvenile Length", main="Marsh Creek, 2011 - Juvenile Length")
legend("bottomright", legend=c("Modeled","Otolith"),lty=1,col=c(1,2),bty='n')
lines(temp_df$fl[which(temp_df$id %in% "MAR.2011.7.C1")], col="red")
lines(temp_df$fl[which(temp_df$id %in% "MAR.2011.7.C2")], col="red")
lines(temp_df$fl[which(temp_df$id %in% "MAR.2011.7.C3")], col="red")
lines(temp_df$fl[which(temp_df$id %in% "MAR.2011.7.C4")], col="red")
lines(temp_df$fl[which(temp_df$id %in% "MAR.2011.7.C5")], col="red")
lines(temp_df$fl[which(temp_df$id %in% "MAR.2011.7.C6")], col="red")
lines(temp_df$fl[which(temp_df$id %in% "MAR.2011.7.C7")], col="red")
lines(temp_df$fl[which(temp_df$id %in% "MAR.2011.7.C8")], col="red")
lines(temp_df$fl[which(temp_df$id %in% "MAR.2011.7.C9")], col="red")
lines(temp_df$fl[which(temp_df$id %in% "MAR.2011.7.C10")], col="red")
lines(temp_df$fl[which(temp_df$id %in% "MAR.2011.7.C11")], col="red")
lines(temp_df$fl[which(temp_df$id %in% "MAR.2011.7.C12")], col="red")
lines(temp_df$fl[which(temp_df$id %in% "MAR.2011.7.C13")], col="red")


