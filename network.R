##################################################################################

#NETWORK DATA

##################################################################################
#INIT

#install.packages("dplyr")
#install.packages("sqldf")
#install.packages("repmis")
#install.packages("RCurl")
library(dplyr)
library(sqldf)
library(repmis)Y
library(RCurl)

#set graph engine default 
old.par <- par(no.readonly = TRUE)

source_data("https://github.com/mwkoomen/network_analysis/blob/master/Network_Data.rdata?raw=true")
script <- ("https://github.com/mwkoomen/network_analysis/raw/master/network.R")

###################################################################################
#MAIN

#clear workspace
rm(list=setdiff(ls(), "old.par"))

#reset graph engine 
par(old.par)

#partial set 
par_set <- `Network Data`

#missing connections for full set
comp_set <- expand.grid(1991:2017, 1:200, 1:2, 1:200, 1:2)
dif <- sqldf("
    select 
      var1 as Year,var2 as Region_A, var3 as Country_A, var4 as Region_B, var5 as Country_B
    from comp_set 
      except 
        select 
          Year, Region_A, Country_A, Region_B, Country_B
        from par_set")
nrow(comp_set)-nrow(par_set)==nrow(dif)
full_set <- sqldf("
    select * from par_set 
      union all 
        select *, 0 from dif")
test_unique <- sqldf("
                     select count(*) as rownumb 
                     from full_set 
                     group by Year, Region_A, Country_A, Region_B, Country_B
                     order by rownumb
                     ")
max(test_unique$rownumb)==1
rm("comp_set", "dif", "test_unique")

#Normalized Intensity
par_set$Intensity_norm = (
  par_set$Intensity-min(par_set$Intensity)) / 
  (max(par_set$Intensity)-min(par_set$Intensity)
  )
  cor(par_set$Intensity, par_set$Intensity_norm) == 1
full_set$Intensity_norm = (
    full_set$Intensity-min(full_set$Intensity)) / 
    (max(full_set$Intensity)-min(full_set$Intensity)
    )
  cor(full_set$Intensity, full_set$Intensity_norm) == 1   

#partial set plot
  pc1c1 <- par_set %>%
    filter(Country_A==1, Country_B==1) %>%
    group_by(Year) %>%
    summarise(mean(Intensity))
  pc1c2 <- par_set %>%
    filter(Country_A==1, Country_B==2) %>%
    group_by(Year) %>%
    summarise(mean(Intensity))
  pc2c2 <- par_set %>%
    filter(Country_A==2, Country_B==2) %>%
    group_by(Year) %>%
    summarise(mean(Intensity))
  pc2c1 <- par_set %>%
    filter(Country_A==2, Country_B==1) %>%
    group_by(Year) %>%
    summarise(mean(Intensity))
  #plot
  plot(pc1c1,type = "o",col = "red", xlab = "Year", ylab = "Intensity", 
       main = "Network Intensity (partial)") 
  lines(pc1c2, type = "o", col = "blue") 
  lines(pc2c2, type = "o", col = "darkgreen")
  lines(pc2c1, type = "o", col = "orange")
  legend(1991, 45000, legend=c("C1:Internal", "C1:External", "C2:Internal", "C2:External"),
         col=c("red", "blue", "darkgreen", "orange"), lty=1:2, cex=0.5)
  
#full set plot
fc1c1 <- full_set %>%
  filter(Country_A==1, Country_B==1) %>%
  group_by(Year) %>%
  summarise(mean(Intensity))
fc1c2 <- full_set %>%
  filter(Country_A==1, Country_B==2) %>%
  group_by(Year) %>%
  summarise(mean(Intensity))
fc2c2 <- full_set %>%
  filter(Country_A==2, Country_B==2) %>%
  group_by(Year) %>%
  summarise(mean(Intensity))
fc2c1 <- full_set %>%
  filter(Country_A==2, Country_B==1) %>%
  group_by(Year) %>%
  summarise(mean(Intensity))
  #plot
  plot(fc1c1,type = "o",col = "red", xlab = "Year", ylab = "Intensity", 
       main = "Network Intensity (full)") 
  lines(fc1c2, type = "o", col = "blue") 
  lines(fc2c2, type = "o", col = "darkgreen")
  lines(fc2c1, type = "o", col = "orange")
  legend(1991, 8500, legend=c("C1:Internal", "C1:External", "C2:Internal", "C2:External"),
         col=c("red", "blue", "darkgreen", "orange"), lty=1:2, cex=0.6)






