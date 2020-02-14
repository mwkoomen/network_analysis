
#NETWORK DATA


# Initialize --------------------------------------------------------------


#Load block, if these packages are not installed, you have to remove the pound sign and 
#run each code. 

#install.packages("dplyr")
#install.packages("sqldf")
#install.packages("repmis")
#install.packages("RCurl")
#install.packages("igraph")
#install.packages("readstata13")

#Then load the packages into the current working session:
library(igraph)
library(dplyr)
library(sqldf)
library(readstata13)

#!!!!! Run seperately the first time: this might prompt you with a quesiton if you want to load 
#external data in a temporary cache or create a seperate folder. Either is fine, 
#but you have to press either Y/n to continue. 
library(repmis)

library(RCurl)

#Code that saves the default settings of the graph engine
old.par <- par(no.readonly = TRUE)


# Load data ---------------------------------------------------------------

#clear workspace and reset graph engine 
rm(list=setdiff(ls(), "old.par"))
par(old.par)

#This imports the Rdata file from github 
#The data is stored as the partial set of data, i.e. null connections are missing
source_data("https://github.com/mwkoomen/network_analysis/blob/master/Network%20Data.rdata?raw=true")
par_set <- `Network Data`

#If the import doesn't work, you can load the data locally (from the original stata file): 
#your_work_directory <- "C:/Users/User/Documents/network_analysis"
#setwd(your_work_directory)
#par_set <- read.dta13("network_data.dta")


# Normalize Intensity -----------------------------------------------------

#This code creates a normalized measure for the link intensity and test the 
#correlation with the raw version of the link intensity. Should return "TRUE"
par_set$Intensity_norm = (
  par_set$Intensity-min(par_set$Intensity)) / 
  (max(par_set$Intensity)-min(par_set$Intensity)
  )
cor(par_set$Intensity, par_set$Intensity_norm) == 1


# Create full set ---------------------------------------------------------

#This block creates a full data set that contains all the missing links 

#This code creates a full (theoretical) set of dimensions 27*200*2*200*2
comp_set <- expand.grid(1991:2017, 1:200, 1:2, 1:200, 1:2)

#This code compares the (theorectical) full set to the partial (real) data set
#The result is a data frame that contains all the combinations of years and regions per country 
#that are not contained in the partial data set
dif <- sqldf("
              select 
                var1 as Year,var2 as Region_A, var3 as Country_A, 
                var4 as Region_B, var5 as Country_B
              from comp_set 
              except 
              select 
              Year, Region_A, Country_A, Region_B, Country_B
              from par_set
             ")
#This test should return "TRUE"
nrow(comp_set)-nrow(par_set)==nrow(dif)

#This code creates a full data set by combining the partial data set with the dif data set 
#created above. The added rows from dif are given "0" as their link intensity 
full_set <- sqldf("
                  select * from par_set 
                    union all 
                  select *, 0, 0 from dif
                  ")

#This test makes sure that there are no duplicate rows, grouped by year, country and region
#The test should return "TRUE"
test_unique <- sqldf("
                     select count(*) as rownumb 
                     from full_set 
                     group by Year, Region_A, Country_A, Region_B, Country_B
                     order by rownumb
                     ")
max(test_unique$rownumb)==1

#This code removes the temporary data frames
rm("comp_set", "dif", "test_unique")


# Partial set link intensity plot -----------------------------------------

#This code block will plot the link intesity per country (internal/external) per year using the partial 
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
plot(pc2c2,type = "o",col = "red", xlab = "Year", ylab = "Intensity", 
       main = "Network Intensity (partial)") 
  lines(pc2c1, type = "o", col = "blue") 
  lines(pc1c1, type = "o", col = "darkgreen")
  lines(pc1c2, type = "o", col = "orange")
  legend(1991, 45000, legend=c("C2:Internal", "C2:External", "C1:Internal", "C1:External"),
         col=c("red", "blue", "darkgreen", "orange"), lty=1:2, cex=0.5)

  
# Full set link intensity plot --------------------------------------------

#This code block will plot the link intesity per country (internal/external) per year using the full set
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
  plot(fc1c1,type = "o",col = "darkgreen", xlab = "Year", ylab = "Intensity", 
       main = "Network Intensity (full)") 
  lines(fc1c2, type = "o", col = "orange") 
  lines(fc2c2, type = "o", col = "red")
  lines(fc2c1, type = "o", col = "blue")
  legend(1991, 8500, legend=c("C2:Internal", "C2:External", "C1:Internal", "C1:External"),
         col=c("red", "blue", "darkgreen", "orange"), lty=1:2, cex=0.6)


# Compute and plot eigenvector centrality ---------------------------------

#This sections computes eigenvector centrality 
#Currently it is only an example of the internal link intensity for country 1 in 1991 in the partial data set
#To Do: expand for all years and combination of countries 

  #Country 1: Internal links 
  for (i in 1991:2017) {
    s <- paste("select * from par_set where Year = ", i, " and Country_A = 1 and Country_B = 1")
    x <- sqldf(s)
    g <- graph(c(x$Region_A, x$Region_B))
    eigenv <- paste("ev_11_", i, sep="")
    assign(eigenv, eigen_centrality(g, weights = x$Intensity_norm))
  }
  #Country 1: External links  
  for (i in 1991:2017) {
    s <- paste("select * from par_set where Year = ", i, " and Country_A = 1 and Country_B = 2")
    x <- sqldf(s)
    g <- graph(c(x$Region_A, x$Region_B))
    eigenv <- paste("ev_12_", i, sep="")
    assign(eigenv, eigen_centrality(g, weights = x$Intensity_norm))
  }
  #Country 2: Internal links  
  for (i in 1991:2017) {
    s <- paste("select * from par_set where Year = ", i, " and Country_A = 2 and Country_B = 2")
    x <- sqldf(s)
    g <- graph(c(x$Region_A, x$Region_B))
    eigenv <- paste("ev_22_", i, sep="")
    assign(eigenv, eigen_centrality(g, weights = x$Intensity_norm))
  }
  #Country 2: External links  
  for (i in 1991:2017) {
    s <- paste("select * from par_set where Year = ", i, " and Country_A = 2 and Country_B = 1")
    x <- sqldf(s)
    g <- graph(c(x$Region_A, x$Region_B))
    eigenv <- paste("ev_21_", i, sep="")
    assign(eigenv, eigen_centrality(g, weights = x$Intensity_norm))
  }
  remove(x,g,i,s,eigenv)

#compile eigenvalues into data frames  
year <- c(1991:2017)
#country 1 internal
eigenvalue <- c(ev_11_1991$value, 
                ev_11_1992$value, 
                ev_11_1993$value, 
                ev_11_1994$value, 
                ev_11_1995$value,
                ev_11_1996$value,
                ev_11_1997$value,
                ev_11_1998$value,
                ev_11_1999$value,
                ev_11_2000$value,
                ev_11_2001$value,
                ev_11_2002$value,
                ev_11_2003$value,
                ev_11_2004$value,
                ev_11_2005$value,
                ev_11_2006$value,
                ev_11_2007$value,
                ev_11_2008$value,
                ev_11_2009$value,
                ev_11_2010$value,
                ev_11_2011$value,
                ev_11_2012$value,
                ev_11_2013$value,
                ev_11_2014$value,
                ev_11_2015$value,
                ev_11_2016$value,
                ev_11_2017$value) 
eigenvalue11 <- data.frame(year, eigenvalue)
remove(
  ev_11_1991,
  ev_11_1992,
  ev_11_1993,
  ev_11_1994,
  ev_11_1995,
  ev_11_1996,
  ev_11_1997,
  ev_11_1998,
  ev_11_1999,
  ev_11_2000,
  ev_11_2001,
  ev_11_2002,
  ev_11_2003,
  ev_11_2004,
  ev_11_2005,
  ev_11_2006,
  ev_11_2007,
  ev_11_2008,
  ev_11_2009,
  ev_11_2010,
  ev_11_2011,
  ev_11_2012,
  ev_11_2013,
  ev_11_2014,
  ev_11_2015,
  ev_11_2016,
  ev_11_2017)
#country 1: external
eigenvalue <- c(ev_12_1991$value, 
                ev_12_1992$value, 
                ev_12_1993$value, 
                ev_12_1994$value, 
                ev_12_1995$value,
                ev_12_1996$value,
                ev_12_1997$value,
                ev_12_1998$value,
                ev_12_1999$value,
                ev_12_2000$value,
                ev_12_2001$value,
                ev_12_2002$value,
                ev_12_2003$value,
                ev_12_2004$value,
                ev_12_2005$value,
                ev_12_2006$value,
                ev_12_2007$value,
                ev_12_2008$value,
                ev_12_2009$value,
                ev_12_2010$value,
                ev_12_2011$value,
                ev_12_2012$value,
                ev_12_2013$value,
                ev_12_2014$value,
                ev_12_2015$value,
                ev_12_2016$value,
                ev_12_2017$value) 
eigenvalue12 <- data.frame(year, eigenvalue)
remove(
  ev_12_1991,
  ev_12_1992,
  ev_12_1993,
  ev_12_1994,
  ev_12_1995,
  ev_12_1996,
  ev_12_1997,
  ev_12_1998,
  ev_12_1999,
  ev_12_2000,
  ev_12_2001,
  ev_12_2002,
  ev_12_2003,
  ev_12_2004,
  ev_12_2005,
  ev_12_2006,
  ev_12_2007,
  ev_12_2008,
  ev_12_2009,
  ev_12_2010,
  ev_12_2011,
  ev_12_2012,
  ev_12_2013,
  ev_12_2014,
  ev_12_2015,
  ev_12_2016,
  ev_12_2017)
#country 2: internal
eigenvalue <- c(ev_22_1991$value, 
                ev_22_1992$value, 
                ev_22_1993$value, 
                ev_22_1994$value, 
                ev_22_1995$value,
                ev_22_1996$value,
                ev_22_1997$value,
                ev_22_1998$value,
                ev_22_1999$value,
                ev_22_2000$value,
                ev_22_2001$value,
                ev_22_2002$value,
                ev_22_2003$value,
                ev_22_2004$value,
                ev_22_2005$value,
                ev_22_2006$value,
                ev_22_2007$value,
                ev_22_2008$value,
                ev_22_2009$value,
                ev_22_2010$value,
                ev_22_2011$value,
                ev_22_2012$value,
                ev_22_2013$value,
                ev_22_2014$value,
                ev_22_2015$value,
                ev_22_2016$value,
                ev_22_2017$value) 
eigenvalue22 <- data.frame(year, eigenvalue)
remove(
  ev_22_1991,
  ev_22_1992,
  ev_22_1993,
  ev_22_1994,
  ev_22_1995,
  ev_22_1996,
  ev_22_1997,
  ev_22_1998,
  ev_22_1999,
  ev_22_2000,
  ev_22_2001,
  ev_22_2002,
  ev_22_2003,
  ev_22_2004,
  ev_22_2005,
  ev_22_2006,
  ev_22_2007,
  ev_22_2008,
  ev_22_2009,
  ev_22_2010,
  ev_22_2011,
  ev_22_2012,
  ev_22_2013,
  ev_22_2014,
  ev_22_2015,
  ev_22_2016,
  ev_22_2017)
#country 2: external 
eigenvalue <- c(ev_21_1991$value, 
                ev_21_1992$value, 
                ev_21_1993$value, 
                ev_21_1994$value, 
                ev_21_1995$value,
                ev_21_1996$value,
                ev_21_1997$value,
                ev_21_1998$value,
                ev_21_1999$value,
                ev_21_2000$value,
                ev_21_2001$value,
                ev_21_2002$value,
                ev_21_2003$value,
                ev_21_2004$value,
                ev_21_2005$value,
                ev_21_2006$value,
                ev_21_2007$value,
                ev_21_2008$value,
                ev_21_2009$value,
                ev_21_2010$value,
                ev_21_2011$value,
                ev_21_2012$value,
                ev_21_2013$value,
                ev_21_2014$value,
                ev_21_2015$value,
                ev_21_2016$value,
                ev_21_2017$value) 
eigenvalue21 <- data.frame(year, eigenvalue)
remove(
  ev_21_1991,
  ev_21_1992,
  ev_21_1993,
  ev_21_1994,
  ev_21_1995,
  ev_21_1996,
  ev_21_1997,
  ev_21_1998,
  ev_21_1999,
  ev_21_2000,
  ev_21_2001,
  ev_21_2002,
  ev_21_2003,
  ev_21_2004,
  ev_21_2005,
  ev_21_2006,
  ev_21_2007,
  ev_21_2008,
  ev_21_2009,
  ev_21_2010,
  ev_21_2011,
  ev_21_2012,
  ev_21_2013,
  ev_21_2014,
  ev_21_2015,
  ev_21_2016,
  ev_21_2017)
remove(year, eigenvalue)
        