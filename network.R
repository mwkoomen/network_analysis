
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
#install.packages("CINNA")

#Then load the packages into the current working session:
library(igraph)
library(dplyr)
library(sqldf)
library(readstata13)
library(CINNA)

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
         col=c("red", "blue", "darkgreen", "orange"), lty=1:2, cex=0.8)

  
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
         col=c("red", "blue", "darkgreen", "orange"), lty=1:2, cex=0.8)


# Construct nework graphs / adjacency matrix per year, country, and ex/int links------------------------------------------------
  nodes <- data.frame(ID=1:200)
  
  #Country 1: Internal links 
  for (i in 1991:2017) {
    s <- paste("select * from par_set where Year = ", i, " and Country_A = 1 and Country_B = 1")
    x <- sqldf(s)
    n <- paste("net_11_", i, sep="")
    links <- data.frame(RegionA=x$Region_A, RegionB=x$Region_B)
    assign(n, graph_from_data_frame(links, nodes, directed = TRUE))
  }
  #Country 1: External links  
  for (i in 1991:2017) {
    s <- paste("select * from par_set where Year = ", i, " and Country_A = 1 and Country_B = 2")
    x <- sqldf(s)
    n <- paste("net_12_", i, sep="")
    links <- data.frame(RegionA=x$Region_A, RegionB=x$Region_B)
    assign(n, graph_from_data_frame(links, nodes, directed = TRUE))
  }
  #Country 2: Internal links  
  for (i in 1991:2017) {
    s <- paste("select * from par_set where Year = ", i, " and Country_A = 2 and Country_B = 2")
    n <- paste("net_22_", i, sep="")
    links <- data.frame(RegionA=x$Region_A, RegionB=x$Region_B)
    assign(n, graph_from_data_frame(links, nodes, directed = TRUE))
  }
  #Country 2: External links  
  for (i in 1991:2017) {
    s <- paste("select * from par_set where Year = ", i, " and Country_A = 2 and Country_B = 1")
    x <- sqldf(s)
    n <- paste("net_21_", i, sep="")
    links <- data.frame(RegionA=x$Region_A, RegionB=x$Region_B)
    assign(n, graph_from_data_frame(links, nodes, directed = TRUE))
  }
  remove(x,i,s,links, nodes, net, n)
  
centralities <- proper_centralities(net_11_1991)
# Compute and plot degree centrality (works) --------------------------------------
  nodes <- data.frame(ID=1:200)

  #Country 1: Internal links 
  for (i in 1991:2017) {
    s <- paste("select * from par_set where Year = ", i, " and Country_A = 1 and Country_B = 1")
    x <- sqldf(s)
    links <- data.frame(RegionA=x$Region_A, RegionB=x$Region_B)
    net <- graph_from_data_frame(links, nodes, directed = TRUE)
    n <- paste("deg.cent_11_", i, sep="")
    assign(n, centr_degree(net, mode = "all", loops=TRUE, normalized = TRUE))
  }
  #Country 1: External links  
  for (i in 1991:2017) {
    s <- paste("select * from par_set where Year = ", i, " and Country_A = 1 and Country_B = 2")
    x <- sqldf(s)
    links <- data.frame(RegionA=x$Region_A, RegionB=x$Region_B)
    net <- graph_from_data_frame(links, nodes, directed = TRUE)
    n <- paste("deg.cent_12_", i, sep="")
    assign(n, centr_degree(net, mode = "all", loops=TRUE, normalized = TRUE))
  }
  #Country 2: Internal links  
  for (i in 1991:2017) {
    s <- paste("select * from par_set where Year = ", i, " and Country_A = 2 and Country_B = 2")
    x <- sqldf(s)
    links <- data.frame(RegionA=x$Region_A, RegionB=x$Region_B)
    net <- graph_from_data_frame(links, nodes, directed = TRUE)
    n <- paste("deg.cent_22_", i, sep="")
    assign(n, centr_degree(net, mode = "all", loops=TRUE, normalized = TRUE))
  }
  #Country 2: External links  
  for (i in 1991:2017) {
    s <- paste("select * from par_set where Year = ", i, " and Country_A = 2 and Country_B = 1")
    x <- sqldf(s)
    links <- data.frame(RegionA=x$Region_A, RegionB=x$Region_B)
    net <- graph_from_data_frame(links, nodes, directed = TRUE)
    n <- paste("deg.cent_21_", i, sep="")
    assign(n, centr_degree(net, mode = "all", loops=TRUE, normalized = TRUE))
  }
  remove(x,i,s,links, nodes, net, n)

  # Compile Degree Centrality (works) -----------------------------------------------
  #compile degree centrality  
  year <- c(1991:2017)
  #country 1: internal
  cn <- c(deg.cent_11_1991$centralization,
          deg.cent_11_1992$centralization,
          deg.cent_11_1993$centralization,
          deg.cent_11_1994$centralization,
          deg.cent_11_1995$centralization,
          deg.cent_11_1996$centralization,
          deg.cent_11_1997$centralization,
          deg.cent_11_1998$centralization,
          deg.cent_11_1999$centralization,
          deg.cent_11_2000$centralization,
          deg.cent_11_2001$centralization,
          deg.cent_11_2002$centralization,
          deg.cent_11_2003$centralization,
          deg.cent_11_2004$centralization,
          deg.cent_11_2005$centralization,
          deg.cent_11_2006$centralization,
          deg.cent_11_2007$centralization,
          deg.cent_11_2008$centralization,
          deg.cent_11_2009$centralization,
          deg.cent_11_2010$centralization,
          deg.cent_11_2011$centralization,
          deg.cent_11_2012$centralization,
          deg.cent_11_2013$centralization,
          deg.cent_11_2014$centralization,
          deg.cent_11_2015$centralization,
          deg.cent_11_2016$centralization,
          deg.cent_11_2017$centralization
  )
  degree_centrality11 <- data.frame(year, cn)
  #country 1: external
  cn <- c(deg.cent_12_1991$centralization,
          deg.cent_12_1992$centralization,
          deg.cent_12_1993$centralization,
          deg.cent_12_1994$centralization,
          deg.cent_12_1995$centralization,
          deg.cent_12_1996$centralization,
          deg.cent_12_1997$centralization,
          deg.cent_12_1998$centralization,
          deg.cent_12_1999$centralization,
          deg.cent_12_2000$centralization,
          deg.cent_12_2001$centralization,
          deg.cent_12_2002$centralization,
          deg.cent_12_2003$centralization,
          deg.cent_12_2004$centralization,
          deg.cent_12_2005$centralization,
          deg.cent_12_2006$centralization,
          deg.cent_12_2007$centralization,
          deg.cent_12_2008$centralization,
          deg.cent_12_2009$centralization,
          deg.cent_12_2010$centralization,
          deg.cent_12_2011$centralization,
          deg.cent_12_2012$centralization,
          deg.cent_12_2013$centralization,
          deg.cent_12_2014$centralization,
          deg.cent_12_2015$centralization,
          deg.cent_12_2016$centralization,
          deg.cent_12_2017$centralization
  )
  degree_centrality12 <- data.frame(year, cn)
  #country 2: internal
  cn <- c(deg.cent_22_1991$centralization,
          deg.cent_22_1992$centralization,
          deg.cent_22_1993$centralization,
          deg.cent_22_1994$centralization,
          deg.cent_22_1995$centralization,
          deg.cent_22_1996$centralization,
          deg.cent_22_1997$centralization,
          deg.cent_22_1998$centralization,
          deg.cent_22_1999$centralization,
          deg.cent_22_2000$centralization,
          deg.cent_22_2001$centralization,
          deg.cent_22_2002$centralization,
          deg.cent_22_2003$centralization,
          deg.cent_22_2004$centralization,
          deg.cent_22_2005$centralization,
          deg.cent_22_2006$centralization,
          deg.cent_22_2007$centralization,
          deg.cent_22_2008$centralization,
          deg.cent_22_2009$centralization,
          deg.cent_22_2010$centralization,
          deg.cent_22_2011$centralization,
          deg.cent_22_2012$centralization,
          deg.cent_22_2013$centralization,
          deg.cent_22_2014$centralization,
          deg.cent_22_2015$centralization,
          deg.cent_22_2016$centralization,
          deg.cent_22_2017$centralization
  )
  degree_centrality22 <- data.frame(year, cn)
  #country 2: external 
  cn <- c(deg.cent_21_1991$centralization,
          deg.cent_21_1992$centralization,
          deg.cent_21_1993$centralization,
          deg.cent_21_1994$centralization,
          deg.cent_21_1995$centralization,
          deg.cent_21_1996$centralization,
          deg.cent_21_1997$centralization,
          deg.cent_21_1998$centralization,
          deg.cent_21_1999$centralization,
          deg.cent_21_2000$centralization,
          deg.cent_21_2001$centralization,
          deg.cent_21_2002$centralization,
          deg.cent_21_2003$centralization,
          deg.cent_21_2004$centralization,
          deg.cent_21_2005$centralization,
          deg.cent_21_2006$centralization,
          deg.cent_21_2007$centralization,
          deg.cent_21_2008$centralization,
          deg.cent_21_2009$centralization,
          deg.cent_21_2010$centralization,
          deg.cent_21_2011$centralization,
          deg.cent_21_2012$centralization,
          deg.cent_21_2013$centralization,
          deg.cent_21_2014$centralization,
          deg.cent_21_2015$centralization,
          deg.cent_21_2016$centralization,
          deg.cent_21_2017$centralization
  )
  degree_centrality21 <- data.frame(year, cn)
  remove(year, cn)
  
  #compile degree centrality / theorectical max into data frames  
  year <- c(1991:2017)
  #country 1: internal
  cn <- c(deg.cent_11_1991$centralization/deg.cent_11_1991$theoretical_max, 
          deg.cent_11_1992$centralization/deg.cent_11_1992$theoretical_max, 
          deg.cent_11_1993$centralization/deg.cent_11_1993$theoretical_max, 
          deg.cent_11_1994$centralization/deg.cent_11_1994$theoretical_max, 
          deg.cent_11_1995$centralization/deg.cent_11_1995$theoretical_max,
          deg.cent_11_1996$centralization/deg.cent_11_1996$theoretical_max,
          deg.cent_11_1997$centralization/deg.cent_11_1997$theoretical_max,
          deg.cent_11_1998$centralization/deg.cent_11_1998$theoretical_max,
          deg.cent_11_1999$centralization/deg.cent_11_1999$theoretical_max,
          deg.cent_11_2000$centralization/deg.cent_11_2000$theoretical_max,
          deg.cent_11_2001$centralization/deg.cent_11_2001$theoretical_max,
          deg.cent_11_2002$centralization/deg.cent_11_2002$theoretical_max,
          deg.cent_11_2003$centralization/deg.cent_11_2003$theoretical_max,
          deg.cent_11_2004$centralization/deg.cent_11_2004$theoretical_max,
          deg.cent_11_2005$centralization/deg.cent_11_2005$theoretical_max,
          deg.cent_11_2006$centralization/deg.cent_11_2006$theoretical_max,
          deg.cent_11_2007$centralization/deg.cent_11_2007$theoretical_max,
          deg.cent_11_2008$centralization/deg.cent_11_2008$theoretical_max,
          deg.cent_11_2009$centralization/deg.cent_11_2009$theoretical_max,
          deg.cent_11_2010$centralization/deg.cent_11_2010$theoretical_max,
          deg.cent_11_2011$centralization/deg.cent_11_2011$theoretical_max,
          deg.cent_11_2012$centralization/deg.cent_11_2012$theoretical_max,
          deg.cent_11_2013$centralization/deg.cent_11_2013$theoretical_max,
          deg.cent_11_2014$centralization/deg.cent_11_2014$theoretical_max,
          deg.cent_11_2015$centralization/deg.cent_11_2015$theoretical_max,
          deg.cent_11_2016$centralization/deg.cent_11_2016$theoretical_max,
          deg.cent_11_2017$centralization/deg.cent_11_2017$theoretical_max
  )
  degree_centrality11_max <- data.frame(year, cn)
  remove(
    deg.cent_11_1991,
    deg.cent_11_1992,
    deg.cent_11_1993,
    deg.cent_11_1994,
    deg.cent_11_1995,
    deg.cent_11_1996,
    deg.cent_11_1997,
    deg.cent_11_1998,
    deg.cent_11_1999,
    deg.cent_11_2000,
    deg.cent_11_2001,
    deg.cent_11_2002,
    deg.cent_11_2003,
    deg.cent_11_2004,
    deg.cent_11_2005,
    deg.cent_11_2006,
    deg.cent_11_2007,
    deg.cent_11_2008,
    deg.cent_11_2009,
    deg.cent_11_2010,
    deg.cent_11_2011,
    deg.cent_11_2012,
    deg.cent_11_2013,
    deg.cent_11_2014,
    deg.cent_11_2015,
    deg.cent_11_2016,
    deg.cent_11_2017)
  #country 1: external
  cn <- c(deg.cent_12_1991$centralization/deg.cent_12_1991$theoretical_max, 
          deg.cent_12_1992$centralization/deg.cent_12_1992$theoretical_max, 
          deg.cent_12_1993$centralization/deg.cent_12_1993$theoretical_max, 
          deg.cent_12_1994$centralization/deg.cent_12_1994$theoretical_max, 
          deg.cent_12_1995$centralization/deg.cent_12_1995$theoretical_max,
          deg.cent_12_1996$centralization/deg.cent_12_1996$theoretical_max,
          deg.cent_12_1997$centralization/deg.cent_12_1997$theoretical_max,
          deg.cent_12_1998$centralization/deg.cent_12_1998$theoretical_max,
          deg.cent_12_1999$centralization/deg.cent_12_1999$theoretical_max,
          deg.cent_12_2000$centralization/deg.cent_12_2000$theoretical_max,
          deg.cent_12_2001$centralization/deg.cent_12_2001$theoretical_max,
          deg.cent_12_2002$centralization/deg.cent_12_2002$theoretical_max,
          deg.cent_12_2003$centralization/deg.cent_12_2003$theoretical_max,
          deg.cent_12_2004$centralization/deg.cent_12_2004$theoretical_max,
          deg.cent_12_2005$centralization/deg.cent_12_2005$theoretical_max,
          deg.cent_12_2006$centralization/deg.cent_12_2006$theoretical_max,
          deg.cent_12_2007$centralization/deg.cent_12_2007$theoretical_max,
          deg.cent_12_2008$centralization/deg.cent_12_2008$theoretical_max,
          deg.cent_12_2009$centralization/deg.cent_12_2009$theoretical_max,
          deg.cent_12_2010$centralization/deg.cent_12_2010$theoretical_max,
          deg.cent_12_2011$centralization/deg.cent_12_2011$theoretical_max,
          deg.cent_12_2012$centralization/deg.cent_12_2012$theoretical_max,
          deg.cent_12_2013$centralization/deg.cent_12_2013$theoretical_max,
          deg.cent_12_2014$centralization/deg.cent_12_2014$theoretical_max,
          deg.cent_12_2015$centralization/deg.cent_12_2015$theoretical_max,
          deg.cent_12_2016$centralization/deg.cent_12_2016$theoretical_max,
          deg.cent_12_2017$centralization/deg.cent_12_2017$theoretical_max
  )
  degree_centrality12_max <- data.frame(year, cn)
  remove(
    deg.cent_12_1991,
    deg.cent_12_1992,
    deg.cent_12_1993,
    deg.cent_12_1994,
    deg.cent_12_1995,
    deg.cent_12_1996,
    deg.cent_12_1997,
    deg.cent_12_1998,
    deg.cent_12_1999,
    deg.cent_12_2000,
    deg.cent_12_2001,
    deg.cent_12_2002,
    deg.cent_12_2003,
    deg.cent_12_2004,
    deg.cent_12_2005,
    deg.cent_12_2006,
    deg.cent_12_2007,
    deg.cent_12_2008,
    deg.cent_12_2009,
    deg.cent_12_2010,
    deg.cent_12_2011,
    deg.cent_12_2012,
    deg.cent_12_2013,
    deg.cent_12_2014,
    deg.cent_12_2015,
    deg.cent_12_2016,
    deg.cent_12_2017)
  #country 2: internal
  cn <- c(deg.cent_22_1991$centralization/deg.cent_22_1991$theoretical_max, 
          deg.cent_22_1992$centralization/deg.cent_22_1992$theoretical_max, 
          deg.cent_22_1993$centralization/deg.cent_22_1993$theoretical_max, 
          deg.cent_22_1994$centralization/deg.cent_22_1994$theoretical_max, 
          deg.cent_22_1995$centralization/deg.cent_22_1995$theoretical_max,
          deg.cent_22_1996$centralization/deg.cent_22_1996$theoretical_max,
          deg.cent_22_1997$centralization/deg.cent_22_1997$theoretical_max,
          deg.cent_22_1998$centralization/deg.cent_22_1998$theoretical_max,
          deg.cent_22_1999$centralization/deg.cent_22_1999$theoretical_max,
          deg.cent_22_2000$centralization/deg.cent_22_2000$theoretical_max,
          deg.cent_22_2001$centralization/deg.cent_22_2001$theoretical_max,
          deg.cent_22_2002$centralization/deg.cent_22_2002$theoretical_max,
          deg.cent_22_2003$centralization/deg.cent_22_2003$theoretical_max,
          deg.cent_22_2004$centralization/deg.cent_22_2004$theoretical_max,
          deg.cent_22_2005$centralization/deg.cent_22_2005$theoretical_max,
          deg.cent_22_2006$centralization/deg.cent_22_2006$theoretical_max,
          deg.cent_22_2007$centralization/deg.cent_22_2007$theoretical_max,
          deg.cent_22_2008$centralization/deg.cent_22_2008$theoretical_max,
          deg.cent_22_2009$centralization/deg.cent_22_2009$theoretical_max,
          deg.cent_22_2010$centralization/deg.cent_22_2010$theoretical_max,
          deg.cent_22_2011$centralization/deg.cent_22_2011$theoretical_max,
          deg.cent_22_2012$centralization/deg.cent_22_2012$theoretical_max,
          deg.cent_22_2013$centralization/deg.cent_22_2013$theoretical_max,
          deg.cent_22_2014$centralization/deg.cent_22_2014$theoretical_max,
          deg.cent_22_2015$centralization/deg.cent_22_2015$theoretical_max,
          deg.cent_22_2016$centralization/deg.cent_22_2016$theoretical_max,
          deg.cent_22_2017$centralization/deg.cent_22_2017$theoretical_max
  )
  degree_centrality22_max <- data.frame(year, cn)
  remove(
    deg.cent_22_1991,
    deg.cent_22_1992,
    deg.cent_22_1993,
    deg.cent_22_1994,
    deg.cent_22_1995,
    deg.cent_22_1996,
    deg.cent_22_1997,
    deg.cent_22_1998,
    deg.cent_22_1999,
    deg.cent_22_2000,
    deg.cent_22_2001,
    deg.cent_22_2002,
    deg.cent_22_2003,
    deg.cent_22_2004,
    deg.cent_22_2005,
    deg.cent_22_2006,
    deg.cent_22_2007,
    deg.cent_22_2008,
    deg.cent_22_2009,
    deg.cent_22_2010,
    deg.cent_22_2011,
    deg.cent_22_2012,
    deg.cent_22_2013,
    deg.cent_22_2014,
    deg.cent_22_2015,
    deg.cent_22_2016,
    deg.cent_22_2017)
  #country 2: external 
  cn <- c(deg.cent_21_1991$centralization/deg.cent_21_1991$theoretical_max, 
          deg.cent_21_1992$centralization/deg.cent_21_1992$theoretical_max, 
          deg.cent_21_1993$centralization/deg.cent_21_1993$theoretical_max, 
          deg.cent_21_1994$centralization/deg.cent_21_1994$theoretical_max, 
          deg.cent_21_1995$centralization/deg.cent_21_1995$theoretical_max,
          deg.cent_21_1996$centralization/deg.cent_21_1996$theoretical_max,
          deg.cent_21_1997$centralization/deg.cent_21_1997$theoretical_max,
          deg.cent_21_1998$centralization/deg.cent_21_1998$theoretical_max,
          deg.cent_21_1999$centralization/deg.cent_21_1999$theoretical_max,
          deg.cent_21_2000$centralization/deg.cent_21_2000$theoretical_max,
          deg.cent_21_2001$centralization/deg.cent_21_2001$theoretical_max,
          deg.cent_21_2002$centralization/deg.cent_21_2002$theoretical_max,
          deg.cent_21_2003$centralization/deg.cent_21_2003$theoretical_max,
          deg.cent_21_2004$centralization/deg.cent_21_2004$theoretical_max,
          deg.cent_21_2005$centralization/deg.cent_21_2005$theoretical_max,
          deg.cent_21_2006$centralization/deg.cent_21_2006$theoretical_max,
          deg.cent_21_2007$centralization/deg.cent_21_2007$theoretical_max,
          deg.cent_21_2008$centralization/deg.cent_21_2008$theoretical_max,
          deg.cent_21_2009$centralization/deg.cent_21_2009$theoretical_max,
          deg.cent_21_2010$centralization/deg.cent_21_2010$theoretical_max,
          deg.cent_21_2011$centralization/deg.cent_21_2011$theoretical_max,
          deg.cent_21_2012$centralization/deg.cent_21_2012$theoretical_max,
          deg.cent_21_2013$centralization/deg.cent_21_2013$theoretical_max,
          deg.cent_21_2014$centralization/deg.cent_21_2014$theoretical_max,
          deg.cent_21_2015$centralization/deg.cent_21_2015$theoretical_max,
          deg.cent_21_2016$centralization/deg.cent_21_2016$theoretical_max,
          deg.cent_21_2017$centralization/deg.cent_21_2017$theoretical_max
  )
  degree_centrality21_max <- data.frame(year, cn)
  remove(
    deg.cent_21_1991,
    deg.cent_21_1992,
    deg.cent_21_1993,
    deg.cent_21_1994,
    deg.cent_21_1995,
    deg.cent_21_1996,
    deg.cent_21_1997,
    deg.cent_21_1998,
    deg.cent_21_1999,
    deg.cent_21_2000,
    deg.cent_21_2001,
    deg.cent_21_2002,
    deg.cent_21_2003,
    deg.cent_21_2004,
    deg.cent_21_2005,
    deg.cent_21_2006,
    deg.cent_21_2007,
    deg.cent_21_2008,
    deg.cent_21_2009,
    deg.cent_21_2010,
    deg.cent_21_2011,
    deg.cent_21_2012,
    deg.cent_21_2013,
    deg.cent_21_2014,
    deg.cent_21_2015,
    deg.cent_21_2016,
    deg.cent_21_2017)
  remove(year, cn)


  # Plot Degree Centrality [integrated] (work-in-progress) --------------------------------------------------------------------

  #plot  
  plot(degree_centrality22$year,degree_centrality22$cn, type = "o",col = "red", xlab = "Year", ylab = "Degree Centrality", 
       main = "Degree centrality", ylim = c(0.05,0.25)) 
  lines(degree_centrality21$year, degree_centrality21$cn, type = "o", col = "blue") 
  lines(degree_centrality11$year, degree_centrality11$cn, type = "o", col = "darkgreen")
  lines(degree_centrality12$year, degree_centrality12$cn, type = "o", col = "orange")
  legend(2010, 0.43, legend=c("C2:Internal", "C2:External", "C1:Internal", "C1:External"),
         col=c("red", "blue", "darkgreen", "orange"), lty=1:2, cex=0.8)  
  
  #plot max 
  plot(degree_centrality22_max$year,degree_centrality22_max$cn, type = "o",col = "red", xlab = "Year", 
       ylab = "Degree Centrality / Theorectical max", 
       main = "Degree centrality", ylim=c(0.000004,0.000012)) 
  lines(degree_centrality21_max$year, degree_centrality21_max$cn, type = "o", col = "blue") 
  lines(degree_centrality11_max$year, degree_centrality11_max$cn, type = "o", col = "darkgreen")
  lines(degree_centrality12_max$year, degree_centrality12_max$cn, type = "o", col = "orange")
  legend(2010, 0.000012, legend=c("C2:Internal", "C2:External", "C1:Internal", "C1:External"),
         col=c("red", "blue", "darkgreen", "orange"), lty=1:2, cex=0.8)  
  



  # Compile Degree Centrality [integrated] (work-in-progress) -----------------------------------------------
  #compile degree centrality  
  year <- c(1991:2017)
  #country 1: internal
  cn <- c(deg.cent_11_1991$centralization,
          deg.cent_11_1992$centralization,
          deg.cent_11_1993$centralization,
          deg.cent_11_1994$centralization,
          deg.cent_11_1995$centralization,
          deg.cent_11_1996$centralization,
          deg.cent_11_1997$centralization,
          deg.cent_11_1998$centralization,
          deg.cent_11_1999$centralization,
          deg.cent_11_2000$centralization,
          deg.cent_11_2001$centralization,
          deg.cent_11_2002$centralization,
          deg.cent_11_2003$centralization,
          deg.cent_11_2004$centralization,
          deg.cent_11_2005$centralization,
          deg.cent_11_2006$centralization,
          deg.cent_11_2007$centralization,
          deg.cent_11_2008$centralization,
          deg.cent_11_2009$centralization,
          deg.cent_11_2010$centralization,
          deg.cent_11_2011$centralization,
          deg.cent_11_2012$centralization,
          deg.cent_11_2013$centralization,
          deg.cent_11_2014$centralization,
          deg.cent_11_2015$centralization,
          deg.cent_11_2016$centralization,
          deg.cent_11_2017$centralization
  )
  degree_centrality11 <- data.frame(year, cn)
  #country 1: external
  cn <- c(deg.cent_12_1991$centralization,
          deg.cent_12_1992$centralization,
          deg.cent_12_1993$centralization,
          deg.cent_12_1994$centralization,
          deg.cent_12_1995$centralization,
          deg.cent_12_1996$centralization,
          deg.cent_12_1997$centralization,
          deg.cent_12_1998$centralization,
          deg.cent_12_1999$centralization,
          deg.cent_12_2000$centralization,
          deg.cent_12_2001$centralization,
          deg.cent_12_2002$centralization,
          deg.cent_12_2003$centralization,
          deg.cent_12_2004$centralization,
          deg.cent_12_2005$centralization,
          deg.cent_12_2006$centralization,
          deg.cent_12_2007$centralization,
          deg.cent_12_2008$centralization,
          deg.cent_12_2009$centralization,
          deg.cent_12_2010$centralization,
          deg.cent_12_2011$centralization,
          deg.cent_12_2012$centralization,
          deg.cent_12_2013$centralization,
          deg.cent_12_2014$centralization,
          deg.cent_12_2015$centralization,
          deg.cent_12_2016$centralization,
          deg.cent_12_2017$centralization
  )
  degree_centrality12 <- data.frame(year, cn)
  #country 2: internal
  cn <- c(deg.cent_22_1991$centralization,
          deg.cent_22_1992$centralization,
          deg.cent_22_1993$centralization,
          deg.cent_22_1994$centralization,
          deg.cent_22_1995$centralization,
          deg.cent_22_1996$centralization,
          deg.cent_22_1997$centralization,
          deg.cent_22_1998$centralization,
          deg.cent_22_1999$centralization,
          deg.cent_22_2000$centralization,
          deg.cent_22_2001$centralization,
          deg.cent_22_2002$centralization,
          deg.cent_22_2003$centralization,
          deg.cent_22_2004$centralization,
          deg.cent_22_2005$centralization,
          deg.cent_22_2006$centralization,
          deg.cent_22_2007$centralization,
          deg.cent_22_2008$centralization,
          deg.cent_22_2009$centralization,
          deg.cent_22_2010$centralization,
          deg.cent_22_2011$centralization,
          deg.cent_22_2012$centralization,
          deg.cent_22_2013$centralization,
          deg.cent_22_2014$centralization,
          deg.cent_22_2015$centralization,
          deg.cent_22_2016$centralization,
          deg.cent_22_2017$centralization
  )
  degree_centrality22 <- data.frame(year, cn)
  #country 2: external 
  cn <- c(deg.cent_21_1991$centralization,
          deg.cent_21_1992$centralization,
          deg.cent_21_1993$centralization,
          deg.cent_21_1994$centralization,
          deg.cent_21_1995$centralization,
          deg.cent_21_1996$centralization,
          deg.cent_21_1997$centralization,
          deg.cent_21_1998$centralization,
          deg.cent_21_1999$centralization,
          deg.cent_21_2000$centralization,
          deg.cent_21_2001$centralization,
          deg.cent_21_2002$centralization,
          deg.cent_21_2003$centralization,
          deg.cent_21_2004$centralization,
          deg.cent_21_2005$centralization,
          deg.cent_21_2006$centralization,
          deg.cent_21_2007$centralization,
          deg.cent_21_2008$centralization,
          deg.cent_21_2009$centralization,
          deg.cent_21_2010$centralization,
          deg.cent_21_2011$centralization,
          deg.cent_21_2012$centralization,
          deg.cent_21_2013$centralization,
          deg.cent_21_2014$centralization,
          deg.cent_21_2015$centralization,
          deg.cent_21_2016$centralization,
          deg.cent_21_2017$centralization
  )
  degree_centrality21 <- data.frame(year, cn)
  remove(year, cn)
  
  #compile degree centrality / theorectical max into data frames  
  year <- c(1991:2017)
  #country 1: internal
  cn <- c(deg.cent_11_1991$centralization/deg.cent_11_1991$theoretical_max, 
          deg.cent_11_1992$centralization/deg.cent_11_1992$theoretical_max, 
          deg.cent_11_1993$centralization/deg.cent_11_1993$theoretical_max, 
          deg.cent_11_1994$centralization/deg.cent_11_1994$theoretical_max, 
          deg.cent_11_1995$centralization/deg.cent_11_1995$theoretical_max,
          deg.cent_11_1996$centralization/deg.cent_11_1996$theoretical_max,
          deg.cent_11_1997$centralization/deg.cent_11_1997$theoretical_max,
          deg.cent_11_1998$centralization/deg.cent_11_1998$theoretical_max,
          deg.cent_11_1999$centralization/deg.cent_11_1999$theoretical_max,
          deg.cent_11_2000$centralization/deg.cent_11_2000$theoretical_max,
          deg.cent_11_2001$centralization/deg.cent_11_2001$theoretical_max,
          deg.cent_11_2002$centralization/deg.cent_11_2002$theoretical_max,
          deg.cent_11_2003$centralization/deg.cent_11_2003$theoretical_max,
          deg.cent_11_2004$centralization/deg.cent_11_2004$theoretical_max,
          deg.cent_11_2005$centralization/deg.cent_11_2005$theoretical_max,
          deg.cent_11_2006$centralization/deg.cent_11_2006$theoretical_max,
          deg.cent_11_2007$centralization/deg.cent_11_2007$theoretical_max,
          deg.cent_11_2008$centralization/deg.cent_11_2008$theoretical_max,
          deg.cent_11_2009$centralization/deg.cent_11_2009$theoretical_max,
          deg.cent_11_2010$centralization/deg.cent_11_2010$theoretical_max,
          deg.cent_11_2011$centralization/deg.cent_11_2011$theoretical_max,
          deg.cent_11_2012$centralization/deg.cent_11_2012$theoretical_max,
          deg.cent_11_2013$centralization/deg.cent_11_2013$theoretical_max,
          deg.cent_11_2014$centralization/deg.cent_11_2014$theoretical_max,
          deg.cent_11_2015$centralization/deg.cent_11_2015$theoretical_max,
          deg.cent_11_2016$centralization/deg.cent_11_2016$theoretical_max,
          deg.cent_11_2017$centralization/deg.cent_11_2017$theoretical_max
  )
  degree_centrality11_max <- data.frame(year, cn)
  remove(
    deg.cent_11_1991,
    deg.cent_11_1992,
    deg.cent_11_1993,
    deg.cent_11_1994,
    deg.cent_11_1995,
    deg.cent_11_1996,
    deg.cent_11_1997,
    deg.cent_11_1998,
    deg.cent_11_1999,
    deg.cent_11_2000,
    deg.cent_11_2001,
    deg.cent_11_2002,
    deg.cent_11_2003,
    deg.cent_11_2004,
    deg.cent_11_2005,
    deg.cent_11_2006,
    deg.cent_11_2007,
    deg.cent_11_2008,
    deg.cent_11_2009,
    deg.cent_11_2010,
    deg.cent_11_2011,
    deg.cent_11_2012,
    deg.cent_11_2013,
    deg.cent_11_2014,
    deg.cent_11_2015,
    deg.cent_11_2016,
    deg.cent_11_2017)
  #country 1: external
  cn <- c(deg.cent_12_1991$centralization/deg.cent_12_1991$theoretical_max, 
          deg.cent_12_1992$centralization/deg.cent_12_1992$theoretical_max, 
          deg.cent_12_1993$centralization/deg.cent_12_1993$theoretical_max, 
          deg.cent_12_1994$centralization/deg.cent_12_1994$theoretical_max, 
          deg.cent_12_1995$centralization/deg.cent_12_1995$theoretical_max,
          deg.cent_12_1996$centralization/deg.cent_12_1996$theoretical_max,
          deg.cent_12_1997$centralization/deg.cent_12_1997$theoretical_max,
          deg.cent_12_1998$centralization/deg.cent_12_1998$theoretical_max,
          deg.cent_12_1999$centralization/deg.cent_12_1999$theoretical_max,
          deg.cent_12_2000$centralization/deg.cent_12_2000$theoretical_max,
          deg.cent_12_2001$centralization/deg.cent_12_2001$theoretical_max,
          deg.cent_12_2002$centralization/deg.cent_12_2002$theoretical_max,
          deg.cent_12_2003$centralization/deg.cent_12_2003$theoretical_max,
          deg.cent_12_2004$centralization/deg.cent_12_2004$theoretical_max,
          deg.cent_12_2005$centralization/deg.cent_12_2005$theoretical_max,
          deg.cent_12_2006$centralization/deg.cent_12_2006$theoretical_max,
          deg.cent_12_2007$centralization/deg.cent_12_2007$theoretical_max,
          deg.cent_12_2008$centralization/deg.cent_12_2008$theoretical_max,
          deg.cent_12_2009$centralization/deg.cent_12_2009$theoretical_max,
          deg.cent_12_2010$centralization/deg.cent_12_2010$theoretical_max,
          deg.cent_12_2011$centralization/deg.cent_12_2011$theoretical_max,
          deg.cent_12_2012$centralization/deg.cent_12_2012$theoretical_max,
          deg.cent_12_2013$centralization/deg.cent_12_2013$theoretical_max,
          deg.cent_12_2014$centralization/deg.cent_12_2014$theoretical_max,
          deg.cent_12_2015$centralization/deg.cent_12_2015$theoretical_max,
          deg.cent_12_2016$centralization/deg.cent_12_2016$theoretical_max,
          deg.cent_12_2017$centralization/deg.cent_12_2017$theoretical_max
  )
  degree_centrality12_max <- data.frame(year, cn)
  remove(
    deg.cent_12_1991,
    deg.cent_12_1992,
    deg.cent_12_1993,
    deg.cent_12_1994,
    deg.cent_12_1995,
    deg.cent_12_1996,
    deg.cent_12_1997,
    deg.cent_12_1998,
    deg.cent_12_1999,
    deg.cent_12_2000,
    deg.cent_12_2001,
    deg.cent_12_2002,
    deg.cent_12_2003,
    deg.cent_12_2004,
    deg.cent_12_2005,
    deg.cent_12_2006,
    deg.cent_12_2007,
    deg.cent_12_2008,
    deg.cent_12_2009,
    deg.cent_12_2010,
    deg.cent_12_2011,
    deg.cent_12_2012,
    deg.cent_12_2013,
    deg.cent_12_2014,
    deg.cent_12_2015,
    deg.cent_12_2016,
    deg.cent_12_2017)
  #country 2: internal
  cn <- c(deg.cent_22_1991$centralization/deg.cent_22_1991$theoretical_max, 
          deg.cent_22_1992$centralization/deg.cent_22_1992$theoretical_max, 
          deg.cent_22_1993$centralization/deg.cent_22_1993$theoretical_max, 
          deg.cent_22_1994$centralization/deg.cent_22_1994$theoretical_max, 
          deg.cent_22_1995$centralization/deg.cent_22_1995$theoretical_max,
          deg.cent_22_1996$centralization/deg.cent_22_1996$theoretical_max,
          deg.cent_22_1997$centralization/deg.cent_22_1997$theoretical_max,
          deg.cent_22_1998$centralization/deg.cent_22_1998$theoretical_max,
          deg.cent_22_1999$centralization/deg.cent_22_1999$theoretical_max,
          deg.cent_22_2000$centralization/deg.cent_22_2000$theoretical_max,
          deg.cent_22_2001$centralization/deg.cent_22_2001$theoretical_max,
          deg.cent_22_2002$centralization/deg.cent_22_2002$theoretical_max,
          deg.cent_22_2003$centralization/deg.cent_22_2003$theoretical_max,
          deg.cent_22_2004$centralization/deg.cent_22_2004$theoretical_max,
          deg.cent_22_2005$centralization/deg.cent_22_2005$theoretical_max,
          deg.cent_22_2006$centralization/deg.cent_22_2006$theoretical_max,
          deg.cent_22_2007$centralization/deg.cent_22_2007$theoretical_max,
          deg.cent_22_2008$centralization/deg.cent_22_2008$theoretical_max,
          deg.cent_22_2009$centralization/deg.cent_22_2009$theoretical_max,
          deg.cent_22_2010$centralization/deg.cent_22_2010$theoretical_max,
          deg.cent_22_2011$centralization/deg.cent_22_2011$theoretical_max,
          deg.cent_22_2012$centralization/deg.cent_22_2012$theoretical_max,
          deg.cent_22_2013$centralization/deg.cent_22_2013$theoretical_max,
          deg.cent_22_2014$centralization/deg.cent_22_2014$theoretical_max,
          deg.cent_22_2015$centralization/deg.cent_22_2015$theoretical_max,
          deg.cent_22_2016$centralization/deg.cent_22_2016$theoretical_max,
          deg.cent_22_2017$centralization/deg.cent_22_2017$theoretical_max
  )
  degree_centrality22_max <- data.frame(year, cn)
  remove(
    deg.cent_22_1991,
    deg.cent_22_1992,
    deg.cent_22_1993,
    deg.cent_22_1994,
    deg.cent_22_1995,
    deg.cent_22_1996,
    deg.cent_22_1997,
    deg.cent_22_1998,
    deg.cent_22_1999,
    deg.cent_22_2000,
    deg.cent_22_2001,
    deg.cent_22_2002,
    deg.cent_22_2003,
    deg.cent_22_2004,
    deg.cent_22_2005,
    deg.cent_22_2006,
    deg.cent_22_2007,
    deg.cent_22_2008,
    deg.cent_22_2009,
    deg.cent_22_2010,
    deg.cent_22_2011,
    deg.cent_22_2012,
    deg.cent_22_2013,
    deg.cent_22_2014,
    deg.cent_22_2015,
    deg.cent_22_2016,
    deg.cent_22_2017)
  #country 2: external 
  cn <- c(deg.cent_21_1991$centralization/deg.cent_21_1991$theoretical_max, 
          deg.cent_21_1992$centralization/deg.cent_21_1992$theoretical_max, 
          deg.cent_21_1993$centralization/deg.cent_21_1993$theoretical_max, 
          deg.cent_21_1994$centralization/deg.cent_21_1994$theoretical_max, 
          deg.cent_21_1995$centralization/deg.cent_21_1995$theoretical_max,
          deg.cent_21_1996$centralization/deg.cent_21_1996$theoretical_max,
          deg.cent_21_1997$centralization/deg.cent_21_1997$theoretical_max,
          deg.cent_21_1998$centralization/deg.cent_21_1998$theoretical_max,
          deg.cent_21_1999$centralization/deg.cent_21_1999$theoretical_max,
          deg.cent_21_2000$centralization/deg.cent_21_2000$theoretical_max,
          deg.cent_21_2001$centralization/deg.cent_21_2001$theoretical_max,
          deg.cent_21_2002$centralization/deg.cent_21_2002$theoretical_max,
          deg.cent_21_2003$centralization/deg.cent_21_2003$theoretical_max,
          deg.cent_21_2004$centralization/deg.cent_21_2004$theoretical_max,
          deg.cent_21_2005$centralization/deg.cent_21_2005$theoretical_max,
          deg.cent_21_2006$centralization/deg.cent_21_2006$theoretical_max,
          deg.cent_21_2007$centralization/deg.cent_21_2007$theoretical_max,
          deg.cent_21_2008$centralization/deg.cent_21_2008$theoretical_max,
          deg.cent_21_2009$centralization/deg.cent_21_2009$theoretical_max,
          deg.cent_21_2010$centralization/deg.cent_21_2010$theoretical_max,
          deg.cent_21_2011$centralization/deg.cent_21_2011$theoretical_max,
          deg.cent_21_2012$centralization/deg.cent_21_2012$theoretical_max,
          deg.cent_21_2013$centralization/deg.cent_21_2013$theoretical_max,
          deg.cent_21_2014$centralization/deg.cent_21_2014$theoretical_max,
          deg.cent_21_2015$centralization/deg.cent_21_2015$theoretical_max,
          deg.cent_21_2016$centralization/deg.cent_21_2016$theoretical_max,
          deg.cent_21_2017$centralization/deg.cent_21_2017$theoretical_max
  )
  degree_centrality21_max <- data.frame(year, cn)
  remove(
    deg.cent_21_1991,
    deg.cent_21_1992,
    deg.cent_21_1993,
    deg.cent_21_1994,
    deg.cent_21_1995,
    deg.cent_21_1996,
    deg.cent_21_1997,
    deg.cent_21_1998,
    deg.cent_21_1999,
    deg.cent_21_2000,
    deg.cent_21_2001,
    deg.cent_21_2002,
    deg.cent_21_2003,
    deg.cent_21_2004,
    deg.cent_21_2005,
    deg.cent_21_2006,
    deg.cent_21_2007,
    deg.cent_21_2008,
    deg.cent_21_2009,
    deg.cent_21_2010,
    deg.cent_21_2011,
    deg.cent_21_2012,
    deg.cent_21_2013,
    deg.cent_21_2014,
    deg.cent_21_2015,
    deg.cent_21_2016,
    deg.cent_21_2017)
  remove(year, cn)
  
  
  # Plot Degree Centrality [integrated] (work-in-progress) --------------------------------------------------------------------
  
  #plot  
  plot(degree_centrality22$year,degree_centrality22$cn, type = "o",col = "red", xlab = "Year", ylab = "Degree Centrality", 
       main = "Degree centrality", ylim = c(0.05,0.25)) 
  lines(degree_centrality21$year, degree_centrality21$cn, type = "o", col = "blue") 
  lines(degree_centrality11$year, degree_centrality11$cn, type = "o", col = "darkgreen")
  lines(degree_centrality12$year, degree_centrality12$cn, type = "o", col = "orange")
  legend(2010, 0.43, legend=c("C2:Internal", "C2:External", "C1:Internal", "C1:External"),
         col=c("red", "blue", "darkgreen", "orange"), lty=1:2, cex=0.8)  
  
  #plot max 
  plot(degree_centrality22_max$year,degree_centrality22_max$cn, type = "o",col = "red", xlab = "Year", 
       ylab = "Degree Centrality / Theorectical max", 
       main = "Degree centrality", ylim=c(0.000004,0.000012)) 
  lines(degree_centrality21_max$year, degree_centrality21_max$cn, type = "o", col = "blue") 
  lines(degree_centrality11_max$year, degree_centrality11_max$cn, type = "o", col = "darkgreen")
  lines(degree_centrality12_max$year, degree_centrality12_max$cn, type = "o", col = "orange")
  legend(2010, 0.000012, legend=c("C2:Internal", "C2:External", "C1:Internal", "C1:External"),
         col=c("red", "blue", "darkgreen", "orange"), lty=1:2, cex=0.8)  
  
  
  
  
  
# Compute Harmonic Centrality (works) ---------------------------------------------
#C1: internal
harm_cent11 <- data.frame()
  for (i in 1991:2017){
  n <- paste("net_11_", i, sep="")
  x <- get(n)
  h <- paste("harm.cent_11_", i, sep="")
  assign(h, harmonic_centrality(x, mode = "all", weights = NULL))
  d <- data.frame(year=i, harmonic_centrality=mean(get(h)))
  harm_cent11 <- rbind(harm_cent11, d)
}
  #C1: external
  harm_cent12 <- data.frame()
  for (i in 1991:2017){
    n <- paste("net_12_", i, sep="")
    x <- get(n)
    h <- paste("harm.cent_12_", i, sep="")
    assign(h, harmonic_centrality(x, mode = "all", weights = NULL))
    d <- data.frame(year=i, harmonic_centrality=mean(get(h)))
    harm_cent12 <- rbind(harm_cent12, d)
  }
    #C2: internal
    harm_cent22 <- data.frame()
    for (i in 1991:2017){
      n <- paste("net_22_", i, sep="")
      x <- get(n)
      h <- paste("harm.cent_22_", i, sep="")
      assign(h, harmonic_centrality(x, mode = "all", weights = NULL))
      d <- data.frame(year=i, harmonic_centrality=mean(get(h)))
      harm_cent22 <- rbind(harm_cent22, d)
    }
      #C2: external      
      harm_cent21 <- data.frame()
      for (i in 1991:2017){
        n <- paste("net_21_", i, sep="")
        x <- get(n)
        h <- paste("harm.cent_21_", i, sep="")
        assign(h, harmonic_centrality(x, mode = "all", weights = NULL))
        d <- data.frame(year=i, harmonic_centrality=mean(get(h)))
        harm_cent21 <- rbind(harm_cent21, d)
      }    
      plot(harm_cent22, type = "o",col = "red", xlab = "Year", ylab = "Harmonic Centrality (Mean)", 
           main = "Harmonic centrality (Mean)") 
      lines(harm_cent21, type = "o", col = "blue") 
      lines(harm_cent11, type = "o", col = "darkgreen")
      lines(harm_cent12, type = "o", col = "orange")
      legend(1991, 140, legend=c("C2:Internal", "C2:External", "C1:Internal", "C1:External"),
             col=c("red", "blue", "darkgreen", "orange"), lty=1:2, cex=0.8) 

      par(mfrow=c(4,7))
      plot(harm.cent_11_1991, xlab="", ylab = "", main = "1991")
      plot(harm.cent_11_1992, xlab="", ylab = "", main = "1992")
      plot(harm.cent_11_1993, xlab="", ylab = "", main = "1993")
      plot(harm.cent_11_1994, xlab="", ylab = "", main = "1994")
      plot(harm.cent_11_1995, xlab="", ylab = "", main = "1995")
      plot(harm.cent_11_1996, xlab="", ylab = "", main = "1996")
      plot(harm.cent_11_1997, xlab="", ylab = "", main = "1997")
      plot(harm.cent_11_1998, xlab="", ylab = "", main = "1998")
      plot(harm.cent_11_1999, xlab="", ylab = "", main = "1999")
      plot(harm.cent_11_2000, xlab="", ylab = "", main = "2000")
      plot(harm.cent_11_2001, xlab="", ylab = "", main = "2001")
      plot(harm.cent_11_2002, xlab="", ylab = "", main = "2002")
      plot(harm.cent_11_2003, xlab="", ylab = "", main = "2003")
      plot(harm.cent_11_2004, xlab="", ylab = "", main = "2004")
      plot(harm.cent_11_2005, xlab="", ylab = "", main = "2005")
      plot(harm.cent_11_2006, xlab="", ylab = "", main = "2006")
      plot(harm.cent_11_2007, xlab="", ylab = "", main = "2007")
      plot(harm.cent_11_2008, xlab="", ylab = "", main = "2008")
      plot(harm.cent_11_2009, xlab="", ylab = "", main = "2009")
      plot(harm.cent_11_2010, xlab="", ylab = "", main = "2010")
      plot(harm.cent_11_2011, xlab="", ylab = "", main = "2011")
      plot(harm.cent_11_2012, xlab="", ylab = "", main = "2012")
      plot(harm.cent_11_2013, xlab="", ylab = "", main = "2013")
      plot(harm.cent_11_2014, xlab="", ylab = "", main = "2014")
      plot(harm.cent_11_2015, xlab="", ylab = "", main = "2015")
      plot(harm.cent_11_2016, xlab="", ylab = "", main = "2016")
      plot(harm.cent_11_2017, xlab="", ylab = "", main = "2017")
     
              
# Compute closeness centrality (not recommended)--------------------------------------
  nodes <- data.frame(ID=1:200)
  
  #Country 1: Internal links 
  for (i in 1991:2017) {
    s <- paste("select * from par_set where Year = ", i, " and Country_A = 1 and Country_B = 1")
    x <- sqldf(s)
    links <- data.frame(RegionA=x$Region_A, RegionB=x$Region_B)
    net <- graph_from_data_frame(links, nodes, directed = FALSE)
    n <- paste("close.cent_11_", i, sep="")
    assign(n, closeness(net, mode = "all", weights = NULL, normalized = TRUE))
  }
  #Country 1: External links  
  for (i in 1991:2017) {
    s <- paste("select * from par_set where Year = ", i, " and Country_A = 1 and Country_B = 2")
    x <- sqldf(s)
    links <- data.frame(RegionA=x$Region_A, RegionB=x$Region_B)
    net <- graph_from_data_frame(links, nodes, directed = FALSE)
    n <- paste("close.cent_12_", i, sep="")
    assign(n, closeness(net, mode = "all", weights = NULL, normalized = TRUE))
  }
  #Country 2: Internal links  
  for (i in 1991:2017) {
    s <- paste("select * from par_set where Year = ", i, " and Country_A = 2 and Country_B = 2")
    x <- sqldf(s)
    links <- data.frame(RegionA=x$Region_A, RegionB=x$Region_B)
    net <- graph_from_data_frame(links, nodes, directed = FALSE)
    n <- paste("close.cent_22_", i, sep="")
    assign(n, closeness(net, mode = "all", weights = NULL, normalized = TRUE))
  }
  #Country 2: External links  
  for (i in 1991:2017) {
    s <- paste("select * from par_set where Year = ", i, " and Country_A = 2 and Country_B = 1")
    x <- sqldf(s)
    links <- data.frame(RegionA=x$Region_A, RegionB=x$Region_B)
    net <- graph_from_data_frame(links, nodes, directed = FALSE)
    n <- paste("close.cent_21_", i, sep="")
    assign(n, closeness(net, mode = "all", weights = NULL, normalized = TRUE))
  }
  remove(x,i,s,links, nodes, net, n)
  
  # Compile CLoseness Centrality (not recommended)-----------------------------------------------

year <- c(1991:2017)
  #country 1: internal
  cn <- c(deg.cent_11_1991$centralization,
          deg.cent_11_1992$centralization,
          deg.cent_11_1993$centralization,
          deg.cent_11_1994$centralization,
          deg.cent_11_1995$centralization,
          deg.cent_11_1996$centralization,
          deg.cent_11_1997$centralization,
          deg.cent_11_1998$centralization,
          deg.cent_11_1999$centralization,
          deg.cent_11_2000$centralization,
          deg.cent_11_2001$centralization,
          deg.cent_11_2002$centralization,
          deg.cent_11_2003$centralization,
          deg.cent_11_2004$centralization,
          deg.cent_11_2005$centralization,
          deg.cent_11_2006$centralization,
          deg.cent_11_2007$centralization,
          deg.cent_11_2008$centralization,
          deg.cent_11_2009$centralization,
          deg.cent_11_2010$centralization,
          deg.cent_11_2011$centralization,
          deg.cent_11_2012$centralization,
          deg.cent_11_2013$centralization,
          deg.cent_11_2014$centralization,
          deg.cent_11_2015$centralization,
          deg.cent_11_2016$centralization,
          deg.cent_11_2017$centralization
  )
  degree_centrality11 <- data.frame(year, cn)
  #country 1: external
  cn <- c(deg.cent_12_1991$centralization,
          deg.cent_12_1992$centralization,
          deg.cent_12_1993$centralization,
          deg.cent_12_1994$centralization,
          deg.cent_12_1995$centralization,
          deg.cent_12_1996$centralization,
          deg.cent_12_1997$centralization,
          deg.cent_12_1998$centralization,
          deg.cent_12_1999$centralization,
          deg.cent_12_2000$centralization,
          deg.cent_12_2001$centralization,
          deg.cent_12_2002$centralization,
          deg.cent_12_2003$centralization,
          deg.cent_12_2004$centralization,
          deg.cent_12_2005$centralization,
          deg.cent_12_2006$centralization,
          deg.cent_12_2007$centralization,
          deg.cent_12_2008$centralization,
          deg.cent_12_2009$centralization,
          deg.cent_12_2010$centralization,
          deg.cent_12_2011$centralization,
          deg.cent_12_2012$centralization,
          deg.cent_12_2013$centralization,
          deg.cent_12_2014$centralization,
          deg.cent_12_2015$centralization,
          deg.cent_12_2016$centralization,
          deg.cent_12_2017$centralization
  )
  degree_centrality12 <- data.frame(year, cn)
  #country 2: internal
  cn <- c(deg.cent_22_1991$centralization,
          deg.cent_22_1992$centralization,
          deg.cent_22_1993$centralization,
          deg.cent_22_1994$centralization,
          deg.cent_22_1995$centralization,
          deg.cent_22_1996$centralization,
          deg.cent_22_1997$centralization,
          deg.cent_22_1998$centralization,
          deg.cent_22_1999$centralization,
          deg.cent_22_2000$centralization,
          deg.cent_22_2001$centralization,
          deg.cent_22_2002$centralization,
          deg.cent_22_2003$centralization,
          deg.cent_22_2004$centralization,
          deg.cent_22_2005$centralization,
          deg.cent_22_2006$centralization,
          deg.cent_22_2007$centralization,
          deg.cent_22_2008$centralization,
          deg.cent_22_2009$centralization,
          deg.cent_22_2010$centralization,
          deg.cent_22_2011$centralization,
          deg.cent_22_2012$centralization,
          deg.cent_22_2013$centralization,
          deg.cent_22_2014$centralization,
          deg.cent_22_2015$centralization,
          deg.cent_22_2016$centralization,
          deg.cent_22_2017$centralization
  )
  degree_centrality22 <- data.frame(year, cn)
  #country 2: external 
  cn <- c(deg.cent_21_1991$centralization,
          deg.cent_21_1992$centralization,
          deg.cent_21_1993$centralization,
          deg.cent_21_1994$centralization,
          deg.cent_21_1995$centralization,
          deg.cent_21_1996$centralization,
          deg.cent_21_1997$centralization,
          deg.cent_21_1998$centralization,
          deg.cent_21_1999$centralization,
          deg.cent_21_2000$centralization,
          deg.cent_21_2001$centralization,
          deg.cent_21_2002$centralization,
          deg.cent_21_2003$centralization,
          deg.cent_21_2004$centralization,
          deg.cent_21_2005$centralization,
          deg.cent_21_2006$centralization,
          deg.cent_21_2007$centralization,
          deg.cent_21_2008$centralization,
          deg.cent_21_2009$centralization,
          deg.cent_21_2010$centralization,
          deg.cent_21_2011$centralization,
          deg.cent_21_2012$centralization,
          deg.cent_21_2013$centralization,
          deg.cent_21_2014$centralization,
          deg.cent_21_2015$centralization,
          deg.cent_21_2016$centralization,
          deg.cent_21_2017$centralization
  )
  degree_centrality21 <- data.frame(year, cn)
  remove(year, cn)
  
  #compile degree centrality / theorectical max into data frames  
  year <- c(1991:2017)
  #country 1: internal
  cn <- c(deg.cent_11_1991$centralization/deg.cent_11_1991$theoretical_max, 
          deg.cent_11_1992$centralization/deg.cent_11_1992$theoretical_max, 
          deg.cent_11_1993$centralization/deg.cent_11_1993$theoretical_max, 
          deg.cent_11_1994$centralization/deg.cent_11_1994$theoretical_max, 
          deg.cent_11_1995$centralization/deg.cent_11_1995$theoretical_max,
          deg.cent_11_1996$centralization/deg.cent_11_1996$theoretical_max,
          deg.cent_11_1997$centralization/deg.cent_11_1997$theoretical_max,
          deg.cent_11_1998$centralization/deg.cent_11_1998$theoretical_max,
          deg.cent_11_1999$centralization/deg.cent_11_1999$theoretical_max,
          deg.cent_11_2000$centralization/deg.cent_11_2000$theoretical_max,
          deg.cent_11_2001$centralization/deg.cent_11_2001$theoretical_max,
          deg.cent_11_2002$centralization/deg.cent_11_2002$theoretical_max,
          deg.cent_11_2003$centralization/deg.cent_11_2003$theoretical_max,
          deg.cent_11_2004$centralization/deg.cent_11_2004$theoretical_max,
          deg.cent_11_2005$centralization/deg.cent_11_2005$theoretical_max,
          deg.cent_11_2006$centralization/deg.cent_11_2006$theoretical_max,
          deg.cent_11_2007$centralization/deg.cent_11_2007$theoretical_max,
          deg.cent_11_2008$centralization/deg.cent_11_2008$theoretical_max,
          deg.cent_11_2009$centralization/deg.cent_11_2009$theoretical_max,
          deg.cent_11_2010$centralization/deg.cent_11_2010$theoretical_max,
          deg.cent_11_2011$centralization/deg.cent_11_2011$theoretical_max,
          deg.cent_11_2012$centralization/deg.cent_11_2012$theoretical_max,
          deg.cent_11_2013$centralization/deg.cent_11_2013$theoretical_max,
          deg.cent_11_2014$centralization/deg.cent_11_2014$theoretical_max,
          deg.cent_11_2015$centralization/deg.cent_11_2015$theoretical_max,
          deg.cent_11_2016$centralization/deg.cent_11_2016$theoretical_max,
          deg.cent_11_2017$centralization/deg.cent_11_2017$theoretical_max
  )
  degree_centrality11_max <- data.frame(year, cn)
  remove(
    deg.cent_11_1991,
    deg.cent_11_1992,
    deg.cent_11_1993,
    deg.cent_11_1994,
    deg.cent_11_1995,
    deg.cent_11_1996,
    deg.cent_11_1997,
    deg.cent_11_1998,
    deg.cent_11_1999,
    deg.cent_11_2000,
    deg.cent_11_2001,
    deg.cent_11_2002,
    deg.cent_11_2003,
    deg.cent_11_2004,
    deg.cent_11_2005,
    deg.cent_11_2006,
    deg.cent_11_2007,
    deg.cent_11_2008,
    deg.cent_11_2009,
    deg.cent_11_2010,
    deg.cent_11_2011,
    deg.cent_11_2012,
    deg.cent_11_2013,
    deg.cent_11_2014,
    deg.cent_11_2015,
    deg.cent_11_2016,
    deg.cent_11_2017)
  #country 1: external
  cn <- c(deg.cent_12_1991$centralization/deg.cent_12_1991$theoretical_max, 
          deg.cent_12_1992$centralization/deg.cent_12_1992$theoretical_max, 
          deg.cent_12_1993$centralization/deg.cent_12_1993$theoretical_max, 
          deg.cent_12_1994$centralization/deg.cent_12_1994$theoretical_max, 
          deg.cent_12_1995$centralization/deg.cent_12_1995$theoretical_max,
          deg.cent_12_1996$centralization/deg.cent_12_1996$theoretical_max,
          deg.cent_12_1997$centralization/deg.cent_12_1997$theoretical_max,
          deg.cent_12_1998$centralization/deg.cent_12_1998$theoretical_max,
          deg.cent_12_1999$centralization/deg.cent_12_1999$theoretical_max,
          deg.cent_12_2000$centralization/deg.cent_12_2000$theoretical_max,
          deg.cent_12_2001$centralization/deg.cent_12_2001$theoretical_max,
          deg.cent_12_2002$centralization/deg.cent_12_2002$theoretical_max,
          deg.cent_12_2003$centralization/deg.cent_12_2003$theoretical_max,
          deg.cent_12_2004$centralization/deg.cent_12_2004$theoretical_max,
          deg.cent_12_2005$centralization/deg.cent_12_2005$theoretical_max,
          deg.cent_12_2006$centralization/deg.cent_12_2006$theoretical_max,
          deg.cent_12_2007$centralization/deg.cent_12_2007$theoretical_max,
          deg.cent_12_2008$centralization/deg.cent_12_2008$theoretical_max,
          deg.cent_12_2009$centralization/deg.cent_12_2009$theoretical_max,
          deg.cent_12_2010$centralization/deg.cent_12_2010$theoretical_max,
          deg.cent_12_2011$centralization/deg.cent_12_2011$theoretical_max,
          deg.cent_12_2012$centralization/deg.cent_12_2012$theoretical_max,
          deg.cent_12_2013$centralization/deg.cent_12_2013$theoretical_max,
          deg.cent_12_2014$centralization/deg.cent_12_2014$theoretical_max,
          deg.cent_12_2015$centralization/deg.cent_12_2015$theoretical_max,
          deg.cent_12_2016$centralization/deg.cent_12_2016$theoretical_max,
          deg.cent_12_2017$centralization/deg.cent_12_2017$theoretical_max
  )
  degree_centrality12_max <- data.frame(year, cn)
  remove(
    deg.cent_12_1991,
    deg.cent_12_1992,
    deg.cent_12_1993,
    deg.cent_12_1994,
    deg.cent_12_1995,
    deg.cent_12_1996,
    deg.cent_12_1997,
    deg.cent_12_1998,
    deg.cent_12_1999,
    deg.cent_12_2000,
    deg.cent_12_2001,
    deg.cent_12_2002,
    deg.cent_12_2003,
    deg.cent_12_2004,
    deg.cent_12_2005,
    deg.cent_12_2006,
    deg.cent_12_2007,
    deg.cent_12_2008,
    deg.cent_12_2009,
    deg.cent_12_2010,
    deg.cent_12_2011,
    deg.cent_12_2012,
    deg.cent_12_2013,
    deg.cent_12_2014,
    deg.cent_12_2015,
    deg.cent_12_2016,
    deg.cent_12_2017)
  #country 2: internal
  cn <- c(deg.cent_22_1991$centralization/deg.cent_22_1991$theoretical_max, 
          deg.cent_22_1992$centralization/deg.cent_22_1992$theoretical_max, 
          deg.cent_22_1993$centralization/deg.cent_22_1993$theoretical_max, 
          deg.cent_22_1994$centralization/deg.cent_22_1994$theoretical_max, 
          deg.cent_22_1995$centralization/deg.cent_22_1995$theoretical_max,
          deg.cent_22_1996$centralization/deg.cent_22_1996$theoretical_max,
          deg.cent_22_1997$centralization/deg.cent_22_1997$theoretical_max,
          deg.cent_22_1998$centralization/deg.cent_22_1998$theoretical_max,
          deg.cent_22_1999$centralization/deg.cent_22_1999$theoretical_max,
          deg.cent_22_2000$centralization/deg.cent_22_2000$theoretical_max,
          deg.cent_22_2001$centralization/deg.cent_22_2001$theoretical_max,
          deg.cent_22_2002$centralization/deg.cent_22_2002$theoretical_max,
          deg.cent_22_2003$centralization/deg.cent_22_2003$theoretical_max,
          deg.cent_22_2004$centralization/deg.cent_22_2004$theoretical_max,
          deg.cent_22_2005$centralization/deg.cent_22_2005$theoretical_max,
          deg.cent_22_2006$centralization/deg.cent_22_2006$theoretical_max,
          deg.cent_22_2007$centralization/deg.cent_22_2007$theoretical_max,
          deg.cent_22_2008$centralization/deg.cent_22_2008$theoretical_max,
          deg.cent_22_2009$centralization/deg.cent_22_2009$theoretical_max,
          deg.cent_22_2010$centralization/deg.cent_22_2010$theoretical_max,
          deg.cent_22_2011$centralization/deg.cent_22_2011$theoretical_max,
          deg.cent_22_2012$centralization/deg.cent_22_2012$theoretical_max,
          deg.cent_22_2013$centralization/deg.cent_22_2013$theoretical_max,
          deg.cent_22_2014$centralization/deg.cent_22_2014$theoretical_max,
          deg.cent_22_2015$centralization/deg.cent_22_2015$theoretical_max,
          deg.cent_22_2016$centralization/deg.cent_22_2016$theoretical_max,
          deg.cent_22_2017$centralization/deg.cent_22_2017$theoretical_max
  )
  degree_centrality22_max <- data.frame(year, cn)
  remove(
    deg.cent_22_1991,
    deg.cent_22_1992,
    deg.cent_22_1993,
    deg.cent_22_1994,
    deg.cent_22_1995,
    deg.cent_22_1996,
    deg.cent_22_1997,
    deg.cent_22_1998,
    deg.cent_22_1999,
    deg.cent_22_2000,
    deg.cent_22_2001,
    deg.cent_22_2002,
    deg.cent_22_2003,
    deg.cent_22_2004,
    deg.cent_22_2005,
    deg.cent_22_2006,
    deg.cent_22_2007,
    deg.cent_22_2008,
    deg.cent_22_2009,
    deg.cent_22_2010,
    deg.cent_22_2011,
    deg.cent_22_2012,
    deg.cent_22_2013,
    deg.cent_22_2014,
    deg.cent_22_2015,
    deg.cent_22_2016,
    deg.cent_22_2017)
  #country 2: external 
  cn <- c(deg.cent_21_1991$centralization/deg.cent_21_1991$theoretical_max, 
          deg.cent_21_1992$centralization/deg.cent_21_1992$theoretical_max, 
          deg.cent_21_1993$centralization/deg.cent_21_1993$theoretical_max, 
          deg.cent_21_1994$centralization/deg.cent_21_1994$theoretical_max, 
          deg.cent_21_1995$centralization/deg.cent_21_1995$theoretical_max,
          deg.cent_21_1996$centralization/deg.cent_21_1996$theoretical_max,
          deg.cent_21_1997$centralization/deg.cent_21_1997$theoretical_max,
          deg.cent_21_1998$centralization/deg.cent_21_1998$theoretical_max,
          deg.cent_21_1999$centralization/deg.cent_21_1999$theoretical_max,
          deg.cent_21_2000$centralization/deg.cent_21_2000$theoretical_max,
          deg.cent_21_2001$centralization/deg.cent_21_2001$theoretical_max,
          deg.cent_21_2002$centralization/deg.cent_21_2002$theoretical_max,
          deg.cent_21_2003$centralization/deg.cent_21_2003$theoretical_max,
          deg.cent_21_2004$centralization/deg.cent_21_2004$theoretical_max,
          deg.cent_21_2005$centralization/deg.cent_21_2005$theoretical_max,
          deg.cent_21_2006$centralization/deg.cent_21_2006$theoretical_max,
          deg.cent_21_2007$centralization/deg.cent_21_2007$theoretical_max,
          deg.cent_21_2008$centralization/deg.cent_21_2008$theoretical_max,
          deg.cent_21_2009$centralization/deg.cent_21_2009$theoretical_max,
          deg.cent_21_2010$centralization/deg.cent_21_2010$theoretical_max,
          deg.cent_21_2011$centralization/deg.cent_21_2011$theoretical_max,
          deg.cent_21_2012$centralization/deg.cent_21_2012$theoretical_max,
          deg.cent_21_2013$centralization/deg.cent_21_2013$theoretical_max,
          deg.cent_21_2014$centralization/deg.cent_21_2014$theoretical_max,
          deg.cent_21_2015$centralization/deg.cent_21_2015$theoretical_max,
          deg.cent_21_2016$centralization/deg.cent_21_2016$theoretical_max,
          deg.cent_21_2017$centralization/deg.cent_21_2017$theoretical_max
  )
  degree_centrality21_max <- data.frame(year, cn)
  remove(
    deg.cent_21_1991,
    deg.cent_21_1992,
    deg.cent_21_1993,
    deg.cent_21_1994,
    deg.cent_21_1995,
    deg.cent_21_1996,
    deg.cent_21_1997,
    deg.cent_21_1998,
    deg.cent_21_1999,
    deg.cent_21_2000,
    deg.cent_21_2001,
    deg.cent_21_2002,
    deg.cent_21_2003,
    deg.cent_21_2004,
    deg.cent_21_2005,
    deg.cent_21_2006,
    deg.cent_21_2007,
    deg.cent_21_2008,
    deg.cent_21_2009,
    deg.cent_21_2010,
    deg.cent_21_2011,
    deg.cent_21_2012,
    deg.cent_21_2013,
    deg.cent_21_2014,
    deg.cent_21_2015,
    deg.cent_21_2016,
    deg.cent_21_2017)
  remove(year, cn)
  
  # Plot (not recommended)--------------------------------------------------------------------
  #plot  
  plot(degree_centrality22$year,degree_centrality22$cn, type = "o",col = "red", xlab = "Year", ylab = "Degree Centrality", 
       main = "Degree centrality", ylim = c(0.15,0.45)) 
  lines(degree_centrality21$year, degree_centrality21$cn, type = "o", col = "blue") 
  lines(degree_centrality11$year, degree_centrality11$cn, type = "o", col = "darkgreen")
  lines(degree_centrality12$year, degree_centrality12$cn, type = "o", col = "orange")
  legend(2010, 0.43, legend=c("C2:Internal", "C2:External", "C1:Internal", "C1:External"),
         col=c("red", "blue", "darkgreen", "orange"), lty=1:2, cex=0.8)  
  #plot max 
  plot(degree_centrality22_max$year,degree_centrality22_max$cn, type = "o",col = "red", xlab = "Year", 
       ylab = "Degree Centrality / Theorectical max", 
       main = "Degree centrality", ylim=c(0.000004,0.000012)) 
  lines(degree_centrality21_max$year, degree_centrality21_max$cn, type = "o", col = "blue") 
  lines(degree_centrality11_max$year, degree_centrality11_max$cn, type = "o", col = "darkgreen")
  lines(degree_centrality12_max$year, degree_centrality12_max$cn, type = "o", col = "orange")
  legend(2010, 0.000012, legend=c("C2:Internal", "C2:External", "C1:Internal", "C1:External"),
         col=c("red", "blue", "darkgreen", "orange"), lty=1:2, cex=0.8)  
  
# Newman-Garvin clustering ------------------------------------------------

  
  
  

  
  
  
# Compute and plot eigenvector centrality (needs revision)---------------------------------
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

  # Compile eigenvector centrality (needs revision)-----------------------------------------------------------------
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

  # Plot eigenvector centrality (needs revision)--------------------------------------------------------------------
#plot
plot(eigenvalue22,type = "o",col = "red", xlab = "Year", ylab = "Centrality", 
     main = "Eigenvector centrality (partial)") 
lines(eigenvalue21, type = "o", col = "blue") 
lines(eigenvalue11, type = "o", col = "darkgreen")
lines(eigenvalue12, type = "o", col = "orange")
legend(1991, 1.5, legend=c("C2:Internal", "C2:External", "C1:Internal", "C1:External"),
       col=c("red", "blue", "darkgreen", "orange"), lty=1:2, cex=0.6)

        