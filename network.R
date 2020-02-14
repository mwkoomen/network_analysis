
#NETWORK DATA


# Initialize --------------------------------------------------------------


#Load block, if these packages are not installed, you have to remove the pound sign and 
#run each code. 

#install.packages("dplyr")
#install.packages("sqldf")
#install.packages("repmis")
#install.packages("RCurl")
#install.packages("igraph")

#Then load the packages into the current working session:
library(igraph)
library(dplyr)
library(sqldf)

#!!!!! Run seperately: this will prompt you with a quesiton if you want to load 
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
source_data("https://github.com/mwkoomen/network_analysis/blob/master/Network_Data.rdata?raw=true")

#The data is stored as the partial set of data, i.e. null connections are missing
par_set <- `Network Data`

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
                  select *, 0 from dif
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
  plot(pc1c1,type = "o",col = "red", xlab = "Year", ylab = "Intensity", 
       main = "Network Intensity (partial)") 
  lines(pc1c2, type = "o", col = "blue") 
  lines(pc2c2, type = "o", col = "darkgreen")
  lines(pc2c1, type = "o", col = "orange")
  legend(1991, 45000, legend=c("C1:Internal", "C1:External", "C2:Internal", "C2:External"),
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
  plot(fc1c1,type = "o",col = "red", xlab = "Year", ylab = "Intensity", 
       main = "Network Intensity (full)") 
  lines(fc1c2, type = "o", col = "blue") 
  lines(fc2c2, type = "o", col = "darkgreen")
  lines(fc2c1, type = "o", col = "orange")
  legend(1991, 8500, legend=c("C1:Internal", "C1:External", "C2:Internal", "C2:External"),
         col=c("red", "blue", "darkgreen", "orange"), lty=1:2, cex=0.6)


# Compute and plot eigenvector centrality ---------------------------------

#This sections computes eigenvector centrality 
#Currently it is only an example of the internal link intensity for country 1 in 1991
#To Do: expand for all years and combination of countries 
pc1c1_91 <- sqldf("select * from par_set where Year = 1991 and Country_A = 1 and Country_B = 1") 
graph_11_91 <- graph(c(pc1c1_91$Region_A, pc1c1_91$Region_B)) 
EV_pc1c1_91 <- eigen_centrality(graph_11_91, weights = pc1c1_91$Intensity_norm)
  
  