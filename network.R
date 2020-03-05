####################################################
#NETWORK DATA ANALYSIS
####################################################

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
#install.packages("ggplot2")
#install.packages("sna")
#install.packages("qgraph")

#Then load the packages into the current working session:
library(igraph)
library(dplyr)
library(sqldf)
library(readstata13)
library(CINNA)
library(ggplot2)
library(sna)
library(qgraph)
library(data.table)

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


# TEST duplicate rows [NO DUPLICATE ROWS] -------------------------------------------------
par_set$edge1 <- paste(par_set$Region_A, "-", par_set$Region_B, sep="")
par_set$edge2 <- paste(par_set$Region_B, "-", par_set$Region_A, sep="")

par_set <- par_set %>%
  mutate(autolink=ifelse(Country_A==Country_B & Region_A==Region_B, 1, 0)) %>%
  mutate(border=ifelse(Country_A!=Country_B, 1, 0))

test_doubles <- data.frame()
for (c1 in 1:2){
  for (c2 in 1:2){
    for (y in 1991:2017) {
      s <- paste("select count(*) 
                 from par_set where Year=",y,
                 " and Country_A=",c1,
                 " and Country_B=",c2, 
                 " and autolink=0", sep="")
      r <- as.numeric(sqldf(s))
      f <- paste("select edge1 from par_set where Year=",y,
                 " and Country_A=",c1,
                 " and Country_B=",c2,
                 " and autolink=0
                 except 
                 select edge2 from par_set where Year=",y,
                 " and Country_A=",c1,
                 " and Country_B=",c2,
                 " and autolink=0", sep="")
      t <- sqldf(f)
      z <- data.frame(Year=y, CountryA=c1, no_doubles=nrow(t)==r)
      test_doubles <- rbind(test_doubles, z)
    }
  }    
}
remove(c1,c2,f,y,r,s,t,z)
paste("Data has duplicate row: ", FALSE %in% test_doubles$no_doubles)  

# Create full set ---------------------------------------------------------
regions <- data.frame()
for (a in 188:2){
  for (i in 1:189){
    c <- expand.grid(i,a)
    regions <- rbind(regions, c)
  }
}

# [OLD] Create full set ---------------------------------------------------------

#This block creates a full data set that contains all the missing links 

#This code creates a full (theoretical) set of dimensions 27*189*2*189*2
comp_set <- expand.grid(1991:2017, 1:189, 1:2, 1:189, 1:2)

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
         col=c("red", "blue", "darkgreen", "orange"), lty=1, cex=0.8, bty="n")

  
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
         col=c("red", "blue", "darkgreen", "orange"), lty=1, cex=0.8, bty="n")
# Count zeros --------------------------------------------
  
  zc1c1 <- full_set %>%
    filter(Country_A==1, Country_B==1 & Intensity==0) %>%
    group_by(Year) %>%
    count(Intensity)
  zc1c1$n <- round(zc1c1$n/189)
  zc1c2 <- full_set %>%
    filter(Country_A==1, Country_B==2 & Intensity==0) %>%
    group_by(Year) %>%
    count(Intensity)
  zc1c2$n <- round(zc1c2$n/189)
  zc2c2 <- full_set %>%
    filter(Country_A==2, Country_B==2 & Intensity==0) %>%
    group_by(Year) %>%
    count(Intensity)
  zc2c2$n <- round(zc2c2$n/189)
  zc2c1 <- full_set %>%
    filter(Country_A==2, Country_B==1 & Intensity==0) %>%
    group_by(Year) %>%
    count(Intensity)
  zc2c1$n <- round(zc2c1$n/189)
  zc1c1$link <- "Country 1: Internal"
  zc1c2$link <- "Country 1: External"
  zc2c2$link <- "Country 2: Internal"
  zc2c1$link <- "Country 2: External"
  zero.vertices <- rbind(zc1c1, zc1c2, zc2c2, zc2c1)  
  remove(zc1c1, zc1c2, zc2c2, zc2c1)
  
  p <- ggplot(zero.vertices, aes(Year, n))
  p +  ggtitle("# vertices with no edges") + 
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_smooth(aes(colour = factor(link)))
  
    
  # plot(zc1c1$Year, zc1c1$n, type = "o", col = "darkgreen", xlab = "Year", ylab = "# vertices", 
  #      main = "# vertices with no edges", ylim=c(140,189)) 
  # lines(zc1c2$Year, zc1c2$n, type = "o", col = "orange") 
  # lines(zc2c2$Year, zc2c2$n, type = "o", col = "red")
  # lines(zc2c1$Year, zc2c1$n, type = "o", col = "blue")
  # legend(1991,155, legend=c("C2:Internal", "C2:External", "C1:Internal", "C1:External"),
  #        col=c("red", "blue", "darkgreen", "orange"), lty=1, cex=0.8, bty="n")


# Construct adjacency matrix per year, country, links [directed]------------------------------------------------
  nodes <- data.frame(ID=1:189)
  
  #Country 1: Internal links 
  for (i in 1991:2017) {
    s <- paste("select * from par_set where Year = ", i, " and Country_A = 1 and Country_B = 1")
    x <- sqldf(s)
    n <- paste("net_11_", i, sep="")
    links <- data.frame(RegionA=x$Region_A, RegionB=x$Region_B, w1=x$Intensity, w2=x$Intensity_norm)
    assign(n, graph_from_data_frame(links, nodes, directed = T))
  }
  #Country 1: External links  
  for (i in 1991:2017) {
    s <- paste("select * from par_set where Year = ", i, " and Country_A = 1 and Country_B = 2")
    x <- sqldf(s)
    n <- paste("net_12_", i, sep="")
    links <- data.frame(RegionA=x$Region_A, RegionB=x$Region_B, w1=x$Intensity, w2=x$Intensity_norm)
    assign(n, graph_from_data_frame(links, nodes, directed = T))
  }
  #Country 2: Internal links  
  for (i in 1991:2017) {
    s <- paste("select * from par_set where Year = ", i, " and Country_A = 2 and Country_B = 2")
    n <- paste("net_22_", i, sep="")
    links <- data.frame(RegionA=x$Region_A, RegionB=x$Region_B, w1=x$Intensity, w2=x$Intensity_norm)
    assign(n, graph_from_data_frame(links, nodes, directed = T))
  }
  #Country 2: External links  
  for (i in 1991:2017) {
    s <- paste("select * from par_set where Year = ", i, " and Country_A = 2 and Country_B = 1")
    x <- sqldf(s)
    n <- paste("net_21_", i, sep="")
    links <- data.frame(RegionA=x$Region_A, RegionB=x$Region_B, w1=x$Intensity, w2=x$Intensity_norm)
    assign(n, graph_from_data_frame(links, nodes, directed = T))
  }
  remove(x,i,s,links, nodes, n)

  #example plot with no loops or unconnected regions
  par(old.par)
  plot(delete.vertices(simplify(net_12_2017), degree(net_12_2017)==0), 
       layout=layout_with_fr, vertex.label = NA, vertex.size=6)  

# Construct adjacency matrix per year, country, links [undirected]------------------------------------------------
  nodes <- data.frame(ID=1:189)
  
  #Country 1: Internal links 
  for (i in 1991:2017) {
    s <- paste("select * from par_set where Year = ", i, " and Country_A = 1 and Country_B = 1")
    x <- sqldf(s)
    n <- paste("net_11_", i, sep="")
    links <- data.frame(RegionA=x$Region_A, RegionB=x$Region_B, w1=x$Intensity, w2=x$Intensity_norm)
    assign(n, graph_from_data_frame(links, nodes, directed = F))
  }
  #Country 1: External links  
  for (i in 1991:2017) {
    s <- paste("select * from par_set where Year = ", i, " and Country_A = 1 and Country_B = 2")
    x <- sqldf(s)
    n <- paste("net_12_", i, sep="")
    links <- data.frame(RegionA=x$Region_A, RegionB=x$Region_B, w1=x$Intensity, w2=x$Intensity_norm)
    assign(n, graph_from_data_frame(links, nodes, directed = F))
  }
  #Country 2: Internal links  
  for (i in 1991:2017) {
    s <- paste("select * from par_set where Year = ", i, " and Country_A = 2 and Country_B = 2")
    n <- paste("net_22_", i, sep="")
    links <- data.frame(RegionA=x$Region_A, RegionB=x$Region_B, w1=x$Intensity, w2=x$Intensity_norm)
    assign(n, graph_from_data_frame(links, nodes, directed = F))
  }
  #Country 2: External links  
  for (i in 1991:2017) {
    s <- paste("select * from par_set where Year = ", i, " and Country_A = 2 and Country_B = 1")
    x <- sqldf(s)
    n <- paste("net_21_", i, sep="")
    links <- data.frame(RegionA=x$Region_A, RegionB=x$Region_B, w1=x$Intensity, w2=x$Intensity_norm)
    assign(n, graph_from_data_frame(links, nodes, directed = F))
  }
  remove(x,i,s,links, nodes, n)
  
  #example plot with no loops or unconnected regions
  par(old.par)
  plot(delete.vertices(simplify(net_12_2017), degree(net_12_2017)==0), 
       layout=layout_with_fr, vertex.label = NA, vertex.size=6)  
  
  
# Degree distibution plots [not normalized] [unweighted] (works) ------------------------------------------------

#C1 external links  
  deg12_1 <- degree(net_12_1991, mode="all")
  deg12_2 <- degree(net_12_1996, mode="all")
  deg12_3 <- degree(net_12_2001, mode="all")
  deg12_4 <- degree(net_12_2006, mode="all")
  deg12_5 <- degree(net_12_2011, mode="all")
  deg12_6 <- degree(net_12_2017, mode="all")  
  deg.dist12_1 <- degree_distribution(net_12_1991, cumulative=T, mode="all")
  deg.dist12_2 <- degree_distribution(net_12_1996, cumulative=T, mode="all")
  deg.dist12_3 <- degree_distribution(net_12_2001, cumulative=T, mode="all")
  deg.dist12_4 <- degree_distribution(net_12_2006, cumulative=T, mode="all")
  deg.dist12_5 <- degree_distribution(net_12_2011, cumulative=T, mode="all")
  deg.dist12_6 <- degree_distribution(net_12_2017, cumulative=T, mode="all")  
#C1 internal links  
  deg11_1 <- degree(net_11_1991, mode="all")
  deg11_2 <- degree(net_11_1996, mode="all")
  deg11_3 <- degree(net_11_2001, mode="all")
  deg11_4 <- degree(net_11_2006, mode="all")
  deg11_5 <- degree(net_11_2011, mode="all")
  deg11_6 <- degree(net_11_2017, mode="all")  
  deg.dist11_1 <- degree_distribution(net_11_1991, cumulative=T, mode="all")
  deg.dist11_2 <- degree_distribution(net_11_1996, cumulative=T, mode="all")
  deg.dist11_3 <- degree_distribution(net_11_2001, cumulative=T, mode="all")
  deg.dist11_4 <- degree_distribution(net_11_2006, cumulative=T, mode="all")
  deg.dist11_5 <- degree_distribution(net_11_2011, cumulative=T, mode="all")
  deg.dist11_6 <- degree_distribution(net_11_2017, cumulative=T, mode="all")
#C2 external links  
  deg21_1 <- degree(net_21_1991, mode="all")
  deg21_2 <- degree(net_21_1996, mode="all")
  deg21_3 <- degree(net_21_2001, mode="all")
  deg21_4 <- degree(net_21_2006, mode="all")
  deg21_5 <- degree(net_21_2011, mode="all")
  deg21_6 <- degree(net_21_2017, mode="all")  
  deg.dist21_1 <- degree_distribution(net_21_1991, cumulative=T, mode="all")
  deg.dist21_2 <- degree_distribution(net_21_1996, cumulative=T, mode="all")
  deg.dist21_3 <- degree_distribution(net_21_2001, cumulative=T, mode="all")
  deg.dist21_4 <- degree_distribution(net_21_2006, cumulative=T, mode="all")
  deg.dist21_5 <- degree_distribution(net_21_2011, cumulative=T, mode="all")
  deg.dist21_6 <- degree_distribution(net_21_2017, cumulative=T, mode="all")  
#C2 internal links  
  deg22_1 <- degree(net_22_1991, mode="all")
  deg22_2 <- degree(net_22_1996, mode="all")
  deg22_3 <- degree(net_22_2001, mode="all")
  deg22_4 <- degree(net_22_2006, mode="all")
  deg22_5 <- degree(net_22_2011, mode="all")
  deg22_6 <- degree(net_22_2017, mode="all")  
  deg.dist22_1 <- degree_distribution(net_22_1991, cumulative=T, mode="all")
  deg.dist22_2 <- degree_distribution(net_22_1996, cumulative=T, mode="all")
  deg.dist22_3 <- degree_distribution(net_22_2001, cumulative=T, mode="all")
  deg.dist22_4 <- degree_distribution(net_22_2006, cumulative=T, mode="all")
  deg.dist22_5 <- degree_distribution(net_22_2011, cumulative=T, mode="all")
  deg.dist22_6 <- degree_distribution(net_22_2017, cumulative=T, mode="all")
  
#plot country 1  
  par(mfrow=c(2,3), oma = c(0,2,2,0))
  plot(x=0:max(deg11_1), y=1-deg.dist11_1, cex=1.2, col=rgb(0,0.4,0.8,alpha=0.4), ylim=c(0,1), xlim=c(0,130), 
        xlab="1991", ylab="Cumulative Frequency", type = "h")
  lines(x=0:max(deg12_1), y=1-deg.dist12_1, cex=1.2, col=rgb(1,0.2,0.2, alpha=0.4), ylim=c(0,1), xlim=c(0,130), 
        type="h")
  legend(1,200,legend=c("Internal", "External"),col=c("blue", "red"), lty=1, cex=0.9, bty="n") 
  plot(x=0:max(deg11_2), y=1-deg.dist11_2, cex=1.2, col=rgb(0,0.4,0.8,alpha=0.4), ylim=c(0,1), xlim=c(0,130), 
        xlab="1996", ylab="", type = "h")
  lines(x=0:max(deg12_2), y=1-deg.dist12_2, cex=1.2, col=rgb(1,0.2,0.2, alpha=0.4), ylim=c(0,1), xlim=c(0,130), 
         type="h")
  plot(x=0:max(deg11_3), y=1-deg.dist11_3, cex=1.2, col=rgb(0,0.4,0.8,alpha=0.4), ylim=c(0,1), xlim=c(0,130), 
        xlab="2001", ylab="", type = "h")
  lines(x=0:max(deg12_3), y=1-deg.dist12_3, cex=1.2, col=rgb(1,0.2,0.2, alpha=0.4), ylim=c(0,1), xlim=c(0,130), 
         type="h")
  plot( x=0:max(deg11_4), y=1-deg.dist11_4, cex=1.2, col=rgb(0,0.4,0.8,alpha=0.4), ylim=c(0,1), xlim=c(0,130), 
        xlab="2006", ylab="Cumulative Frequency", type = "h")
  lines(x=0:max(deg12_4), y=1-deg.dist12_4, cex=1.2, col=rgb(1,0.2,0.2, alpha=0.4), ylim=c(0,1), xlim=c(0,130), 
         type="h") 
  plot( x=0:max(deg11_5), y=1-deg.dist11_5, cex=1.2, col=rgb(0,0.4,0.8,alpha=0.4), ylim=c(0,1), xlim=c(0,130), 
        xlab="2011", ylab="", type = "h")
  lines(x=0:max(deg12_5), y=1-deg.dist12_5, cex=1.2, col=rgb(1,0.2,0.2, alpha=0.4), ylim=c(0,1), xlim=c(0,130), 
        type="h")
  plot(x=0:max(deg11_6), y=1-deg.dist11_6, cex=1.2, col=rgb(0,0.4,0.8,alpha=0.4), ylim=c(0,1), xlim=c(0,130), 
        xlab="2017", ylab="", type = "h")
  lines(x=0:max(deg12_6), y=1-deg.dist12_6, cex=1.2, col=rgb(1,0.2,0.2, alpha=0.4), ylim=c(0,1), xlim=c(0,130), 
        type="h")
  mtext("Country 1: Degree centraliy distribution", outer=TRUE, cex=1, font=2)

  #plot country 2  
  par(mfrow=c(2,3), oma = c(0,2,2,0))
  plot( x=0:max(deg22_1), y=1-deg.dist22_1, cex=1.2, col=rgb(0,0.4,0.8,alpha=0.4), ylim=c(0,1), xlim=c(0,130), 
        xlab="1991", ylab="Cumulative Frequency", type = "h")
  lines( x=0:max(deg21_1), y=1-deg.dist21_1, cex=1.2, col=rgb(1,0.2,0.2, alpha=0.4), ylim=c(0,1), xlim=c(0,130), 
        type="h")
  plot( x=0:max(deg22_2), y=1-deg.dist22_2,  cex=1.2, col=rgb(0,0.4,0.8,alpha=0.4), ylim=c(0,1), xlim=c(0,130), 
        xlab="1996", ylab="", type = "h")
  lines( x=0:max(deg21_2), y=1-deg.dist21_2, cex=1.2, col=rgb(1,0.2,0.2, alpha=0.4), ylim=c(0,1), xlim=c(0,130), 
        type="h")
  plot( x=0:max(deg22_3), y=1-deg.dist22_3, cex=1.2, col=rgb(0,0.4,0.8,alpha=0.4), ylim=c(0,1), xlim=c(0,130), 
        xlab="2001", ylab="", type = "h")
  lines( x=0:max(deg21_3), y=1-deg.dist21_3, cex=1.2, col=rgb(1,0.2,0.2, alpha=0.4), ylim=c(0,1), xlim=c(0,130), 
        type="h")
  plot( x=0:max(deg22_4), y=1-deg.dist22_4, cex=1.2, col=rgb(0,0.4,0.8,alpha=0.4), ylim=c(0,1), xlim=c(0,130), 
        xlab="2006", ylab="Cumulative Frequency", type = "h")
  lines( x=0:max(deg21_4), y=1-deg.dist21_4, cex=1.2, col=rgb(1,0.2,0.2, alpha=0.4), ylim=c(0,1), xlim=c(0,130), 
        type="h") 
  plot( x=0:max(deg22_5), y=1-deg.dist22_5, cex=1.2, col=rgb(0,0.4,0.8,alpha=0.4), ylim=c(0,1), xlim=c(0,130), 
        xlab="2011", ylab="", type = "h")
  lines( x=0:max(deg21_5), y=1-deg.dist21_5, cex=1.2, col=rgb(1,0.2,0.2, alpha=0.4), ylim=c(0,1), xlim=c(0,130), 
        type="h")
  plot( x=0:max(deg22_6), y=1-deg.dist22_6, cex=1.2, col=rgb(0,0.4,0.8,alpha=0.4), ylim=c(0,1), xlim=c(0,130), 
        xlab="2017", ylab="", type = "h")
  lines( x=0:max(deg21_6), y=1-deg.dist21_6, cex=1.2, col=rgb(1,0.2,0.2, alpha=0.4), ylim=c(0,1), xlim=c(0,130), 
        type="h")
  mtext("Country 2: Degree centraliy distribution", outer=TRUE, cex=1, font=2)
  
centralities <- proper_centralities(net_11_1991)
# Compute and plot degree centrality [not normalized] [unweighted] (works) --------------------------------------
nodes <- data.frame(ID=1:189)

#Country 1: Internal links 
for (i in 1991:2017) {
  s <- paste("select * from par_set where Year = ", i, " and Country_A = 1 and Country_B = 1")
  x <- sqldf(s)
  links <- data.frame(RegionA=x$Region_A, RegionB=x$Region_B)
  net <- graph_from_data_frame(links, nodes, directed = TRUE)
  n <- paste("deg.cent_11_", i, sep="")
  assign(n, centr_degree(net, mode = "all", loops=TRUE, normalized = F))
}
#Country 1: External links  
for (i in 1991:2017) {
  s <- paste("select * from par_set where Year = ", i, " and Country_A = 1 and Country_B = 2")
  x <- sqldf(s)
  links <- data.frame(RegionA=x$Region_A, RegionB=x$Region_B)
  net <- graph_from_data_frame(links, nodes, directed = TRUE)
  n <- paste("deg.cent_12_", i, sep="")
  assign(n, centr_degree(net, mode = "all", loops=TRUE, normalized = F))
}
#Country 2: Internal links  
for (i in 1991:2017) {
  s <- paste("select * from par_set where Year = ", i, " and Country_A = 2 and Country_B = 2")
  x <- sqldf(s)
  links <- data.frame(RegionA=x$Region_A, RegionB=x$Region_B)
  net <- graph_from_data_frame(links, nodes, directed = TRUE)
  n <- paste("deg.cent_22_", i, sep="")
  assign(n, centr_degree(net, mode = "all", loops=TRUE, normalized = F))
}
#Country 2: External links  
for (i in 1991:2017) {
  s <- paste("select * from par_set where Year = ", i, " and Country_A = 2 and Country_B = 1")
  x <- sqldf(s)
  links <- data.frame(RegionA=x$Region_A, RegionB=x$Region_B)
  net <- graph_from_data_frame(links, nodes, directed = TRUE)
  n <- paste("deg.cent_21_", i, sep="")
  assign(n, centr_degree(net, mode = "all", loops=TRUE, normalized = F))
}
remove(x,i,s,links, nodes, net, n)

# Compute and plot degree centrality [normalized] [unweighted] (works) --------------------------------------
  nodes <- data.frame(ID=1:189)

  #Country 1: Internal links 
  for (i in 1991:2017) {
    s <- paste("select * from par_set where Year = ", i, " and Country_A = 1 and Country_B = 1")
    x <- sqldf(s)
    links <- data.frame(RegionA=x$Region_A, RegionB=x$Region_B)
    net <- graph_from_data_frame(links, nodes, directed = F)
    n <- paste("deg.cent_11_", i, sep="")
    assign(n, centr_degree(net, mode = "all", loops=T, normalized = T))
  }
  #Country 1: External links  
  for (i in 1991:2017) {
    s <- paste("select * from par_set where Year = ", i, " and Country_A = 1 and Country_B = 2")
    x <- sqldf(s)
    links <- data.frame(RegionA=x$Region_A, RegionB=x$Region_B)
    net <- graph_from_data_frame(links, nodes, directed = F)
    n <- paste("deg.cent_12_", i, sep="")
    assign(n, centr_degree(net, mode = "all", loops=T, normalized = T))
  }
  #Country 2: Internal links  
  for (i in 1991:2017) {
    s <- paste("select * from par_set where Year = ", i, " and Country_A = 2 and Country_B = 2")
    x <- sqldf(s)
    links <- data.frame(RegionA=x$Region_A, RegionB=x$Region_B)
    net <- graph_from_data_frame(links, nodes, directed = F)
    n <- paste("deg.cent_22_", i, sep="")
    assign(n, centr_degree(net, mode = "all", loops=T, normalized = T))
  }
  #Country 2: External links  
  for (i in 1991:2017) {
    s <- paste("select * from par_set where Year = ", i, " and Country_A = 2 and Country_B = 1")
    x <- sqldf(s)
    links <- data.frame(RegionA=x$Region_A, RegionB=x$Region_B)
    net <- graph_from_data_frame(links, nodes, directed = F)
    n <- paste("deg.cent_21_", i, sep="")
    assign(n, centr_degree(net, mode = "all", loops=T, normalized = T))
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


  # Plot Degree Centrality [not normalized] (works) --------------------------------------------------------------------
  par(old.par)
  plot(degree_centrality22$year,degree_centrality22$cn, type = "o",col = "red", xlab = "Year", ylab = "Degree Centrality (Mean)", 
       main = "Degree centrality", ylim=c(4000,18000)) 
  lines(degree_centrality21$year, degree_centrality21$cn, type = "o", col = "blue") 
  lines(degree_centrality11$year, degree_centrality11$cn, type = "o", col = "darkgreen")
  lines(degree_centrality12$year, degree_centrality12$cn, type = "o", col = "orange")
  legend(2008, 0.25, legend=c("C2:Internal", "C2:External", "C1:Internal", "C1:External"),
         col=c("red", "blue", "darkgreen", "orange"), lty=1:2, cex=0.8)    
  
  # Plot Degree Centrality [normalized] (works) --------------------------------------------------------------------

  #plot  
  par(old.par)
  plot(degree_centrality22$year,degree_centrality22$cn, type = "o",col = "red", xlab = "Year", ylab = "Degree Centrality (Mean)", 
       main = "Degree centrality", ylim = c(0.15,0.45)) 
  lines(degree_centrality21$year, degree_centrality21$cn, type = "o", col = "blue") 
  lines(degree_centrality11$year, degree_centrality11$cn, type = "o", col = "darkgreen")
  lines(degree_centrality12$year, degree_centrality12$cn, type = "o", col = "orange")
  legend(2008, 0.25, legend=c("C2:Internal", "C2:External", "C1:Internal", "C1:External"),
         col=c("red", "blue", "darkgreen", "orange"), lty=1:2, cex=0.8)  
 
# Compute and plot degree centrality [not normalized] [weighted] (works) --------------------------------------
  nodes <- data.frame(ID=1:189)
  
  #Country 1: Internal links 
  for (i in 1991:2017) {
    s <- paste("select * from par_set where Year = ", i, " and Country_A = 1 and Country_B = 1")
    x <- sqldf(s)
    links <- data.frame(RegionA=x$Region_A, RegionB=x$Region_B, w=x$Intensity_norm)
    net <- graph_from_data_frame(links, nodes, directed = TRUE)
    n <- paste("deg.cent_11_", i, sep="")
    assign(n, strength(net, mode = "all", loops=TRUE, weights = E(net)$w))
  }
  #Country 1: External links  
  for (i in 1991:2017) {
    s <- paste("select * from par_set where Year = ", i, " and Country_A = 1 and Country_B = 2")
    x <- sqldf(s)
    links <- data.frame(RegionA=x$Region_A, RegionB=x$Region_B, w=x$Intensity_norm)
    net <- graph_from_data_frame(links, nodes, directed = TRUE)
    n <- paste("deg.cent_12_", i, sep="")
    assign(n, strength(net, mode = "all", loops=TRUE, weights = E(net)$w))
  }
  #Country 2: Internal links  
  for (i in 1991:2017) {
    s <- paste("select * from par_set where Year = ", i, " and Country_A = 2 and Country_B = 2")
    x <- sqldf(s)
    links <- data.frame(RegionA=x$Region_A, RegionB=x$Region_B, w=x$Intensity_norm)
    net <- graph_from_data_frame(links, nodes, directed = TRUE)
    n <- paste("deg.cent_22_", i, sep="")
    assign(n, strength(net, mode = "all", loops=TRUE, weights = E(net)$w))
  }
  #Country 2: External links  
  for (i in 1991:2017) {
    s <- paste("select * from par_set where Year = ", i, " and Country_A = 2 and Country_B = 1")
    x <- sqldf(s)
    links <- data.frame(RegionA=x$Region_A, RegionB=x$Region_B, w=x$Intensity_norm)
    net <- graph_from_data_frame(links, nodes, directed = TRUE)
    n <- paste("deg.cent_21_", i, sep="")
    assign(n, strength(net, mode = "all", loops=TRUE, weights = E(net)$w))
  }
  remove(x,i,s,links, nodes, net, n)
  
  
  # Compile Degree Centrality [weighted] (works) -----------------------------------------------
  year <- c(1991:2017)
  #country 1: internal
  cn <- c(mean(deg.cent_11_1991),
          mean(deg.cent_11_1992),
          mean(deg.cent_11_1993),
          mean(deg.cent_11_1994),
          mean(deg.cent_11_1995),
          mean(deg.cent_11_1996),
          mean(deg.cent_11_1997),
          mean(deg.cent_11_1998),
          mean(deg.cent_11_1999),
          mean(deg.cent_11_2000),
          mean(deg.cent_11_2001),
          mean(deg.cent_11_2002),
          mean(deg.cent_11_2003),
          mean(deg.cent_11_2004),
          mean(deg.cent_11_2005),
          mean(deg.cent_11_2006),
          mean(deg.cent_11_2007),
          mean(deg.cent_11_2008),
          mean(deg.cent_11_2009),
          mean(deg.cent_11_2010),
          mean(deg.cent_11_2011),
          mean(deg.cent_11_2012),
          mean(deg.cent_11_2013),
          mean(deg.cent_11_2014),
          mean(deg.cent_11_2015),
          mean(deg.cent_11_2016),
          mean(deg.cent_11_2017)
  )
  degree_centrality11 <- data.frame(year, cn)
  #country 1: external
  cn <- c(mean(deg.cent_12_1991),
          mean(deg.cent_12_1992),
          mean(deg.cent_12_1993),
          mean(deg.cent_12_1994),
          mean(deg.cent_12_1995),
          mean(deg.cent_12_1996),
          mean(deg.cent_12_1997),
          mean(deg.cent_12_1998),
          mean(deg.cent_12_1999),
          mean(deg.cent_12_2000),
          mean(deg.cent_12_2001),
          mean(deg.cent_12_2002),
          mean(deg.cent_12_2003),
          mean(deg.cent_12_2004),
          mean(deg.cent_12_2005),
          mean(deg.cent_12_2006),
          mean(deg.cent_12_2007),
          mean(deg.cent_12_2008),
          mean(deg.cent_12_2009),
          mean(deg.cent_12_2010),
          mean(deg.cent_12_2011),
          mean(deg.cent_12_2012),
          mean(deg.cent_12_2013),
          mean(deg.cent_12_2014),
          mean(deg.cent_12_2015),
          mean(deg.cent_12_2016),
          mean(deg.cent_12_2017)
  )
  degree_centrality12 <- data.frame(year, cn)
  #country 2: internal
  cn <- c(mean(deg.cent_22_1991),
          mean(deg.cent_22_1992),
          mean(deg.cent_22_1993),
          mean(deg.cent_22_1994),
          mean(deg.cent_22_1995),
          mean(deg.cent_22_1996),
          mean(deg.cent_22_1997),
          mean(deg.cent_22_1998),
          mean(deg.cent_22_1999),
          mean(deg.cent_22_2000),
          mean(deg.cent_22_2001),
          mean(deg.cent_22_2002),
          mean(deg.cent_22_2003),
          mean(deg.cent_22_2004),
          mean(deg.cent_22_2005),
          mean(deg.cent_22_2006),
          mean(deg.cent_22_2007),
          mean(deg.cent_22_2008),
          mean(deg.cent_22_2009),
          mean(deg.cent_22_2010),
          mean(deg.cent_22_2011),
          mean(deg.cent_22_2012),
          mean(deg.cent_22_2013),
          mean(deg.cent_22_2014),
          mean(deg.cent_22_2015),
          mean(deg.cent_22_2016),
          mean(deg.cent_22_2017)
  )
  degree_centrality22 <- data.frame(year, cn)
  #country 2: external 
  cn <- c(mean(deg.cent_21_1991),
          mean(deg.cent_21_1992),
          mean(deg.cent_21_1993),
          mean(deg.cent_21_1994),
          mean(deg.cent_21_1995),
          mean(deg.cent_21_1996),
          mean(deg.cent_21_1997),
          mean(deg.cent_21_1998),
          mean(deg.cent_21_1999),
          mean(deg.cent_21_2000),
          mean(deg.cent_21_2001),
          mean(deg.cent_21_2002),
          mean(deg.cent_21_2003),
          mean(deg.cent_21_2004),
          mean(deg.cent_21_2005),
          mean(deg.cent_21_2006),
          mean(deg.cent_21_2007),
          mean(deg.cent_21_2008),
          mean(deg.cent_21_2009),
          mean(deg.cent_21_2010),
          mean(deg.cent_21_2011),
          mean(deg.cent_21_2012),
          mean(deg.cent_21_2013),
          mean(deg.cent_21_2014),
          mean(deg.cent_21_2015),
          mean(deg.cent_21_2016),
          mean(deg.cent_21_2017)
  )
  degree_centrality21 <- data.frame(year, cn)
  remove(year, cn)
  
  # Plot Degree Centrality [non normalized] [weighted] (works) --------------------------------------------------------------------
  par(old.par)
  plot(degree_centrality22$year,degree_centrality22$cn, type = "o",col = "red", xlab = "Year", ylab = "Degree Centrality (Mean)", 
       main = "Degree centrality [weighted]", ylim=c(0,0.12)) 
  lines(degree_centrality21$year, degree_centrality21$cn, type = "o", col = "blue") 
  lines(degree_centrality11$year, degree_centrality11$cn, type = "o", col = "darkgreen")
  lines(degree_centrality12$year, degree_centrality12$cn, type = "o", col = "orange")
  legend(1991, 0.12, legend=c("C2:Internal", "C2:External", "C1:Internal", "C1:External"),
         col=c("red", "blue", "darkgreen", "orange"), lty=1, cex=0.8, bty="n")  
  
# Construct nework graphs / adjacency matrix per year, country (works)------------------------------------------------
  nodes <- data.frame(ID=1:378)
  
  #Country 1:  
  for (i in 1991:2017) {
    s <- paste("select * from par_set where Year = ", i, " and Country_A = 1 and (
               Country_B = 1 or Country_B = 2)")
    x <- sqldf(s)
    x <- x %>%
      mutate(Region_B = ifelse(Country_B==2, Region_B+189, Region_B))
    n <- paste("net_1_", i, sep="")
    links <- data.frame(RegionA=x$Region_A, RegionB=x$Region_B)
    assign(n, graph_from_data_frame(links, nodes, directed = FALSE))
  }
  #Country 2:   
  for (i in 1991:2017) {
    s <- paste("select * from par_set where Year = ", i, " and Country_A = 2 and (
               Country_B = 1 or Country_B = 2)")
    x <- sqldf(s)
    x <- x %>%
      mutate(Region_B = ifelse(Country_B==1, Region_B+189, Region_B))
    n <- paste("net_2_", i, sep="")
    links <- data.frame(RegionA=x$Region_A, RegionB=x$Region_B)
    assign(n, graph_from_data_frame(links, nodes, directed = FALSE))
  }
  remove(x,i,s,links, nodes, n)
  # Compile Degree Centrality [integrated] (works) -----------------------------------------------
  #Country 1
  for (i in 1991:2017){
    n <- paste("net_1_", i, sep = "")
    x <- paste("dc_1_", i, sep="")
    assign(x, centr_degree(get(n), mode="total"))
  }
  #Country 2
  for (i in 1991:2017){
    n <- paste("net_2_", i, sep = "")
    x <- paste("dc_2_", i, sep="")
    assign(x, centr_degree(get(n), mode="total"))
  } 
  remove(i, n, x)
  #compile degree centrality  
  year <- c(1991:2017)
  #country 1: 
  cn <- c(dc_1_1991$centralization,
          dc_1_1992$centralization,
          dc_1_1993$centralization,
          dc_1_1994$centralization,
          dc_1_1995$centralization,
          dc_1_1996$centralization,
          dc_1_1997$centralization,
          dc_1_1998$centralization,
          dc_1_1999$centralization,
          dc_1_2000$centralization,
          dc_1_2001$centralization,
          dc_1_2002$centralization,
          dc_1_2003$centralization,
          dc_1_2004$centralization,
          dc_1_2005$centralization,
          dc_1_2006$centralization,
          dc_1_2007$centralization,
          dc_1_2008$centralization,
          dc_1_2009$centralization,
          dc_1_2010$centralization,
          dc_1_2011$centralization,
          dc_1_2012$centralization,
          dc_1_2013$centralization,
          dc_1_2014$centralization,
          dc_1_2015$centralization,
          dc_1_2016$centralization,
          dc_1_2017$centralization
  )
  dc_int1 <- data.frame(year, cn)
  remove(
    dc_1_1991,
    dc_1_1992,
    dc_1_1993,
    dc_1_1994,
    dc_1_1995,
    dc_1_1996,
    dc_1_1997,
    dc_1_1998,
    dc_1_1999,
    dc_1_2000,
    dc_1_2001,
    dc_1_2002,
    dc_1_2003,
    dc_1_2004,
    dc_1_2005,
    dc_1_2006,
    dc_1_2007,
    dc_1_2008,
    dc_1_2009,
    dc_1_2010,
    dc_1_2011,
    dc_1_2012,
    dc_1_2013,
    dc_1_2014,
    dc_1_2015,
    dc_1_2016,
    dc_1_2017)
  #country 2:  
  cn <- c(dc_2_1991$centralization,
          dc_2_1992$centralization,
          dc_2_1993$centralization,
          dc_2_1994$centralization,
          dc_2_1995$centralization,
          dc_2_1996$centralization,
          dc_2_1997$centralization,
          dc_2_1998$centralization,
          dc_2_1999$centralization,
          dc_2_2000$centralization,
          dc_2_2001$centralization,
          dc_2_2002$centralization,
          dc_2_2003$centralization,
          dc_2_2004$centralization,
          dc_2_2005$centralization,
          dc_2_2006$centralization,
          dc_2_2007$centralization,
          dc_2_2008$centralization,
          dc_2_2009$centralization,
          dc_2_2010$centralization,
          dc_2_2011$centralization,
          dc_2_2012$centralization,
          dc_2_2013$centralization,
          dc_2_2014$centralization,
          dc_2_2015$centralization,
          dc_2_2016$centralization,
          dc_2_2017$centralization
  )
  dc_int2 <- data.frame(year, cn)
  remove(year, cn)
  
  remove(
    dc_2_1991,
    dc_2_1992,
    dc_2_1993,
    dc_2_1994,
    dc_2_1995,
    dc_2_1996,
    dc_2_1997,
    dc_2_1998,
    dc_2_1999,
    dc_2_2000,
    dc_2_2001,
    dc_2_2002,
    dc_2_2003,
    dc_2_2004,
    dc_2_2005,
    dc_2_2006,
    dc_2_2007,
    dc_2_2008,
    dc_2_2009,
    dc_2_2010,
    dc_2_2011,
    dc_2_2012,
    dc_2_2013,
    dc_2_2014,
    dc_2_2015,
    dc_2_2016,
    dc_2_2017)

  # Plot Degree Centrality [integrated] (works) --------------------------------------------------------------------
  
  #plot  
  plot(dc_int2, type = "o",col = "red", xlab = "Year", ylab = "Degree Centrality", 
       main = "Degree centrality (integrated)", ylim = c(0.05,0.25)) 
  lines(dc_int1, type = "o", col = "blue") 
  legend(1991, 0.25, legend=c("C2:Internal + External", "C1:Internal + External"),
         col=c("red", "blue"), lty=1:2, cex=0.8)  
  
  
# Compute Harmonic Centrality [unweighted] (works) ---------------------------------------------
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
      par(old.par)
      plot(harm_cent22, type = "o",col = "red", xlab = "Year", ylab = "Harmonic Centrality (Mean)", 
           main = "Harmonic centrality (Mean)") 
      lines(harm_cent21, type = "o", col = "blue") 
      lines(harm_cent11, type = "o", col = "darkgreen")
      lines(harm_cent12, type = "o", col = "orange")
      legend(1991, 140, legend=c("C2:Internal", "C2:External", "C1:Internal", "C1:External"),
             col=c("red", "blue", "darkgreen", "orange"), lty=1:2, cex=0.8) 

      par(mfrow=c(2,3), oma = c(0, 0, 2, 0))
      hist(harm.cent_11_1991, xlim=c(0,120), xlab="", ylab="", main="1991", col="darkred")
      hist(harm.cent_11_1996, xlim=c(0,120), xlab="", ylab="", main="1996", col="darkred")
      hist(harm.cent_11_2001, xlim=c(0,120), xlab="", ylab="", main="2001", col="darkred")
      hist(harm.cent_11_2006, xlim=c(0,120), xlab="", ylab="", main="2006", col="darkred")
      hist(harm.cent_11_2011, xlim=c(0,120), xlab="", ylab="", main="2011", col="darkred")
      hist(harm.cent_11_2017, xlim=c(0,120), xlab="", ylab="", main="2017", col="darkred")
      mtext("Country 1: Internal links (Harmonic centrality)", 
            outer=TRUE, cex=1, font=2)
      
      par(mfrow=c(2,3), oma = c(0, 0, 2, 0))
      hist(harm.cent_12_1991, xlim=c(0,120), xlab="", ylab="", main="1991", col="darkred")
      hist(harm.cent_12_1996, xlim=c(0,120), xlab="", ylab="", main="1996", col="darkred")
      hist(harm.cent_12_2001, xlim=c(0,120), xlab="", ylab="", main="2001", col="darkred")
      hist(harm.cent_12_2006, xlim=c(0,120), xlab="", ylab="", main="2006", col="darkred")
      hist(harm.cent_12_2011, xlim=c(0,120), xlab="", ylab="", main="2011", col="darkred")
      hist(harm.cent_12_2017, xlim=c(0,120), xlab="", ylab="", main="2017", col="darkred")
      mtext("Country 1: External links (Harmonic centrality)", 
            outer=TRUE, cex=1, font=2)      

      par(mfrow=c(2,3), oma = c(0, 0, 2, 0))
      hist(harm.cent_22_1991, xlim=c(0,120), xlab="", ylab="", main="1991", col="darkred")
      hist(harm.cent_22_1996, xlim=c(0,120), xlab="", ylab="", main="1996", col="darkred")
      hist(harm.cent_22_2001, xlim=c(0,120), xlab="", ylab="", main="2001", col="darkred")
      hist(harm.cent_22_2006, xlim=c(0,120), xlab="", ylab="", main="2006", col="darkred")
      hist(harm.cent_22_2011, xlim=c(0,120), xlab="", ylab="", main="2011", col="darkred")
      hist(harm.cent_22_2017, xlim=c(0,120), xlab="", ylab="", main="2017", col="darkred")
      mtext("Country 2: Internal links (Harmonic centrality)", 
            outer=TRUE, cex=1, font=2)      
      
      par(mfrow=c(2,3), oma = c(0, 0, 2, 0))
      hist(harm.cent_21_1991, xlim=c(0,120), xlab="", ylab="", main="1991", col="darkred")
      hist(harm.cent_21_1996, xlim=c(0,120), xlab="", ylab="", main="1996", col="darkred")
      hist(harm.cent_21_2001, xlim=c(0,120), xlab="", ylab="", main="2001", col="darkred")
      hist(harm.cent_21_2006, xlim=c(0,120), xlab="", ylab="", main="2006", col="darkred")
      hist(harm.cent_21_2011, xlim=c(0,120), xlab="", ylab="", main="2011", col="darkred")
      hist(harm.cent_21_2017, xlim=c(0,120), xlab="", ylab="", main="2017", col="darkred")
      mtext("Country 2: External links (Harmonic centrality)", 
            outer=TRUE, cex=1, font=2)      
      
                    
# Compute Harmonic Centrality [weighted] (works) ---------------------------------------------
      #C1: internal
      harm_cent11 <- data.frame()
      for (i in 1991:2017){
        n <- paste("net_11_", i, sep="")
        x <- get(n)
        h <- paste("harm.cent_11_", i, sep="")
        assign(h, harmonic_centrality(x, mode = "all", weights = E(x)$w1))
        d <- data.frame(year=i, harmonic_centrality=mean(get(h)))
        harm_cent11 <- rbind(harm_cent11, d)
      }
      #C1: external
      harm_cent12 <- data.frame()
      for (i in 1991:2017){
        n <- paste("net_12_", i, sep="")
        x <- get(n)
        h <- paste("harm.cent_12_", i, sep="")
        assign(h, harmonic_centrality(x, mode = "all", weights = E(x)$w1))
        d <- data.frame(year=i, harmonic_centrality=mean(get(h)))
        harm_cent12 <- rbind(harm_cent12, d)
      }
      #C2: internal
      harm_cent22 <- data.frame()
      for (i in 1991:2017){
        n <- paste("net_22_", i, sep="")
        x <- get(n)
        h <- paste("harm.cent_22_", i, sep="")
        assign(h, harmonic_centrality(x, mode = "all", weights = E(x)$w1))
        d <- data.frame(year=i, harmonic_centrality=mean(get(h)))
        harm_cent22 <- rbind(harm_cent22, d)
      }
      #C2: external      
      harm_cent21 <- data.frame()
      for (i in 1991:2017){
        n <- paste("net_21_", i, sep="")
        x <- get(n)
        h <- paste("harm.cent_21_", i, sep="")
        assign(h, harmonic_centrality(x, mode = "all", weights = E(x)$w1))
        d <- data.frame(year=i, harmonic_centrality=mean(get(h)))
        harm_cent21 <- rbind(harm_cent21, d)
      } 
      harm_cent11$link <- "Country 1: Internal"
      harm_cent12$link <- "Country 1: External"      
      harm_cent22$link <- "Country 2: Internal"
      harm_cent21$link <- "Country 2: External"
      harm_cent <- rbind(harm_cent11, harm_cent12, harm_cent22, harm_cent21)
      
      p <- ggplot(harm_cent, aes(year, harmonic_centrality))
      p +  ggtitle("Harmonic Centrality [weigted]") + 
        theme(plot.title = element_text(hjust = 0.5)) +
        geom_point(aes(colour = factor(link)))
      
      par(old.par)
      plot(harm_cent22, type = "o",col = "red", xlab = "Year", ylab = "Harmonic Centrality (Mean)", 
           main = "Harmonic centrality (Mean)", ylim=c(0,30)) 
      lines(harm_cent21, type = "o", col = "blue") 
      lines(harm_cent11, type = "o", col = "darkgreen")
      lines(harm_cent12, type = "o", col = "orange")
      legend(1991, 30, legend=c("C2:Internal", "C2:External", "C1:Internal", "C1:External"),
             col=c("red", "blue", "darkgreen", "orange"), lty=1:2, cex=0.8) 

#density plots             
      par(mfrow=c(2,2), oma=c(0,0,2,0))
      plot(density(harm.cent_11_1991), xlab="", ylab="density", main="Country 1: Internal links", cex=0.8, 
           col="red", xlim = c(0,60),ylim=c(0,0.4)) 
      polygon(density(harm.cent_11_1991), col=adjustcolor("red", alpha.f = 0.2), border="black")
      lines(density(harm.cent_11_2001), col="blue")
      polygon(density(harm.cent_11_2001), col=adjustcolor("blue", alpha.f = 0.2), border="black")
      lines(density(harm.cent_11_2017), col="green")
      polygon(density(harm.cent_11_2017), col=adjustcolor("green", alpha.f = 0.2), border="black")
      
      plot(density(harm.cent_12_1991), xlab="", ylab="density", main="Country 1: External links", cex=0.8, 
           col="red", xlim = c(0,60),ylim=c(0,0.4)) 
      polygon(density(harm.cent_12_1991), col=adjustcolor("red", alpha.f = 0.2), border="black")
      lines(density(harm.cent_12_2001), col="blue")
      polygon(density(harm.cent_12_2001), col=adjustcolor("blue", alpha.f = 0.2), border="black")
      lines(density(harm.cent_12_2017), col="green")
      polygon(density(harm.cent_12_2017), col=adjustcolor("green", alpha.f = 0.2), border="black")
      
      plot(density(harm.cent_22_1991), xlab="", ylab="density", main="Country 2: Internal links", cex=0.8, 
           col="red", xlim = c(0,60),ylim=c(0,0.4)) 
      polygon(density(harm.cent_22_1991), col=adjustcolor("red", alpha.f = 0.2), border="black")
      lines(density(harm.cent_22_2001), col="blue")
      polygon(density(harm.cent_22_2001), col=adjustcolor("blue", alpha.f = 0.2), border="black")
      lines(density(harm.cent_22_2017), col="green")
      polygon(density(harm.cent_22_2017), col=adjustcolor("green", alpha.f = 0.2), border="black")
      
      plot(density(harm.cent_21_1991), xlab="", ylab="density", main="Country 2: External links", cex=0.8, 
           col="red", xlim = c(0,60),ylim=c(0,0.4)) 
      polygon(density(harm.cent_21_1991), col=adjustcolor("red", alpha.f = 0.2), border="black")
      lines(density(harm.cent_21_2001), col="blue")
      polygon(density(harm.cent_21_2001), col=adjustcolor("blue", alpha.f = 0.2), border="black")
      lines(density(harm.cent_21_2017), col="green")
      polygon(density(harm.cent_21_2017), col=adjustcolor("green", alpha.f = 0.2), border="black")
   
      mtext("Harmonic Centrality [Weighted]", 
            outer=TRUE, cex=1, font=2)

#number plot 
      #country 1
      t1a <- sort(harm.cent_11_1991)
      t1a <- data.frame(num=1:189, harmonic_centrality=t1a)
      t1b <- sort(harm.cent_12_1991)
      t1b <- data.frame(num=1:189, harmonic_centrality=t1b)      
      t2a <- sort(harm.cent_11_1996)
      t2a <- data.frame(num=1:189, harmonic_centrality=t2a)
      t2b <- sort(harm.cent_12_1996)
      t2b <- data.frame(num=1:189, harmonic_centrality=t2b)      
      t3a <- sort(harm.cent_11_2001)
      t3a <- data.frame(num=1:189, harmonic_centrality=t3a)
      t3b <- sort(harm.cent_12_2001)
      t3b <- data.frame(num=1:189, harmonic_centrality=t3b)      
      t4a <- sort(harm.cent_11_2006)
      t4a <- data.frame(num=1:189, harmonic_centrality=t4a)
      t4b <- sort(harm.cent_12_2006)
      t4b <- data.frame(num=1:189, harmonic_centrality=t4b)      
      t5a <- sort(harm.cent_11_2011)
      t5a <- data.frame(num=1:189, harmonic_centrality=t5a)
      t5b <- sort(harm.cent_12_2011)
      t5b <- data.frame(num=1:189, harmonic_centrality=t5b)      
      t6a <- sort(harm.cent_11_2017)
      t6a <- data.frame(num=1:189, harmonic_centrality=t6a)
      t6b <- sort(harm.cent_12_2017)
      t6b <- data.frame(num=1:189, harmonic_centrality=t6b)      
      
      par(mfrow=c(2,3), oma = c(0,2,2,0))
      plot(t1a, ylim=c(0,210), xlab="1991", ylab="Harmonic Centrality", type="h", col=adjustcolor("blue", alpha.f = 0.2))
      lines(t1b, ylim=c(0,210), xlab="", ylab="", type="h", col=adjustcolor("red", alpha.f = 0.2)) 
      legend(1,200,legend=c("Internal", "External"),col=c("blue", "red"), lty=1, cex=0.9, bty="n") 
      
      plot(t2a, ylim=c(0,210), xlab="1996", ylab="", type="h", col=adjustcolor("blue", alpha.f = 0.2))
      lines(t2b, ylim=c(0,210), xlab="", ylab="", type="h", col=adjustcolor("red", alpha.f = 0.2)) 
      
      plot(t3a, ylim=c(0,210), xlab="2001", ylab="", type="h", col=adjustcolor("blue", alpha.f = 0.2))
      lines(t3b, ylim=c(0,210), xlab="", ylab="", type="h", col=adjustcolor("red", alpha.f = 0.2)) 
      
      plot(t4a, ylim=c(0,210), xlab="2006", ylab="Harmonic Centrality", type="h", col=adjustcolor("blue", alpha.f = 0.2))
      lines(t4b, ylim=c(0,210), xlab="", ylab="", type="h", col=adjustcolor("red", alpha.f = 0.2)) 
      
      plot(t5a, ylim=c(0,210), xlab="2011", ylab="", type="h", col=adjustcolor("blue", alpha.f = 0.2))
      lines(t5b, ylim=c(0,210), xlab="", ylab="", type="h", col=adjustcolor("red", alpha.f = 0.2)) 
      
      plot(t6a, ylim=c(0,210), xlab="2017", ylab="", type="h", col=adjustcolor("blue", alpha.f = 0.2))
      lines(t6b, ylim=c(0,210), xlab="", ylab="", type="h", col=adjustcolor("red", alpha.f = 0.2)) 
      mtext("Country 1: Harmonic Centrality [weighted] per vertex", 
            outer=TRUE, cex=0.9, font=2)
      
      #country 2
      t1a <- sort(harm.cent_21_1991)
      t1a <- data.frame(num=1:189, harmonic_centrality=t1a)
      t1b <- sort(harm.cent_22_1991)
      t1b <- data.frame(num=1:189, harmonic_centrality=t1b)      
      t2a <- sort(harm.cent_21_1996)
      t2a <- data.frame(num=1:189, harmonic_centrality=t2a)
      t2b <- sort(harm.cent_22_1996)
      t2b <- data.frame(num=1:189, harmonic_centrality=t2b)      
      t3a <- sort(harm.cent_21_2001)
      t3a <- data.frame(num=1:189, harmonic_centrality=t3a)
      t3b <- sort(harm.cent_22_2001)
      t3b <- data.frame(num=1:189, harmonic_centrality=t3b)      
      t4a <- sort(harm.cent_21_2006)
      t4a <- data.frame(num=1:189, harmonic_centrality=t4a)
      t4b <- sort(harm.cent_22_2006)
      t4b <- data.frame(num=1:189, harmonic_centrality=t4b)      
      t5a <- sort(harm.cent_21_2011)
      t5a <- data.frame(num=1:189, harmonic_centrality=t5a)
      t5b <- sort(harm.cent_22_2011)
      t5b <- data.frame(num=1:189, harmonic_centrality=t5b)      
      t6a <- sort(harm.cent_21_2017)
      t6a <- data.frame(num=1:189, harmonic_centrality=t6a)
      t6b <- sort(harm.cent_22_2017)
      t6b <- data.frame(num=1:189, harmonic_centrality=t6b)      
      
      par(mfrow=c(2,3), oma = c(0, 2, 2, 0))
      plot(t1a, ylim=c(0,100), xlab="1991", ylab="Harmonic Centrality", type="h", col=adjustcolor("blue", alpha.f = 0.2))
      lines(t1b, ylim=c(0,100), xlab="", ylab="", type="h", col=adjustcolor("red", alpha.f = 0.2)) 
      legend(1,100,legend=c("Internal", "External"),col=c("blue", "red"), lty=1, cex=0.9, bty="n") 
      
      plot(t2a, ylim=c(0,100), xlab="1996", ylab="", type="h", col=adjustcolor("blue", alpha.f = 0.2))
      lines(t2b, ylim=c(0,100), xlab="", ylab="", type="h", col=adjustcolor("red", alpha.f = 0.2)) 
      
      plot(t3a, ylim=c(0,100), xlab="2001", ylab="", type="h", col=adjustcolor("blue", alpha.f = 0.2))
      lines(t3b, ylim=c(0,100), xlab="", ylab="", type="h", col=adjustcolor("red", alpha.f = 0.2)) 
      
      plot(t4a, ylim=c(0,100), xlab="2006", ylab="Harmonic Centrality", type="h", col=adjustcolor("blue", alpha.f = 0.2))
      lines(t4b, ylim=c(0,100), xlab="", ylab="", type="h", col=adjustcolor("red", alpha.f = 0.2)) 
      
      plot(t5a, ylim=c(0,100), xlab="2011", ylab="", type="h", col=adjustcolor("blue", alpha.f = 0.2))
      lines(t5b, ylim=c(0,100), xlab="", ylab="", type="h", col=adjustcolor("red", alpha.f = 0.2)) 
      
      plot(t6a, ylim=c(0,100), xlab="2017", ylab="", type="h", col=adjustcolor("blue", alpha.f = 0.2))
      lines(t6b, ylim=c(0,100), xlab="", ylab="", type="h", col=adjustcolor("red", alpha.f = 0.2)) 
      mtext("Country 2: Harmonic Centrality [weighted] per vertex", 
            outer=TRUE, cex=0.9, font=2)
      
      

#histogram      
      par(mfrow=c(2,3), oma = c(0, 0, 2, 0))
      hist(harm.cent_11_1991, xlab="", ylab="", main="1991", col="darkred")
      hist(harm.cent_11_1996, xlab="", ylab="", main="1996", col="darkred")
      hist(harm.cent_11_2001, xlab="", ylab="", main="2001", col="darkred")
      hist(harm.cent_11_2006, xlab="", ylab="", main="2006", col="darkred")
      hist(harm.cent_11_2011, xlab="", ylab="", main="2011", col="darkred")
      hist(harm.cent_11_2017, xlab="", ylab="", main="2017", col="darkred")
      mtext("Country 1: Internal links (Harmonic centrality)", 
            outer=TRUE, cex=1, font=2)      
      
            
      par(mfrow=c(2,3), oma = c(0, 0, 2, 0))
      hist(harm.cent_12_1991, xlab="", ylab="", main="1991", col="darkred")
      hist(harm.cent_12_1996, xlab="", ylab="", main="1996", col="darkred")
      hist(harm.cent_12_2001, xlab="", ylab="", main="2001", col="darkred")
      hist(harm.cent_12_2006, xlab="", ylab="", main="2006", col="darkred")
      hist(harm.cent_12_2011, xlab="", ylab="", main="2011", col="darkred")
      hist(harm.cent_12_2017, xlab="", ylab="", main="2017", col="darkred")
      mtext("Country 1: External links (Harmonic centrality)", 
            outer=TRUE, cex=1, font=2)      
      
      par(mfrow=c(2,3), oma = c(0, 0, 2, 0))
      hist(harm.cent_22_1991, xlab="", ylab="", main="1991", col="darkred")
      hist(harm.cent_22_1996, xlab="", ylab="", main="1996", col="darkred")
      hist(harm.cent_22_2001, xlab="", ylab="", main="2001", col="darkred")
      hist(harm.cent_22_2006, xlab="", ylab="", main="2006", col="darkred")
      hist(harm.cent_22_2011, xlab="", ylab="", main="2011", col="darkred")
      hist(harm.cent_22_2017, xlab="", ylab="", main="2017", col="darkred")
      mtext("Country 2: Internal links (Harmonic centrality)", 
            outer=TRUE, cex=1, font=2)      
      
      par(mfrow=c(2,3), oma = c(0, 0, 2, 0))
      hist(harm.cent_21_1991, xlab="", ylab="", main="1991", col="darkred")
      hist(harm.cent_21_1996, xlab="", ylab="", main="1996", col="darkred")
      hist(harm.cent_21_2001, xlab="", ylab="", main="2001", col="darkred")
      hist(harm.cent_21_2006, xlab="", ylab="", main="2006", col="darkred")
      hist(harm.cent_21_2011, xlab="", ylab="", main="2011", col="darkred")
      hist(harm.cent_21_2017, xlab="", ylab="", main="2017", col="darkred")
      mtext("Country 2: External links (Harmonic centrality)", 
            outer=TRUE, cex=1, font=2)      
      
      
# Compute closeness centrality (not recommended)--------------------------------------
  nodes <- data.frame(ID=1:189)
  
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
  
# Newman-Garvin clustering (not very usefull)------------------------------------------------
par(old.par)
  ceb <- cluster_edge_betweenness(net_11_1991)
  #dendPlot(ceb, mode="hclust")
  #plot(ceb, net_11_1991)
  modularity(ceb)
  #membership(ceb)

  ceb2 <- cluster_edge_betweenness(net_11_2017)
  #dendPlot(ceb2, mode="hclust")
  #plot(ceb, net_11_1991)
  modularity(ceb2)
  #membership(ceb2)
  
  ceb3 <- cluster_edge_betweenness(net_12_1991)
  #dendPlot(ceb, mode="hclust")
  #plot(ceb, net_11_1991)
  modularity(ceb)
  #membership(ceb)
  
  ceb4 <- cluster_edge_betweenness(net_12_2017)
  #dendPlot(ceb2, mode="hclust")
  #plot(ceb, net_11_1991)
  modularity(ceb2)
  #membership(ceb2)
  
  ceb5 <- cluster_edge_betweenness(net_22_1991)
  #dendPlot(ceb, mode="hclust")
  #plot(ceb, net_11_1991)
  modularity(ceb)
  #membership(ceb)
  
  ceb6 <- cluster_edge_betweenness(net_22_2017)
  #dendPlot(ceb2, mode="hclust")
  #plot(ceb, net_11_1991)
  modularity(ceb2)
  #membership(ceb2)
  
  ceb3 <- cluster_edge_betweenness(net_12_1991)
  #dendPlot(ceb, mode="hclust")
  #plot(ceb, net_11_1991)
  modularity(ceb)
  #membership(ceb)
  
  ceb4 <- cluster_edge_betweenness(net_12_2017)
  #dendPlot(ceb2, mode="hclust")
  #plot(ceb, net_11_1991)
  modularity(ceb2)
  #membership(ceb2)
  
  
    
  
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

        