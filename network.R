####################################################
#NETWORK DATA ANALYSIS
####################################################

#### INITIALIZE-------------------------------------------------------------------------------------
  # Install / Load packages --------------------------------------------------------------
#Load block, if these packages are not installed, you have to remove the pound sign and run each code. 
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
#install.packages("plm")
#install.packages("lmtest")

#Then load the packages into the current working session:
library(igraph)
library(dplyr)
library(sqldf)
library(readstata13)
library(CINNA)
library(ggplot2)
#library(sna)
library(qgraph)
library(data.table)
library(plm)
library(lmtest)
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


#### DATA MANAGEMENT---------------------------------------------------------------------------------
  # Create partial data set -------------------------------------------------
  par_set <- `Network Data`

#normalize Intensitiy measure (range between 0-1)
  par_set$Intensity_norm = (
    par_set$Intensity-min(par_set$Intensity)) / 
    (max(par_set$Intensity)-min(par_set$Intensity)
    )
  cor(par_set$Intensity, par_set$Intensity_norm) == 1
  
#Order regions so that they are consistant across years
  t1 <- par_set %>%
    filter(Region_A <= Region_B) %>%
    select(
      Year, 
      Region_A, 
      Country_A, 
      Region_B, 
      Country_B, 
      Intensity, 
      Intensity_norm)
  t2 <- par_set %>%
    filter(Region_A > Region_B) %>%
    select(
      Year, 
      Region_B=Region_A, 
      Country_B=Country_A, 
      Region_A=Region_B, 
      Country_A=Country_B, 
      Intensity, 
      Intensity_norm)
  par_set <- rbind(t1, t2)
  remove(t1,t2)
  
#Add dummies for auto-link (region connected to itself) 
  par_set <- par_set %>%
    mutate(autolink=ifelse(Country_A==Country_B & Region_A==Region_B, 1, 0)) 

#Add dummies for links across borders    
  par_set <- par_set %>%
    mutate(border=ifelse(Country_A!=Country_B, 1, 0))

#Add dummies for countries 1&2 per link (and not per region)
par_set <- par_set %>% 
  mutate(c1=ifelse(Country_A==1|Country_B==1, 1,0)) %>%
  mutate(c2=ifelse(Country_A==2|Country_B==2, 1,0))

#Add ID for link 
  par_set$edge1 <- paste(par_set$Region_A, "-", par_set$Region_B, sep="")
  
  par_set$edge2 <- paste(par_set$Region_B, "-", par_set$Region_A, sep="")
  
#Create unique row id
  par_set$rowid <- 
    paste(
      par_set$Year,"-",
      par_set$edge1,"-", sep="")
  
#Test for duplicate rows
  test <- sqldf("select rowid, count(rowid) 
                from par_set group by rowid 
                having count(rowid)>1")
  paste("Data has no duplicate rows:",nrow(test)==0)
  remove(test)  
  
  # Create full data set ---------------------------------------------------------

  #IDs 
  #1 thru 22 are Country 2
  # 23 thru 99 are Country 1 
  # 100 thru 140 are Country 2  
  # 141 thru 189 are Country 1
  regions_c1 <- c(23:99, 141:189)
  regions_c2 <- c(1:22, 100:140) 
  regions_c1c2 <- c(regions_c1,regions_c2)
  
#Test number of autolinks per year: 
  autolink_count <- data.frame()
  for (i in 1991:2017){
    t <- paste("select Year, Country_A, count(Country_A) from par_set where autolink = 1 and Year=",i,
               " group by Country_A, Year",sep="")
    test <- sqldf(t)
    autolink_count <- rbind(autolink_count, test)
  } 
  remove(test,t,i)
  
#Which autolinks are missing:  
  autolink_check <- data.frame()
  for (y in 1991:2017){
    t <- paste("select * from par_set where autolink = 1 and Year=",y,sep="")
    test <- sqldf(t)
    for (r in 1:189){
      if (r %in% test$Region_A==F){
        s1 <- data.frame(Year=y, Region_miss=r)
        autolink_check <- rbind(autolink_check, s1)
      }
    }
  }
  remove(y,t,test,r,s1)
  
#Autolink for region 176 is missing for year 1991-1996
  #add as links with intensitiy 0 
  autolink_add <- data.frame(Year=1991:1996, Region_A=176, c1=1, 
                             Region_B=176, c2=0, Intensity=0)
  remove(autolink_check, autolink_count)
  
#Country 1: internal region connection list
  regions_int_c1 <- data.frame()
  for (a in regions_c1){
    for (i in regions_c1){
      c <- expand.grid(i,a)
        if(c$Var1 < c$Var2){
          regions_int_c1 <- rbind(regions_int_c1, c)        
        }
      }
    }
  remove(a,i,c)

#Country 2: internal region connection list
  regions_int_c2 <- data.frame()
  for (a in regions_c2){
    for (i in regions_c2){
      c <- expand.grid(i,a)
      if(c$Var1 < c$Var2){
        regions_int_c2 <- rbind(regions_int_c2, c)        
      }
    }
  }
  remove(a,i,c)
  
#Country 1-2: external region list for cross-border conncetions
  regions_ext <- data.frame()
  for (a in regions_c2){
    for (i in regions_c1){
      c <- expand.grid(a,i)
          regions_ext <- rbind(regions_ext, c)        
        }
      }
  remove(a,i,c)
  #Reorder external regions list
    t1 <- regions_ext %>%
      filter(Var1 <= Var2) %>%
      select(Var1, Var2)
    t2 <- regions_ext %>%
      filter(Var1 > Var2) %>%
      select(Var2=Var1, Var1=Var2)
    regions_ext <- rbind(t1, t2)
    remove(t1,t2)
  
  
#test for duplicates in dummy lists (Should return TRUE)
  regions_int_c1$edge1 <- paste(regions_int_c1$Var1, "-", regions_int_c1$Var2, sep="")
  regions_int_c1$edge2 <- paste(regions_int_c1$Var2, "-", regions_int_c1$Var1, sep="")
  test1 <- sqldf("select edge1 from regions_int_c1 except select edge2 from regions_int_c1")
    regions_int_c2$edge1 <- paste(regions_int_c2$Var1, "-", regions_int_c2$Var2, sep="")
    regions_int_c2$edge2 <- paste(regions_int_c2$Var2, "-", regions_int_c2$Var1, sep="")
    test2 <- sqldf("select edge1 from regions_int_c2 except select edge2 from regions_int_c2")
      regions_ext$edge1 <- paste(regions_ext$Var1, "-", regions_ext$Var2, sep="")
      regions_ext$edge2 <- paste(regions_ext$Var2, "-", regions_ext$Var1, sep="")
      test3 <- sqldf("select edge1 from regions_ext except select edge2 from regions_ext")
        paste("Internal region list c1 has no duplicate rows:", nrow(regions_int_c1)==nrow(test1))
        paste("Internal region list c2 has no duplicate rows:", nrow(regions_int_c2)==nrow(test2))
        paste("External region list has no duplicate rows:", nrow(regions_ext)==nrow(test3))
        remove(test1,test2,test3)

#add rows to par_set
  full_set <- select (par_set,c(Year, Region_A, c1, Region_B, c2,Intensity))
#add autolinks 
  full_set <- rbind(full_set, autolink_add)
  remove(autolink_add)
#country 1: internal
  for (y in 1991:2017) {
    t <- paste("select edge1 from regions_int_c1  
                except 
                  select edge1 from par_set where Year=",y," and c1=1 and c2=0 
                    union 
                    select edge2 as edge1 from par_set 
                      where Year=",y," and c1=1 and c2=0",sep="")
    test <- sqldf(t)
    a <- paste("select ",y," as Year, Var1 as Region_A, 1 as c1, 
                  Var2 as Region_B, 0 as c2,0 as Intensity
                  from regions_int_c1 where edge1 in 
                  (select edge1 from test)")
    addrow <- sqldf(a)
    full_set <- rbind(full_set, addrow)
  }
  #country 2: internal
    for (y in 1991:2017) {
      t <- paste("select edge1 from regions_int_c2  
                  except 
                    select edge1 from par_set where Year=",y," and c1=0 and c2=1 
                      union 
                      select edge2 as edge1 from par_set
                        where Year=",y," and c1=0 and c2=1",sep="")
      test <- sqldf(t)
      a <- paste("select ",y," as Year, Var1 as Region_A, 0 as c1, 
                    Var2 as Region_B, 1 as c2,0 as Intensity
                    from regions_int_c2 where edge1 in 
                    (select edge1 from test)")
      addrow <- sqldf(a)
      full_set <- rbind(full_set, addrow)
    }
    #country 1-2 external
      for (y in 1991:2017) {
        t <- paste("select edge1 from regions_ext  
                    except 
                      select edge1 from par_set where Year=",y," and c1=1 and c2=1 
                        union 
                        select edge2 as edge1 from par_set
                          where Year=",y," and c1=1 and c2=1",sep="")
        test <- sqldf(t)
        a <- paste("select ",y," as Year, Var1 as Region_A, 1 as c1, 
                      Var2 as Region_B, 1 as c2,0 as Intensity
                      from regions_ext where edge1 in 
                      (select edge1 from test)")
        addrow <- sqldf(a)
        full_set <- rbind(full_set, addrow)
      }
      remove(t,y,a,test,addrow, regions_int_c1, regions_int_c2,regions_ext)

#add dummies for autolink (region connected to itself) and links across borders
  full_set <- full_set %>%
    mutate(autolink=ifelse(c1+c2<=1 & Region_A==Region_B, 1, 0)) %>%
    mutate(border=ifelse(c1+c2>=2, 1, 0))

#create new unique region ID (rid1 / rid2) 
  full_set$edge1 <- paste(full_set$Region_A, "-", full_set$Region_B, sep="")
  full_set$edge2 <- paste(full_set$Region_B, "-", full_set$Region_A, sep="")

#create unique row id
  full_set$rowid <- 
    paste(full_set$Year,"-",full_set$edge1,"-",full_set$c1,"-",full_set$c2,"-",full_set$autolink, sep="")

#test duplicates 
  test <- sqldf("select rowid, count(rowid) from full_set group by rowid having count(rowid)>1")
  paste("Data has no duplicate rows:",nrow(test)==0)
  remove(test)

#normalize Intensitiy measure (range between 0-1)
  full_set$Intensity_norm = (
    full_set$Intensity-min(full_set$Intensity)) / 
    (max(full_set$Intensity)-min(full_set$Intensity)
    )
  cor(full_set$Intensity, full_set$Intensity_norm) == 1

  # Construct adjacency matrix per year, country, links [undirected] ------------------------------------------------
  nodes_c1 <- data.frame(ID=regions_c1)
  nodes_c2 <- data.frame(ID=regions_c2)
  nodes_c1c2 <- data.frame(ID=regions_c1c2)
  
  #Country 1: Internal links 
  for (i in 1991:2017) {
    s <- paste("select * from par_set where Year = ", i, "and c1=1 and c2=0")
    x <- sqldf(s)
    n <- paste("net_11_", i, sep="")
    links <- data.frame(RegionA=x$Region_A, RegionB=x$Region_B, w1=x$Intensity, w2=x$Intensity_norm)
    assign(n, graph_from_data_frame(links, nodes_c1, directed = F))
  }
  #Country 1-2: External links  
  for (i in 1991:2017) {
    s <- paste("select * from par_set where Year = ", i, " and c1=1 and c2=1")
    x <- sqldf(s)
    n <- paste("net_12_", i, sep="")
    links <- data.frame(RegionA=x$Region_A, RegionB=x$Region_B, w1=x$Intensity, w2=x$Intensity_norm)
    assign(n, graph_from_data_frame(links, nodes_c1c2, directed = F))
  }
  #Country 2: Internal links  
  for (i in 1991:2017) {
    s <- paste("select * from par_set where Year = ", i, " and c1=0 and c2=1")
    x <- sqldf(s)
    n <- paste("net_22_", i, sep="")
    links <- data.frame(RegionA=x$Region_A, RegionB=x$Region_B, w1=x$Intensity, w2=x$Intensity_norm)
    assign(n, graph_from_data_frame(links, nodes_c2, directed = F))
  }
  remove(x,i,s,links, nodes_c1, nodes_c2, nodes_c1c2, n)
  
  # #example plot with no loops or unconnected regions
  # par(old.par)
  # plot(delete.vertices(simplify(net_11_2017), degree(net_11_2017)==0), 
  #      layout=layout_with_fr, vertex.label = NA, vertex.size=6)  
  
  
  
  

#### DATA ANALYSIS -----------------------------------------------------------------------------------
  # Partial set link intensity plot -----------------------------------------
#This code block will plot the link intesity per country (internal/external) per year using partial data
  pc1c1 <- par_set %>%
    filter(c1==1, c2==0) %>%
    group_by(Year) %>%
    summarise(mean(Intensity))
  pc1c2 <- par_set %>%
    filter(c1==1, c2==1) %>%
    group_by(Year) %>%
    summarise(mean(Intensity))
  pc2c2 <- par_set %>%
    filter(c1==0, c2==1) %>%
    group_by(Year) %>%
    summarise(mean(Intensity))
#plot
par(mfrow=c(1,2), oma = c(0, 0, 0, 0))
plot(pc2c2,type = "l",col = "red", xlab = "Year", ylab = "Intensity", 
       main = "Network Intensity (partial)", ylim=c(0,60000), lwd=2) 
  lines(pc1c1, type = "l", col = "darkgreen", lwd=2)
  lines(pc1c2, type = "l", col = "blue", lwd=2)
  legend(1991, 60000, legend=c("C2:Internal", "C1:Internal", "C1-C2:External"),
         col=c("red", "darkgreen", "blue"), lty=1, cex=0.8, bty="n", lwd=2)
  remove(pc1c1, pc1c2, pc2c2)

  # Full set link intensity plot --------------------------------------------
#This code block will plot the link intesity per country (internal/external) per year using the full set
  pc1c1 <- par_set %>%
    filter(c1==1, c2==0) %>%
    group_by(Year) %>%
    summarise(mean(Intensity))
  pc1c2 <- par_set %>%
    filter(c1==1, c2==1) %>%
    group_by(Year) %>%
    summarise(mean(Intensity))
  pc2c2 <- par_set %>%
    filter(c1==0, c2==1) %>%
    group_by(Year) %>%
    summarise(mean(Intensity))
  fc1c1 <- full_set %>%
    filter(c1==1, c2==0) %>%
    group_by(Year) %>%
    summarise(mean(Intensity))
  fc1c2 <- full_set %>%
    filter(c1==1, c2==1) %>%
    group_by(Year) %>%
    summarise(mean(Intensity))
  fc2c2 <- full_set %>%
    filter(c1==0, c2==1) %>%
    group_by(Year) %>%
    summarise(mean(Intensity))
  #plot
    par(mfrow=c(1,2), oma = c(0, 0, 0, 0))
    plot(pc2c2,type = "l",col = "red", xlab = "Year", ylab = "Intensity", 
         main = "Network Intensity (partial)", ylim=c(0,60000), lwd=2) 
    lines(pc1c1, type = "l", col = "darkgreen", lwd=2)
    lines(pc1c2, type = "l", col = "blue", lwd=2)
    legend(1991, 60000, legend=c("C2:Internal", "C1:Internal", "C1-C2:External"),
           col=c("red", "darkgreen", "blue"), lty=1, cex=0.8, bty="n", lwd=2)
        plot(fc1c1,type = "l",col = "darkgreen", xlab = "Year", ylab = "", 
             main = "Network Intensity (full)", ylim=c(0, 60000), lty=2, lwd=2)
        lines(pc1c1, type = "l", col = "darkgreen", lty=1, lwd=2)
        lines(fc1c2, type = "l", col = "blue", lty=2, lwd=2)
        lines(pc1c2, type = "l", col = "blue", lty=1, lwd=2)
        lines(fc2c2, type = "l", col = "red", lty=2, lwd=2)
        lines(pc2c2, type = "l", col = "red", lty=1, lwd=2)
        legend(1991, 60000, legend=c("C2:Internal", "C1:Internal", "C1-C2:External"),
               col=c("red", "darkgreen", "blue"), lty=2, cex=0.8, bty="n", lwd=2)
        remove(fc1c1, fc2c2, fc1c2, pc1c1, pc1c2, pc2c2)
        
        
  
  # Count zeros --------------------------------------------
  zc1c1 <- full_set %>%
    filter(c1==1, c2==0 & Intensity==0) %>%
    group_by(Year) %>%
    count(Intensity)
  zc1c1$share <- round(zc1c1$n/189)/189
  zc1c2 <- full_set %>%
    filter(c1==1, c2==1 & Intensity==0) %>%
    group_by(Year) %>%
    count(Intensity)
  zc1c2$share <- round(zc1c2$n/189)/189
  zc2c2 <- full_set %>%
    filter(c1==0, c2==1 & Intensity==0) %>%
    group_by(Year) %>%
    count(Intensity)
  zc2c2$share <- round(zc2c2$n/189)/189
  zc1c1$link <- "Country 1: Internal"
  zc1c2$link <- "Country 1-2: External"
  zc2c2$link <- "Country 2: Internal"
  zero.vertices <- rbind(zc1c1, zc1c2, zc2c2)  
  remove(zc1c1, zc1c2, zc2c2)
  
   p <- ggplot(zero.vertices, aes(Year, share))
   p +  ggtitle("Share of regions with no connections") + 
     theme(plot.title = element_text(hjust = 0.5)) +
     geom_smooth(aes(colour = factor(link))) +
     coord_cartesian(xlim = c(1991, 2017), ylim = c(0, 0.25))
   remove(p, zero.vertices)
  
  # Degree distibution plots [not normalized] [unweighted] ------------------------------------------------
#C1-2 external links 
  deg12_1 <- degree(net_12_1991)
  deg12_2 <- degree(net_12_1996)
  deg12_3 <- degree(net_12_2001)
  deg12_4 <- degree(net_12_2006)
  deg12_5 <- degree(net_12_2011)
  deg12_6 <- degree(net_12_2017)  
  deg.dist12_1 <- degree_distribution(net_12_1991, cumulative=T)
  deg.dist12_2 <- degree_distribution(net_12_1996, cumulative=T)
  deg.dist12_3 <- degree_distribution(net_12_2001, cumulative=T)
  deg.dist12_4 <- degree_distribution(net_12_2006, cumulative=T)
  deg.dist12_5 <- degree_distribution(net_12_2011, cumulative=T)
  deg.dist12_6 <- degree_distribution(net_12_2017, cumulative=T)  
#C1 internal links  
  deg11_1 <- degree(net_11_1991)
  deg11_2 <- degree(net_11_1996)
  deg11_3 <- degree(net_11_2001)
  deg11_4 <- degree(net_11_2006)
  deg11_5 <- degree(net_11_2011)
  deg11_6 <- degree(net_11_2017)
  deg.dist11_1 <- degree_distribution(net_11_1991, cumulative=T)
  deg.dist11_2 <- degree_distribution(net_11_1996, cumulative=T)
  deg.dist11_3 <- degree_distribution(net_11_2001, cumulative=T)
  deg.dist11_4 <- degree_distribution(net_11_2006, cumulative=T)
  deg.dist11_5 <- degree_distribution(net_11_2011, cumulative=T)
  deg.dist11_6 <- degree_distribution(net_11_2017, cumulative=T)
#C2 internal links  
  deg22_1 <- degree(net_22_1991)
  deg22_2 <- degree(net_22_1996)
  deg22_3 <- degree(net_22_2001)
  deg22_4 <- degree(net_22_2006)
  deg22_5 <- degree(net_22_2011)
  deg22_6 <- degree(net_22_2017)
  deg.dist22_1 <- degree_distribution(net_22_1991, cumulative=T)
  deg.dist22_2 <- degree_distribution(net_22_1996, cumulative=T)
  deg.dist22_3 <- degree_distribution(net_22_2001, cumulative=T)
  deg.dist22_4 <- degree_distribution(net_22_2006, cumulative=T)
  deg.dist22_5 <- degree_distribution(net_22_2011, cumulative=T)
  deg.dist22_6 <- degree_distribution(net_22_2017, cumulative=T)

#plot c1  
  par(mfrow=c(2,3), oma = c(0,2,2,0))
  plot(x=0:max(deg11_1), y=1-deg.dist11_1, cex=1.2, col=rgb(0,0.4,0.8,alpha=1), ylim=c(0,1), xlim=c(0,127), 
        xlab="1991", ylab="Cumulative Frequency", type = "h")
  plot(x=0:max(deg11_2), y=1-deg.dist11_2, cex=1.2, col=rgb(0,0.4,0.8,alpha=1), ylim=c(0,1), xlim=c(0,127), 
        xlab="1996", ylab="", type = "h")
  plot(x=0:max(deg11_3), y=1-deg.dist11_3, cex=1.2, col=rgb(0,0.4,0.8,alpha=1), ylim=c(0,1), xlim=c(0,127), 
        xlab="2001", ylab="", type = "h")
  plot( x=0:max(deg11_4), y=1-deg.dist11_4, cex=1.2, col=rgb(0,0.4,0.8,alpha=1), ylim=c(0,1), xlim=c(0,127), 
        xlab="2006", ylab="Cumulative Frequency", type = "h")
  plot( x=0:max(deg11_5), y=1-deg.dist11_5, cex=1.2, col=rgb(0,0.4,0.8,alpha=1), ylim=c(0,1), xlim=c(0,127), 
        xlab="2011", ylab="", type = "h")
  plot(x=0:max(deg11_6), y=1-deg.dist11_6, cex=1.2, col=rgb(0,0.4,0.8,alpha=1), ylim=c(0,1), xlim=c(0,127), 
        xlab="2017", ylab="", type = "h")
  mtext("Country 1: Degree centraliy distribution", outer=TRUE, cex=1, font=2)

#plot c2  
  par(mfrow=c(2,3), oma = c(0,2,2,0))
  plot( x=0:max(deg22_1), y=1-deg.dist22_1, cex=1.2, col=rgb(0,0.5,0,alpha=1), ylim=c(0,1), xlim=c(0,64), 
        xlab="1991", ylab="Cumulative Frequency", type = "h")
  plot( x=0:max(deg22_2), y=1-deg.dist22_2,  cex=1.2, col=rgb(0,0.5,0,alpha=1), ylim=c(0,1), xlim=c(0,64), 
        xlab="1996", ylab="", type = "h")
  plot( x=0:max(deg22_3), y=1-deg.dist22_3, cex=1.2, col=rgb(0,0.5,0,alpha=1), ylim=c(0,1), xlim=c(0,64), 
        xlab="2001", ylab="", type = "h")
  plot( x=0:max(deg22_4), y=1-deg.dist22_4, cex=1.2, col=rgb(0,0.5,0,alpha=1), ylim=c(0,1), xlim=c(0,64), 
        xlab="2006", ylab="Cumulative Frequency", type = "h")
  plot( x=0:max(deg22_5), y=1-deg.dist22_5, cex=1.2, col=rgb(0,0.5,0,alpha=1), ylim=c(0,1), xlim=c(0,64), 
        xlab="2011", ylab="", type = "h")
  plot( x=0:max(deg22_6), y=1-deg.dist22_6, cex=1.2, col=rgb(0,0.5,0,alpha=1), ylim=c(0,1), xlim=c(0,64), 
        xlab="2017", ylab="", type = "h")
  mtext("Country 2: Degree centraliy distribution", outer=TRUE, cex=1, font=2)

#plot c1c2
  par(mfrow=c(2,3), oma = c(0,2,2,0))
    plot(x=0:max(deg12_1), y=1-deg.dist12_1, cex=1.2, col=rgb(1,0.2,0.2, alpha=1), ylim=c(0,1), xlim=c(0,126), 
        type="h", xlab="1991", ylab="Cumulative Frequency")
    plot(x=0:max(deg12_2), y=1-deg.dist12_2, cex=1.2, col=rgb(1,0.2,0.2, alpha=1), ylim=c(0,1), xlim=c(0,126), 
          type="h", xlab="1996", ylab="")
    plot(x=0:max(deg12_3), y=1-deg.dist12_3, cex=1.2, col=rgb(1,0.2,0.2, alpha=1), ylim=c(0,1), xlim=c(0,126), 
          type="h", xlab="2001", ylab="")
    plot(x=0:max(deg12_4), y=1-deg.dist12_4, cex=1.2, col=rgb(1,0.2,0.2, alpha=1), ylim=c(0,1), xlim=c(0,126), 
          type="h", xlab="2006", ylab="Cumulative Frequency") 
    plot(x=0:max(deg12_5), y=1-deg.dist12_5, cex=1.2, col=rgb(1,0.2,0.2, alpha=1), ylim=c(0,1), xlim=c(0,126), 
          type="h", xlab="2011", ylab="")
    plot(x=0:max(deg12_6), y=1-deg.dist12_6, cex=1.2, col=rgb(1,0.2,0.2, alpha=1), ylim=c(0,1), xlim=c(0,126), 
          type="h", xlab="2017", ylab="")    
    mtext("Cross-country 1-2: Degree centraliy distribution", outer=TRUE, cex=1, font=2)

remove(deg11_1,deg11_2,deg11_3,deg11_4,deg11_5,deg11_6,
       deg22_1,deg22_2,deg22_3,deg22_4,deg22_5,deg22_6,
       deg12_1,deg12_2,deg12_3,deg12_4,deg12_5,deg12_6,
       deg.dist11_1,deg.dist11_2,deg.dist11_3,deg.dist11_4,deg.dist11_5,deg.dist11_6,
       deg.dist12_1,deg.dist12_2,deg.dist12_3,deg.dist12_4,deg.dist12_5,deg.dist12_6,
       deg.dist22_1,deg.dist22_2,deg.dist22_3,deg.dist22_4,deg.dist22_5,deg.dist22_6       
      )
        
  # Compute and plot degree centrality [normalized] [unweighted] --------------------------------------
    nodes_c1 <- data.frame(ID=regions_c1)
    nodes_c2 <- data.frame(ID=regions_c2)
    nodes_c1c2 <- data.frame(ID=regions_c1c2)
    dc <- data.frame()
    dcb <- data.frame()
    
  #Country 1: Internal links 
  for (i in 1991:2017) {
    s <- paste("select * from par_set where Year = ", i, " and Country_A = 1 and Country_B = 1")
    x <- sqldf(s)
    links <- data.frame(RegionA=x$Region_A, RegionB=x$Region_B)
    net <- graph_from_data_frame(links, nodes_c1, directed = F)
    n <- paste("deg.cent_11_", i, sep="")
    assign(n, centr_degree(net, mode = "all", loops=T, normalized = T))
    d <- data.frame(ID=regions_c1, Year=i, dc=get(n)$res/127) #regions +1  because of autolinks (loops)
    dc <- rbind(dc, d)
  }
  #Country 2: Internal links  
  for (i in 1991:2017) {
    s <- paste("select * from par_set where Year = ", i, " and Country_A = 2 and Country_B = 2")
    x <- sqldf(s)
    links <- data.frame(RegionA=x$Region_A, RegionB=x$Region_B)
    net <- graph_from_data_frame(links, nodes_c2, directed = F)
    n <- paste("deg.cent_22_", i, sep="")
    assign(n, centr_degree(net, mode = "all", loops=T, normalized = T))
    d <- data.frame(ID=regions_c2, Year=i, dc=get(n)$res/64)
    dc <- rbind(dc, d)
  }
  #Country 1-2: External links  
    for (i in 1991:2017) {
      s <- paste("select * from par_set where Year = ", i, " and border = 1")
      x <- sqldf(s)
      links <- data.frame(RegionA=x$Region_A, RegionB=x$Region_B)
      net <- graph_from_data_frame(links, nodes_c1c2, directed = F)
      n <- paste("deg.cent_12_", i, sep="")
      assign(n, centr_degree(net, mode = "all", loops=T, normalized = T))
      d <- data.frame(ID=regions_c1c2, Year=i, dcb=get(n)$res)
      dcb <- rbind(dcb, d)
    }
    dc <- left_join(dc, dcb, by=c("ID","Year"))
    dc <- dc %>%
      mutate(dcb = case_when(
        ID %in% regions_c1 ~ dcb/63,
        ID %in% regions_c2 ~ dcb/126
        )
      )
    remove(x,i,s,links,nodes_c1,nodes_c2,nodes_c1c2,net,n,d,dcb)
  
    # Compile Degree Centrality -----------------------------------------------
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
  #country 1-2: external
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
  remove(year, cn)

    # Plot Degree Centrality [normalized] --------------------------------------------------------------------
  #plot  
  par(old.par)
  plot(degree_centrality11$year,1-degree_centrality11$cn, type = "l",col = "red", xlab = "Year", ylab = "Degree Centrality (Mean)", 
       main = "Degree centrality [normalized / unweighted]", ylim = c(0.4,1),lwd=2) 
  lines(degree_centrality12$year, 1-degree_centrality12$cn, type = "l", col = "blue",lwd=2) 
  lines(degree_centrality22$year, 1-degree_centrality22$cn, type = "l", col = "darkgreen",lwd=2)
  legend(1991, 1, legend=c("C1:Internal", "C2:Internal", "C1-C2:External"),
         col=c("red", "darkgreen","blue"), lty=1, cex=0.8, bty = "n",lwd=2)  
 
  # Compute and plot degree centrality [normalized] [weighted] --------------------------------------
  nodes_c1 <- data.frame(ID=regions_c1)
  nodes_c2 <- data.frame(ID=regions_c2)
  nodes_c1c2 <- data.frame(ID=regions_c1c2)
  dc_w <- data.frame()
  dcb_w <- data.frame()
  
  #Country 1: Internal links 
  for (i in 1991:2017) {
    s <- paste("select * from par_set where Year = ", i, " and c1=1 and c2=0")
    x <- sqldf(s)
    links <- data.frame(RegionA=x$Region_A, RegionB=x$Region_B, w=1+x$Intensity_norm)
    net <- graph_from_data_frame(links, nodes_c1, directed = F)
    n <- paste("deg.cent_11_", i, sep="")
    assign(n, strength(net, mode = "all", loops=T, weights = E(net)$w))
    d <- data.frame(ID=regions_c1, Year=i, dcw=get(n)/127) 
    dc_w <- rbind(dc_w, d)
  }
  #Country 2: Internal links
  nodes <- data.frame(ID=190:378)
  for (i in 1991:2017) {
    s <- paste("select * from par_set where Year = ", i, " and c1=0 and c2=1")
    x <- sqldf(s)
    links <- data.frame(RegionA=x$Region_A, RegionB=x$Region_B, w=1+x$Intensity_norm)
    net <- graph_from_data_frame(links, nodes_c2, directed = F)
    n <- paste("deg.cent_22_", i, sep="")
    assign(n, strength(net, mode = "all", loops=T, weights = E(net)$w))
    d <- data.frame(ID=regions_c2, Year=i, dcw=get(n)/64)
    dc_w <- rbind(dc_w, d)
  }
  #Country 1-2: External links
  nodes <- data.frame(ID=1:378)
  for (i in 1991:2017) {
    s <- paste("select * from par_set where Year = ", i, " and c1=1 and c2=1")
    x <- sqldf(s)
    links <- data.frame(RegionA=x$Region_A, RegionB=x$Region_B, w=1+x$Intensity_norm)
    net <- graph_from_data_frame(links, nodes_c1c2, directed = F)
    n <- paste("deg.cent_12_", i, sep="")
    assign(n, strength(net, mode = "all", loops=T, weights = E(net)$w))
    d <- data.frame(ID=regions_c1c2, Year=i, dcbw=get(n))
    dcb_w <- rbind(dcb_w, d)
  }
  dc_w <- left_join(dc_w, dcb_w, by=c("ID","Year"))
  dc_w <- dc_w %>%
    mutate(dcbw = case_when(
      ID %in% regions_c1 ~ dcbw/63,
      ID %in% regions_c2 ~ dcbw/126
      )
    )  
  remove(i,s,x,links,net,n,nodes_c1,nodes_c2,nodes_c1c2,d,dcb_w)  

    # Compile Degree Centrality [normalized / weighted] -----------------------------------------------
  year <- c(1991:2017)
  #country 1: internal
  degree_centrality11 <- data.frame()
  for (y in year){
    m <- paste("deg.cent_11_",y,sep="" )
    cn <- mean(get(m))/126
    d <- data.frame(Year=y, Degree=cn)
    degree_centrality11 <- rbind(degree_centrality11,d)
  }
  remove(y,m,cn,d)
  #country 1-2: external
  degree_centrality12 <- data.frame()
  for (y in year){
    m <- paste("deg.cent_12_",y,sep="" )
    cn <- mean(get(m))/189
    d <- data.frame(Year=y, Degree=cn)
    degree_centrality12 <- rbind(degree_centrality12,d)
  }
  remove(y,m,cn,d)
  #country 2: internal
  degree_centrality22 <- data.frame()
  for (y in year){
    m <- paste("deg.cent_22_",y,sep="" )
    cn <- mean(get(m))/63
    d <- data.frame(Year=y, Degree=cn)
    degree_centrality22 <- rbind(degree_centrality22,d)
  }
  remove(y,m,cn,d,year)
  
    # Plot Degree Centrality [normalized / weighted] --------------------------------------------------------------------
  par(old.par)
  plot(degree_centrality22$Year,degree_centrality22$Degree, type = "l",col = "red", xlab = "Year", 
       ylab = "Degree Centrality (Mean)", main = "Degree centrality [normalized / weighted]", ylim=c(0,1), lwd=2) 
  lines(degree_centrality11$Year, degree_centrality11$Degree, type = "l", col = "darkgreen", lwd=2)
  lines(degree_centrality12$Year, degree_centrality12$Degree, type = "l", col = "blue", lwd=2)
  legend(1990, 1, legend=c("C1:Internal", "C2:Internal", "C1-C2:External"),
         col=c("darkgreen","red", "blue"), lty=1, cex=0.8, bty="n", lwd=2)  
  
  # Compute Harmonic Centrality [normalized / unweighted] ---------------------------------------------
  #C1: internal
  hc <- data.frame()
  hcb <- data.frame()
  
  harm_cent11 <- data.frame()
    for (i in 1991:2017){
    n <- paste("net_11_", i, sep="")
    x <- get(n)
    h <- paste("harm.cent_11_", i, sep="")
    assign(h, harmonic_centrality(x, mode = "all", weights = NULL)/126)
    d <- data.frame(year=i, harmonic_centrality= mean(get(h)))
    harm_cent11 <- rbind(harm_cent11, d)
    d1 <- data.frame(ID=regions_c1,Year=i, hc=get(h))
    hc <- rbind(hc,d1)
  }
    #C1-C2: external
    harm_cent12 <- data.frame()
    for (i in 1991:2017){
      n <- paste("net_12_", i, sep="")
      x <- get(n)
      h <- paste("harm.cent_12_", i, sep="")
      assign(h, harmonic_centrality(x, mode = "all", weights = NULL))
      d <- data.frame(year=i, harmonic_centrality=mean(get(h)))
      harm_cent12 <- rbind(harm_cent12, d)
      d1 <- data.frame(ID=regions_c1c2,Year=i, hcb=get(h))
      hcb <- rbind(hcb,d1)      
    }
      #C2: internal
      harm_cent22 <- data.frame()
      for (i in 1991:2017){
        n <- paste("net_22_", i, sep="")
        x <- get(n)
        h <- paste("harm.cent_22_", i, sep="")
        assign(h, harmonic_centrality(x, mode = "all", weights = NULL)/63)
        d <- data.frame(year=i, harmonic_centrality=mean(get(h)))
        harm_cent22 <- rbind(harm_cent22, d)
        d1 <- data.frame(ID=regions_c2,Year=i, hc=get(h))
        hc <- rbind(hc,d1)        
      }
      hc <- left_join(hc, hcb, by=c("ID","Year"))
      hc <- hc %>%
        mutate(hcb = case_when(
          ID %in% regions_c1 ~ hcb/189, #same weight because harmonic centraliy contains indirect conncetions
          ID %in% regions_c2 ~ hcb/189
          )
        )      
      remove(x,i,n,h,d,d1,hcb)

      #plot 1  
      par(old.par)
        plot(harm_cent22, type = "l",col = "red", xlab = "Year", ylab = "Harmonic Centrality (Mean)", 
             main = "Harmonic centrality [normalized / unweighted]", ylim=c(0.4,1),lwd=2) 
        lines(harm_cent11, type = "l", col = "darkgreen",lwd=2)
        lines(harm_cent12, type = "l", col = "blue",lwd=2)
        legend(1990, 1, legend=c("C1:Internal","C2:Internal", "C1-C2:External"),
               col=c("darkgreen","red", "blue"), lty=NA, cex=0.8, bty="n", lwd=10, pch=21) 

        #density plots             
        b <- adjustcolor("blue", alpha.f = 0.2)
        r <- adjustcolor("red", alpha.f = 0.2)
        g <- adjustcolor("green", alpha.f = 0.2)
        par(mfrow=c(2,2), oma=c(0,0,2,0))
        plot(density(harm.cent_11_1991), xlab="", ylab="# nodes", main="Country 1: Internal links", cex=0.8, 
             col="red", xlim = c(0,1),ylim=c(0,40)) 
        polygon(density(harm.cent_11_1991), col=r, border="black")
        lines(density(harm.cent_11_2001), col="blue")
        polygon(density(harm.cent_11_2001), col=b, border="black")
        lines(density(harm.cent_11_2017), col="green")
        polygon(density(harm.cent_11_2017), col=g, border="black")
        
        plot(density(harm.cent_12_1991), xlab="", ylab="", main="Country 1-2: External links", cex=0.8, 
             col="red", xlim = c(0,1),ylim=c(0,40)) 
        polygon(density(harm.cent_12_1991), col=r, border="black")
        lines(density(harm.cent_12_2001), col="blue")
        polygon(density(harm.cent_12_2001), col=b, border="black")
        lines(density(harm.cent_12_2017), col="green")
        polygon(density(harm.cent_12_2017), col=g, border="black")
        
        plot(density(harm.cent_22_1991), xlab="", ylab="# nodes", main="Country 2: Internal links", cex=0.8, 
             col="red", xlim = c(0,1),ylim=c(0,40)) 
        polygon(density(harm.cent_22_1991), col=r, border="black")
        lines(density(harm.cent_22_2001), col="blue")
        polygon(density(harm.cent_22_2001), col=b, border="black")
        lines(density(harm.cent_22_2017), col="green")
        polygon(density(harm.cent_22_2017), col=g, border="black")
        legend(1.2, 40, legend=c("1991","2001", "2017"),
               col=c(r,b,g), lty=NA, pch=21, cex=1.2, bty = "n", lwd=32, xpd = "NA") 
        mtext("Harmonic Centrality [normalized / unweighted]", 
              outer=TRUE, cex=1, font=2)
        
  # Compute Harmonic Centrality [normalized / weighted] ---------------------------------------------
        #C1: internal
        hc_w <- data.frame()
        hcb_w <- data.frame()
        
        harm_cent11 <- data.frame()
        for (y in 1991:2017){
          n <- paste("net_11_", y, sep="")
          x <- get(n)
          h <- paste("harm.cent_11_", y, sep="")
          assign(h, harmonic_centrality(x, mode = "all", weights =1+E(x)$w2)/126)
          d <- data.frame(year=y, harmonic_centrality=mean(get(h)))
          harm_cent11 <- rbind(harm_cent11, d)
          d1 <- data.frame(ID=regions_c1,Year=y, hcw=get(h))
          hc_w <- rbind(hc_w,d1)          
        }
        #C1-C2: external
        harm_cent12 <- data.frame()
        for (y in 1991:2017){
          n <- paste("net_12_", y, sep="")
          x <- get(n)
          h <- paste("harm.cent_12_", y, sep="")
          assign(h, harmonic_centrality(x, mode = "all", weights =1+E(x)$w2))
          d <- data.frame(year=y, harmonic_centrality=mean(get(h)))
          harm_cent12 <- rbind(harm_cent12, d)
          d1 <- data.frame(ID=regions_c1c2,Year=y, hcbw=get(h))
          hcb_w <- rbind(hcb_w,d1)            
        }
        #C2: internal
        harm_cent22 <- data.frame()
        for (y in 1991:2017){
          n <- paste("net_22_", y, sep="")
          x <- get(n)
          h <- paste("harm.cent_22_", y, sep="")
          assign(h, harmonic_centrality(x, mode = "all", weights =1+E(x)$w2)/63)
          d <- data.frame(year=y, harmonic_centrality=mean(get(h)))
          harm_cent22 <- rbind(harm_cent22, d)
          d1 <- data.frame(ID=regions_c2,Year=y, hcw=get(h))
          hc_w <- rbind(hc_w,d1)            
        }
        hc_w <- left_join(hc_w, hcb_w, by=c("ID","Year"))
        hc_w <- hc_w %>%
          mutate(hcbw = case_when(
            ID %in% regions_c1 ~ hcbw/189,
            ID %in% regions_c2 ~ hcbw/189
            )
          )      
        remove(x,y,n,h,d,d1,hcb_w)
        
        #plot 1  
        par(old.par)
        plot(harm_cent22, type = "l",col = "red", xlab = "Year", ylab = "Harmonic Centrality (Mean)", 
             main = "Harmonic centrality [normalized / weighted]", ylim=c(0.4,1),lwd=2) 
        lines(harm_cent11, type = "l", col = "darkgreen",lwd=2)
        lines(harm_cent12, type = "l", col = "blue",lwd=2)
        legend(1990, 1, legend=c("C1:Internal","C2:Internal", "C1-C2:External"),
               col=c("darkgreen","red", "blue"), lty=NA, cex=0.8, bty="n", lwd=10, pch=21) 
        
        #density plots             
        b <- adjustcolor("blue", alpha.f = 0.2)
        r <- adjustcolor("red", alpha.f = 0.2)
        g <- adjustcolor("green", alpha.f = 0.2)
        par(mfrow=c(2,2), oma=c(0,0,2,0))
        plot(density(harm.cent_11_1991), xlab="", ylab="# nodes", main="Country 1: Internal links", cex=0.8, 
             col="red", xlim = c(0,1),ylim=c(0,40)) 
        polygon(density(harm.cent_11_1991), col=r, border="black")
        lines(density(harm.cent_11_2001), col="blue")
        polygon(density(harm.cent_11_2001), col=b, border="black")
        lines(density(harm.cent_11_2017), col="green")
        polygon(density(harm.cent_11_2017), col=g, border="black")
        
        plot(density(harm.cent_12_1991), xlab="", ylab="", main="Country 1-2: External links", cex=0.8, 
             col="red", xlim = c(0,1),ylim=c(0,40)) 
        polygon(density(harm.cent_12_1991), col=r, border="black")
        lines(density(harm.cent_12_2001), col="blue")
        polygon(density(harm.cent_12_2001), col=b, border="black")
        lines(density(harm.cent_12_2017), col="green")
        polygon(density(harm.cent_12_2017), col=g, border="black")
        
        plot(density(harm.cent_22_1991), xlab="", ylab="# nodes", main="Country 2: Internal links", cex=0.8, 
             col="red", xlim = c(0,1),ylim=c(0,40)) 
        polygon(density(harm.cent_22_1991), col=r, border="black")
        lines(density(harm.cent_22_2001), col="blue")
        polygon(density(harm.cent_22_2001), col=b, border="black")
        lines(density(harm.cent_22_2017), col="green")
        polygon(density(harm.cent_22_2017), col=g, border="black")
        legend(1.2, 40, legend=c("1991","2001", "2017"),
               col=c(r,b,g), lty=NA, pch=21, cex=1.2, bty = "n", lwd=32, xpd = "NA") 
        mtext("Harmonic Centrality [normalized / weighted]", 
              outer=TRUE, cex=1, font=2)
        
  # Fixed effects regressions ------------------------------------------------
    #Intensity ~ border (with year and link) 
        #partial set 
        pyl <- plm(Intensity_norm ~ border,
                  data = par_set, 
                  index = c("edge1", "Year"),
                  model = "within")
        coeftest(pyl)
        #full set
        fyl <- plm(Intensity_norm ~ border, 
                  data = full_set,
                  index = c("edge1", "Year"), 
                  model = "within")
        coeftest(fyl)

    #Degree centrality ~ border [unweighted] 
        dc_py <- plm(dcb ~ dc, 
                   data = dc,
                   index = c("ID","Year"),
                   model = "within")
      coeftest(dc_py)
      
    #Degree centrality ~ border [weighted]
      dc_w_py <- plm(dcbw ~ dcw, 
                  data = dc_w, 
                  index = c("ID","Year"),
                  model = "within")
      coeftest(dc_w_py)     

#fixed effect regression on harmonic centrality works but the results can not easily or even 
#meaninfully be interpreted.
    #Harmonic centrality ~ border [unweighted] 
      hc_py <- plm(hcb ~ hc, 
                   data = hc,
                   index = c("ID","Year"),
                   model = "within")
      coeftest(hc_py)
      
    #Degree centrality ~ border [weighted]
      hc_w_py <- plm(hcbw ~ hcw, 
                     data = hc_w, 
                     index = c("Year"),
                     model = "within")
      coeftest(hc_w_py)        


#### THE CODE GRAVEYARD / SCRIPTS OF THE DEAD -----------------------------------------------------------
      # TEST duplicate rows full set [No DUPLICATES!!!] -------------------------------------------------
    test_doubles <- data.frame()
    for (c1 in 1:2){
      for (c2 in 1:2){
        for (y in 1991:2017) {
          s <- paste("select count(*) 
                 from full_set where Year=",y,
                     " and Country_A=",c1,
                     " and Country_B=",c2, 
                     " and autolink=0", sep="")
          r <- as.numeric(sqldf(s))
          f <- paste("select edge1 from full_set where Year=",y,
                     " and Country_A=",c1,
                     " and Country_B=",c2,
                     " and autolink=0
                 except 
                 select edge2 from full_set where Year=",y,
                     " and Country_A=",c1,
                     " and Country_B=",c2,
                     " and autolink=0", sep="")
          t <- sqldf(f)
          z <- data.frame(Year=y, CountryA=c1, no_doubles=nrow(t)==r)
          test_doubles <- rbind(test_doubles, z)
        }
      }    
    }
    paste("Data has duplicate rows: ", FALSE %in% test_doubles$no_doubles)  
    remove(c1,c2,f,y,r,s,t,z, test_doubles)
    
    
      # TEST duplicate rows partial set [NO DUPLICATE ROWS] -------------------------------------------------
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
    paste("Data has duplicate row: ", FALSE %in% test_doubles$no_doubles)  
    remove(c1,c2,f,y,r,s,t,z, test_doubles)
    
    
      # Create full set ---------------------------------------------------------
    regions <- data.frame()
    for (a in 1:189){
      for (i in 1:189){
        c <- expand.grid(i,a)
        if(c$Var1 != c$Var2) {
          if(c$Var1 > c$Var2){
            regions <- rbind(regions, c)        
          }
        }
      }
    }
    remove(a,i,c)
    #test duplicates
    regions$edge1 <- paste(regions$Var1, "-", regions$Var2, sep="")
    regions$edge2 <- paste(regions$Var2, "-", regions$Var1, sep="")
    test <- sqldf("select edge1 from regions except select edge2 from regions")
    paste("region list has no duplicate rows:", nrow(regions)==nrow(test))
    remove(test)
    
    #add rows to par_set
    full_set <- select (par_set,c(Year, Region_A, c1, Region_B, c2, 
                                  Intensity))
    for (y in 1991:2017) {
      t <- paste("select edge1 from regions  
              except 
                select edge1 from par_set 
                  where Year=",y, 
                 " and c1=1 and c2=0 
                  union 
                  select edge2 as edge1 from par_set
                    where Year=",y,
                 " and c1=1 and c2=0",sep="")
      test <- sqldf(t)
      a <- paste("select ",y," as Year, Var1 as Region_A, 1 as c1, 
                Var2 as Region_B, 0 as c2,0 as Intensity
                from regions where edge1 in 
                (select edge1 from test)")
      addrow <- sqldf(a)
      full_set <- rbind(full_set, addrow)
    }
    for (y in 1991:2017) {
      t <- paste("select edge1 from regions  
              except 
                select edge1 from par_set 
                  where Year=",y, 
                 " and c1=0 and c2=1 
                  union 
                  select edge2 as edge1 from par_set
                    where Year=",y,
                 " and c1=0 and c2=1",sep="")
      test <- sqldf(t)
      a <- paste("select ",y," as Year, Var1 as Region_A, 0 as c1, 
                Var2 as Region_B, 1 as c2,0 as Intensity
                from regions where edge1 in 
                (select edge1 from test)")
      addrow <- sqldf(a)
      full_set <- rbind(full_set, addrow)
    }
    for (y in 1991:2017) {
      t <- paste("select edge1 from regions  
              except 
                select edge1 from par_set 
                  where Year=",y, 
                 " and c1=1 and c2=1 
                  union 
                  select edge2 as edge1 from par_set
                    where Year=",y,
                 " and c1=1 and c2=1",sep="")
      test <- sqldf(t)
      a <- paste("select ",y," as Year, Var1 as Region_A, 1 as c1, 
                Var2 as Region_B, 1 as c2,0 as Intensity
                from regions where edge1 in 
                (select edge1 from test)")
      addrow <- sqldf(a)
      full_set <- rbind(full_set, addrow)
    }
    remove(t,y,a,test,addrow)
    
    full_set$edge1 <- paste(full_set$Region_A, "-", full_set$Region_B, sep="")
    full_set$edge2 <- paste(full_set$Region_B, "-", full_set$Region_A, sep="")
    
    full_set <- full_set %>%
      mutate(autolink=ifelse(c1+c2<=1, 1, 0)) %>%
      mutate(border=ifelse(c1+c2>=2, 1, 0))
    
    full_set$Intensity_norm = (
      full_set$Intensity-min(full_set$Intensity)) / 
      (max(full_set$Intensity)-min(full_set$Intensity)
      )
    cor(full_set$Intensity, full_set$Intensity_norm) == 1
    
      # Create full set with country---------------------------------------------------------
    regions <- data.frame()
    for (a in 1:189){
      for (i in 1:189){
        c <- expand.grid(i,a)
        if(c$Var1 != c$Var2) {
          if(c$Var1 > c$Var2){
            regions <- rbind(regions, c)        
          }
        }
      }
    }
    remove(a,i,c)
    #test duplicates
    regions$edge1 <- paste(regions$Var1, "-", regions$Var2, sep="")
    regions$edge2 <- paste(regions$Var2, "-", regions$Var1, sep="")
    test <- sqldf("select edge1 from regions except select edge2 from regions")
    paste("region list has no duplicate rows:", nrow(regions)==nrow(test))
    remove(test)
    
    #add rows to par_set
    full_set <- select (par_set,c(Year, Region_A, Country_A, Region_B, Country_B, Intensity))
    #country 1 internal
    for (y in 1991:2017) {
      t <- paste("select edge1 from regions  
              except 
                select edge1 from par_set 
                  where Year=",y, 
                 " and Country_A=1 and Country_B=1 
                  union 
                  select edge2 as edge1 from par_set
                    where Year=",y,
                 " and Country_A=1 and Country_B=1",sep="")
      test <- sqldf(t)
      a <- paste("select ",y," as Year, Var1 as Region_A, 1 as Country_A, 
                Var2 as Region_B, 1 as Country_B,0 as Intensity
                from regions where edge1 in 
                (select edge1 from test)")
      addrow <- sqldf(a)
      full_set <- rbind(full_set, addrow)
    }
    #country 2 internal
    for (y in 1991:2017) {
      t <- paste("select edge1 from regions  
              except 
                select edge1 from par_set 
                  where Year=",y, 
                 " and Country_A=2 and Country_B=2 
                  union 
                  select edge2 as edge1 from par_set
                    where Year=",y,
                 " and Country_A=2 and Country_B=2",sep="")
      test <- sqldf(t)
      a <- paste("select ",y," as Year, Var1 as Region_A, 2 as Country_A, 
                Var2 as Region_B, 2 as Country_B,0 as Intensity
                from regions where edge1 in 
                (select edge1 from test)")
      addrow <- sqldf(a)
      full_set <- rbind(full_set, addrow)
    }
    #country 1-2 external
    for (y in 1991:2017) {
      t <- paste("select edge1 from regions  
              except 
                select edge1 from par_set 
                  where Year=",y, 
                 " and Country_A=1 and Country_B=2 
                  union 
                  select edge2 as edge1 from par_set
                    where Year=",y,
                 " and Country_A=1 and Country_B=2",sep="")
      test <- sqldf(t)
      a <- paste("select ",y," as Year, Var1 as Region_A, 1 as Country_A, 
                Var2 as Region_B, 2 as Country_B,0 as Intensity
                from regions where edge1 in 
                (select edge1 from test)")
      addrow <- sqldf(a)
      full_set <- rbind(full_set, addrow)
    }
    #country 2-1 external
    for (y in 1991:2017) {
      t <- paste("select edge1 from regions  
              except 
                select edge1 from par_set 
                  where Year=",y, 
                 " and Country_A=2 and Country_B=1 
                  union 
                  select edge2 as edge1 from par_set
                    where Year=",y,
                 " and Country_A=2 and Country_B=1",sep="")
      test <- sqldf(t)
      a <- paste("select ",y," as Year, Var1 as Region_A, 2 as Country_A, 
                Var2 as Region_B, 1 as Country_B,0 as Intensity
                from regions where edge1 in 
                (select edge1 from test)")
      addrow <- sqldf(a)
      full_set <- rbind(full_set, addrow)
    }
    remove(t,y,a,test,addrow)
    
    full_set <- full_set %>%
      mutate(autolink=ifelse(Country_A==Country_B & Region_A==Region_B, 1, 0)) %>%
      mutate(border=ifelse(Country_A!=Country_B, 1, 0))
    
    full_set$edge1 <- paste(full_set$Region_A, "-", full_set$Region_B, sep="")
    full_set$edge2 <- paste(full_set$Region_B, "-", full_set$Region_A, sep="")
    
    #test duplicates
    test_doubles <- data.frame()
    for (c1 in 1:2){
      for (c2 in 1:2){
        for (y in 1991:2017) {
          s <- paste("select count(*) 
                 from full_set where Year=",y,
                     " and Country_A=",c1,
                     " and Country_B=",c2, 
                     " and autolink=0", sep="")
          r <- as.numeric(sqldf(s))
          f <- paste("select edge1 from full_set where Year=",y,
                     " and Country_A=",c1,
                     " and Country_B=",c2,
                     " and autolink=0
                 except 
                 select edge2 from full_set where Year=",y,
                     " and Country_A=",c1,
                     " and Country_B=",c2,
                     " and autolink=0", sep="")
          t <- sqldf(f)
          z <- data.frame(Year=y, CountryA=c1, no_doubles=nrow(t)==r)
          test_doubles <- rbind(test_doubles, z)
        }
      }    
    }
    paste("Data has duplicate row: ", FALSE %in% test_doubles$no_doubles)  
    remove(c1,c2,f,y,r,s,t,z, test_doubles)
    
    (27*(189*(189)/2))*3
    
    
    full_set$Intensity_norm = (
      full_set$Intensity-min(full_set$Intensity)) / 
      (max(full_set$Intensity)-min(full_set$Intensity)
      )
    cor(full_set$Intensity, full_set$Intensity_norm) == 1
    
    
      # Construct adjacency matrix per year, country, links [directed] ------------------------------------------------
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
      
      #For all links
      nodes <- data.frame(ID=1:378)
      dc <- data.frame()
      for (b in 0:1){
        for (i in 1991:2017){
          s <- paste("select * from par_set where Year = ", i,"and border=",b)
          x <- sqldf(s)
          links <- data.frame(RegionA=x$rid1, RegionB=x$rid2)
          net <- graph_from_data_frame(links, nodes, directed = F)
          n <- paste("deg.cent_", i, sep="")
          assign(n, strength(net, mode = "all", loops=T))
          t <- data.frame(rid1=1:378, rid2=1:378, Year=i, dc_w=get(n), border=b)
          dc <- rbind(dc, t)
          remove(i,s,x,links,net,n,t)
        }
      }  
      remove(nodes)  
      
      dc_set <- left_join(par_set, dc, by=c("rid1", "Year", "border"))
      dc_set$rid2 <- dc_set$rid2.x
      dc_set <- left_join(dc_set, dc, by=c("rid2", "Year", "border"))
      dc_set <- dc_set %>% 
        mutate(dc=ifelse(autolink==0, dc_w.x+dc_w.y, dc_w.x))
      
      
      
      
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
      
      # Compute closeness centrality (not recommended) --------------------------------------
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
      
      # Compile CLoseness Centrality (not recommended) -----------------------------------------------
      
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
      
      # Plot (not recommended) --------------------------------------------------------------------
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
      
      # Newman-Garvin clustering (not very usefull) ------------------------------------------------
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
      
      
      
      
      # Compute and plot eigenvector centrality (needs revision) ---------------------------------
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
      
      # Compile eigenvector centrality (needs revision) -----------------------------------------------------------------
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
      
      # Plot eigenvector centrality (needs revision) --------------------------------------------------------------------
      #plot
      plot(eigenvalue22,type = "o",col = "red", xlab = "Year", ylab = "Centrality", 
           main = "Eigenvector centrality (partial)") 
      lines(eigenvalue21, type = "o", col = "blue") 
      lines(eigenvalue11, type = "o", col = "darkgreen")
      lines(eigenvalue12, type = "o", col = "orange")
      legend(1991, 1.5, legend=c("C2:Internal", "C2:External", "C1:Internal", "C1:External"),
             col=c("red", "blue", "darkgreen", "orange"), lty=1:2, cex=0.6)
      
      
    
      # Construct nework graphs / adjacency matrix per year, country (works) ------------------------------------------------
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
      
      
      # Other stuff--------------------------------------------------
      #plot 2  
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
      #sorted histogram plot 
      #country 1
      t1a <- sort(harm.cent_11_1991)
      t1a <- data.frame(num=1:126, harmonic_centrality=t1a)
      t2a <- sort(harm.cent_11_1996)
      t2a <- data.frame(num=1:126, harmonic_centrality=t2a)
      t3a <- sort(harm.cent_11_2001)
      t3a <- data.frame(num=1:126, harmonic_centrality=t3a)
      t4a <- sort(harm.cent_11_2006)
      t4a <- data.frame(num=1:126, harmonic_centrality=t4a)
      t5a <- sort(harm.cent_11_2011)
      t5a <- data.frame(num=1:126, harmonic_centrality=t5a)
      t6a <- sort(harm.cent_11_2017)
      t6a <- data.frame(num=1:126, harmonic_centrality=t6a)
      
      par(mfrow=c(2,3), oma = c(0,2,2,0))
      plot(t1a, ylim=c(0,210), xlab="1991", ylab="Harmonic Centrality", type="h", col=b)
      plot(t2a, ylim=c(0,210), xlab="1996", ylab="", type="h", col=adjustcolor("blue", alpha.f = 0.2))
      plot(t3a, ylim=c(0,210), xlab="2001", ylab="", type="h", col=adjustcolor("blue", alpha.f = 0.2))
      plot(t4a, ylim=c(0,210), xlab="2006", ylab="Harmonic Centrality", type="h", col=adjustcolor("blue", alpha.f = 0.2))
      plot(t5a, ylim=c(0,210), xlab="2011", ylab="", type="h", col=adjustcolor("blue", alpha.f = 0.2))
      plot(t6a, ylim=c(0,210), xlab="2017", ylab="", type="h", col=adjustcolor("blue", alpha.f = 0.2))
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
      
      