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
  
  #This block of code checks and prepares the raw data for analyis. 
  #We call this set a partial data set because regions with Intenstiy==0 are still missing
  
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
  
  #add dummies per regions 
    for (i in 1:189){
      par_set[[paste0("r_",i)]] <- ifelse(par_set$Region_A==i|par_set$Region_B==i,1,0)
    }

  # Create full data set ---------------------------------------------------------
  
  #This block of code takes the partial data set and checks which combinations of regions are missing 
  #It adds those links to create a full data set 
  
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
  
  #add dummies per regions 
    for (i in 1:189){
      full_set[[paste0("r_",i)]] <- ifelse(full_set$Region_A==i|full_set$Region_B==i,1,0)
    }
      
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
  #This code block will plot the link intesity per country (internal/external) per year using the full and partial data set
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
          par(old.par)
          plot(fc1c1,type = "l",col = "darkgreen", xlab = "Year", ylab = "", 
               main = "Link Intensity (Mean)", ylim=c(0, 60000), lty=2, lwd=2)
          lines(pc1c1, type = "l", col = "darkgreen", lty=1, lwd=2)
          lines(fc1c2, type = "l", col = "blue", lty=2, lwd=2)
          lines(pc1c2, type = "l", col = "blue", lty=1, lwd=2)
          lines(fc2c2, type = "l", col = "red", lty=2, lwd=2)
          lines(pc2c2, type = "l", col = "red", lty=1, lwd=2)
          legend(1991, 60000, legend=c("C1:Internal (partial)", "C1:Internal (full)",
                                       "C2:Internal (partial)", "C2:Internal (full)",
                                       "C1-C2:External (partial)","C1-C2:External (full"),
                 col=c("darkgreen","darkgreen","red","red","blue","blue"), lty=1:2, cex=0.8, bty="n", lwd=2)
                  remove(fc1c1, fc2c2, fc1c2, pc1c1, pc1c2, pc2c2)
        
  # Count zeros --------------------------------------------
  # This block of code counts the number of links that have Intensity==0 and plots them per year
                
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

   #plot
   plot(zc2c2$Year, zc2c2$share, type = "l",col = "red", xlab = "Year", ylab = "Share", 
        main = "Share of regions with no connections", ylim=c(0,0.2), lwd=2) 
   lines(zc1c1$Year,zc1c1$share, type = "l", col = "darkgreen", lwd=2)
   lines(zc1c2$Year,zc1c2$share, type = "l", col = "blue", lwd=2)
   legend(2008, 0.2, legend=c("C1:Internal", "C2:Internal", "C1-C2:External"),
          col=c("darkgreen", "red", "blue"), lty=1, cex=0.8, bty="n", lwd=2)
   remove(zc1c1, zc1c2, zc2c2)
   
  # Degree distibution plots [not normalized] [unweighted] ------------------------------------------------
  # This block of code computes regional degree centraliy and plots the distribution across time
   
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
  
  # plots
    par(mfrow=c(2,3), oma = c(0,2,2,0))  
    plot(sort(deg11_1/126), xlab="regions", ylab="degree centrality", ylim=c(0,1), col="darkgreen", main="1991")
    plot(sort(deg11_2/126), xlab="regions", ylab="", ylim=c(0,1), col="darkgreen", main="1996")
    plot(sort(deg11_3/126), xlab="regions", ylab="", ylim=c(0,1), col="darkgreen", main="2001")
    plot(sort(deg11_4/126), xlab="regions", ylab="degree centrality", ylim=c(0,1), col="darkgreen", main="2006")
    plot(sort(deg11_5/126), xlab="regions", ylab="", ylim=c(0,1), col="darkgreen", main="2011")
    plot(sort(deg11_6/126), xlab="regions", ylab="", ylim=c(0,1), col="darkgreen", main="2017")
    mtext("", outer=TRUE, cex=1, font=2)
    
    par(mfrow=c(2,3), oma = c(0,2,2,0))  
    plot(sort(deg22_1/63), xlab="regions", ylab="degree centrality", ylim=c(0,1), col="red", main="1991")
    plot(sort(deg22_2/63), xlab="regions", ylab="", ylim=c(0,1), col="red", main="1996")
    plot(sort(deg22_3/63), xlab="regions", ylab="", ylim=c(0,1), col="red", main="2001")
    plot(sort(deg22_4/63), xlab="regions", ylab="degree centrality", ylim=c(0,1), col="red", main="2006")
    plot(sort(deg22_5/63), xlab="regions", ylab="", ylim=c(0,1), col="red", main="2011")
    plot(sort(deg22_6/63), xlab="regions", ylab="", ylim=c(0,1), col="red", main="2017")
    mtext("", outer=TRUE, cex=1, font=2)
    
    intx <- c(rep(63,22),rep(126,77),rep(63,41), rep(126,49))
    
    intf1 <- data.frame(ID=1:189, dc=deg12_1/intx)
    intf1$colour="darkgreen"
    intf1$colour[intf1$ID %in% regions_c2]="red"
    intf2 <- data.frame(ID=1:189, dc=deg12_2/intx)
    intf2$colour="darkgreen"
    intf2$colour[intf1$ID %in% regions_c2]="red"  
    intf3 <- data.frame(ID=1:189, dc=deg12_3/intx)
    intf3$colour="darkgreen"
    intf3$colour[intf1$ID %in% regions_c2]="red"
    intf4 <- data.frame(ID=1:189, dc=deg12_4/intx)
    intf4$colour="darkgreen"
    intf4$colour[intf1$ID %in% regions_c2]="red"
    intf5 <- data.frame(ID=1:189, dc=deg12_5/intx)
    intf5$colour="darkgreen"
    intf5$colour[intf1$ID %in% regions_c2]="red"
    intf6 <- data.frame(ID=1:189, dc=deg12_6/intx)
    intf6$colour="darkgreen"
    intf6$colour[intf1$ID %in% regions_c2]="red"
    
    
    par(mfrow=c(2,3), oma = c(0,2,2,0))  
    plot(sort(intf1$dc), xlab="", ylab="degree centrality", ylim=c(0,1), col=intf1$colour, main="1991")
    legend(0, 1, legend=c("Country 1","Country 2"),
           col=c("darkgreen","red"), lty=NA, cex=1, bty="n", lwd=10, pch=21) 
    plot(sort(intf2$dc), xlab="regions", ylab="", ylim=c(0,1), col=intf1$colour, main="1996")
    plot(sort(intf3$dc), xlab="regions", ylab="", ylim=c(0,1), col=intf1$colour, main="2001")
    plot(sort(intf4$dc), xlab="regions", ylab="degree centrality", ylim=c(0,1), col=intf1$colour, main="2006")
    plot(sort(intf5$dc), xlab="regions", ylab="", ylim=c(0,1), col=intf1$colour, main="2011")
    plot(sort(intf6$dc), xlab="regions", ylab="", ylim=c(0,1), col=intf1$colour, main="2017")
    mtext("", outer=TRUE, cex=1, font=2)
    
    remove(intx)
    
    remove(deg11_1,deg11_2,deg11_3,deg11_4,deg11_5,deg11_6,
           deg22_1,deg22_2,deg22_3,deg22_4,deg22_5,deg22_6,
           deg12_1,deg12_2,deg12_3,deg12_4,deg12_5,deg12_6,
           deg.dist11_1,deg.dist11_2,deg.dist11_3,deg.dist11_4,deg.dist11_5,deg.dist11_6,
           deg.dist12_1,deg.dist12_2,deg.dist12_3,deg.dist12_4,deg.dist12_5,deg.dist12_6,
           deg.dist22_1,deg.dist22_2,deg.dist22_3,deg.dist22_4,deg.dist22_5,deg.dist22_6       
          )
        
  # Compute and plot degree centrality [normalized] [unweighted] --------------------------------------
  #This block of code computes the normalized degree centrality per network per year    
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
    #This block of code compiles the above calculated degree centrality into a single data object
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
    #This code plots the compiled normalized and unweighted degree centrality per year 
    par(old.par)
    plot(degree_centrality11$year,1-degree_centrality11$cn, type = "l",col = "darkgreen", xlab = "Year", ylab = "Degree Centrality (Mean)", 
         main = "Degree centrality [normalized / unweighted]", ylim = c(0,1),lwd=2) 
    lines(degree_centrality12$year, 1-degree_centrality12$cn, type = "l", col = "blue",lwd=2) 
    lines(degree_centrality22$year, 1-degree_centrality22$cn, type = "l", col = "red",lwd=2)
    legend(1991, 1, legend=c("C1:Internal", "C2:Internal", "C1-C2:External"),
           col=c("darkgreen", "red","blue"), lty=1, cex=0.8, bty = "n",lwd=2)  
 
  # Compute and plot degree centrality [normalized] [weighted] --------------------------------------
  #This block of code computes the normalized degree centrality per network per year    
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
    #This block of code compiles the above calculated degree centrality into a single data object
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
    #This code plots the compiled normalized and weighted degree centrality per year 
    par(old.par)
    plot(degree_centrality22$Year,degree_centrality22$Degree, type = "l",col = "red", xlab = "Year", 
         ylab = "Degree Centrality (Mean)", main = "Degree centrality [normalized / weighted]", ylim=c(0,1), lwd=2) 
    lines(degree_centrality11$Year, degree_centrality11$Degree, type = "l", col = "darkgreen", lwd=2)
    lines(degree_centrality12$Year, degree_centrality12$Degree, type = "l", col = "blue", lwd=2)
    legend(1990, 1, legend=c("C1:Internal", "C2:Internal", "C1-C2:External"),
           col=c("darkgreen","red", "blue"), lty=1, cex=0.8, bty="n", lwd=2)  
  
  # Compute Harmonic Centrality [normalized / unweighted] ---------------------------------------------
  #This block of code computes the normalized harmoic centrality per network per year    
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
      assign(h, harmonic_centrality(x, mode = "all", weights = NULL)/189)
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
      remove(x,i,n,h,d,d1,hcb)

      #This code plots the compiled normalized and unweighted harmonic centrality per year 
      par(old.par)
        plot(harm_cent22, type = "l",col = "red", xlab = "Year", ylab = "Harmonic Centrality (Mean)", 
             main = "Harmonic centrality [normalized / unweighted]", ylim=c(0.4,1),lwd=2) 
        lines(harm_cent11, type = "l", col = "darkgreen",lwd=2)
        lines(harm_cent12, type = "l", col = "blue",lwd=2)
        legend(1990, 1, legend=c("C1:Internal","C2:Internal", "C1-C2:External"),
               col=c("darkgreen","red", "blue"), lty=NA, cex=0.8, bty="n", lwd=10, pch=21) 

        #This code create density plots for the harmonic centraliy in 1991,2001, and 2017            
        b <- adjustcolor("blue", alpha.f = 0.2)
        r <- adjustcolor("red", alpha.f = 0.2)
        g <- adjustcolor("green", alpha.f = 0.2)
        par(mfrow=c(3,1), oma=c(0,0,2,0))
        plot(density(harm.cent_11_1991), xlab="", ylab="# of nodes", main="Country 1: Internal links", cex=0.8, 
             col="red", xlim = c(0,1),ylim=c(0,40)) 
        legend(0, 53, legend=c("1991","2001", "2017"),
               col=c(r,b,g), lty=NA, pch=21, cex=1, bty = "n", lwd=20, xpd = "NA") 
        polygon(density(harm.cent_11_1991), col=r, border="black")
        lines(density(harm.cent_11_2001), col="blue")
        polygon(density(harm.cent_11_2001), col=b, border="black")
        lines(density(harm.cent_11_2017), col="green")
        polygon(density(harm.cent_11_2017), col=g, border="black")
        
        plot(density(harm.cent_22_1991), xlab="", ylab="# of nodes", main="Country 2: Internal links", cex=0.8, 
             col="red", xlim = c(0,1),ylim=c(0,40))
        legend(0, 53, legend=c("1991","2001", "2017"),
               col=c(r,b,g), lty=NA, pch=21, cex=1, bty = "n", lwd=20, xpd = "NA") 
        polygon(density(harm.cent_22_1991), col=r, border="black")
        lines(density(harm.cent_22_2001), col="blue")
        polygon(density(harm.cent_22_2001), col=b, border="black")
        lines(density(harm.cent_22_2017), col="green")
        polygon(density(harm.cent_22_2017), col=g, border="black")
        
        plot(density(harm.cent_12_1991), xlab="", ylab="# of nodes", main="Country 1-2: External links", cex=0.8, 
             col="red", xlim = c(0,1),ylim=c(0,40))
         legend(0, 53, legend=c("1991","2001", "2017"),
                col=c(r,b,g), lty=NA, pch=21, cex=1, bty = "n", lwd=20, xpd = "NA") 
        polygon(density(harm.cent_12_1991), col=r, border="black")
        lines(density(harm.cent_12_2001), col="blue")
        polygon(density(harm.cent_12_2001), col=b, border="black")
        lines(density(harm.cent_12_2017), col="green")
        polygon(density(harm.cent_12_2017), col=g, border="black")
                
        mtext("Harmonic Centrality [normalized / unweighted]", 
              outer=TRUE, cex=1, font=2)
        
  # Compute Harmonic Centrality [normalized / weighted] ---------------------------------------------
    #This block of code computes the normalized harmoic centrality per network per year    
        
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
          assign(h, harmonic_centrality(x, mode = "all", weights =1+E(x)$w2)/189)
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
        remove(x,y,n,h,d,d1,hcb_w)
        
        #This code plots the compiled normalized and weighted harmonic centrality per year 
        par(old.par)
        plot(harm_cent22, type = "l",col = "red", xlab = "Year", ylab = "Harmonic Centrality (Mean)", 
             main = "Harmonic centrality [normalized / weighted]", ylim=c(0.4,1),lwd=2) 
        lines(harm_cent11, type = "l", col = "darkgreen",lwd=2)
        lines(harm_cent12, type = "l", col = "blue",lwd=2)
        legend(1990, 1, legend=c("C1:Internal","C2:Internal", "C1-C2:External"),
               col=c("darkgreen","red", "blue"), lty=NA, cex=0.8, bty="n", lwd=10, pch=21) 
        
        #This code create density plots for the harmonic centraliy in 1991,2001, and 2017            
        b <- adjustcolor("blue", alpha.f = 0.2)
        r <- adjustcolor("red", alpha.f = 0.2)
        g <- adjustcolor("green", alpha.f = 0.2)
        par(mfrow=c(3,1), oma=c(0,0,2,0))
        plot(density(harm.cent_11_1991), xlab="", ylab="# nodes", main="Country 1: Internal links", cex=0.8, 
             col="red", xlim = c(0,1),ylim=c(0,45)) 
        legend(0, 58, legend=c("1991","2001", "2017"),
               col=c(r,b,g), lty=NA, pch=21, cex=1, bty = "n", lwd=20, xpd = "NA") 
        polygon(density(harm.cent_11_1991), col=r, border="black")
        lines(density(harm.cent_11_2001), col="blue")
        polygon(density(harm.cent_11_2001), col=b, border="black")
        lines(density(harm.cent_11_2017), col="green")
        polygon(density(harm.cent_11_2017), col=g, border="black")
        
        plot(density(harm.cent_22_1991), xlab="", ylab="# nodes", main="Country 2: Internal links", cex=0.8, 
             col="red", xlim = c(0,1),ylim=c(0,45))
        legend(0, 58, legend=c("1991","2001", "2017"),
               col=c(r,b,g), lty=NA, pch=21, cex=1, bty = "n", lwd=20, xpd = "NA") 
        polygon(density(harm.cent_22_1991), col=r, border="black")
        lines(density(harm.cent_22_2001), col="blue")
        polygon(density(harm.cent_22_2001), col=b, border="black")
        lines(density(harm.cent_22_2017), col="green")
        polygon(density(harm.cent_22_2017), col=g, border="black")

        plot(density(harm.cent_12_1991), xlab="", ylab="", main="Country 1-2: External links", cex=0.8, 
             col="red", xlim = c(0,1),ylim=c(0,45))
        legend(0, 58, legend=c("1991","2001", "2017"),
               col=c(r,b,g), lty=NA, pch=21, cex=1, bty = "n", lwd=20, xpd = "NA") 
        polygon(density(harm.cent_12_1991), col=r, border="black")
        lines(density(harm.cent_12_2001), col="blue")
        polygon(density(harm.cent_12_2001), col=b, border="black")
        lines(density(harm.cent_12_2017), col="green")
        polygon(density(harm.cent_12_2017), col=g, border="black")
        
        mtext("Harmonic Centrality [normalized / weighted]", 
              outer=TRUE, cex=1, font=2)
        
  # Fixed effects regressions ------------------------------------------------
    #Reg Model 1: Intensity ~ border (with year and region fixed effects) 
        #On partial data
        pyl <- plm(log(1+Intensity) ~ border + r_1 + r_2 + r_3 +
                     r_4	+ r_5	+ r_6	+  r_7	+ r_8	+ r_9	+ r_10	+ r_11	+ r_12	+ r_13	+ r_14	+ r_15	+ r_16	+
                     r_17	+ r_18	+ r_19	+ r_20	+ r_21	+ r_22	+ r_23	+ r_24	+ r_25	+ r_26	+ r_27	+ r_28	+ r_29	+
                     r_30	+ r_31	+ r_32	+ r_33	+ r_34	+ r_35	+ r_36	+ r_37	+ r_38	+ r_39	+ r_40	+ r_41	+ r_42	+
                     r_43	+ r_44	+ r_45	+ r_46	+ r_47	+ r_48	+ r_49	+ r_50	+ r_51	+ r_52	+ r_53	+ r_54	+ r_55	+
                     r_56	+ r_57	+ r_58	+ r_59	+ r_60	+ r_61	+ r_62	+ r_63	+ r_64	+ r_65	+ r_66	+ r_67	+ r_68	+
                     r_69	+ r_70	+ r_71	+ r_72	+ r_73	+ r_74	+ r_75	+ r_76	+ r_77	+ r_78	+ r_79	+ r_80	+ r_81	+
                     r_82	+ r_83	+ r_84	+ r_85	+ r_86	+ r_87	+ r_88	+ r_89	+ r_90	+ r_91	+ r_92	+ r_93	+ r_94	+
                     r_95	+ r_96	+ r_97	+ r_98	+ r_99	+ r_100 + r_101	+ r_102	+ r_103	+ r_104	+ r_105	+ 
                     r_106	+
                     r_107	+ r_108	+ r_109	+ r_110	+ r_111	+ r_112	+ r_113	+ r_114	+ r_115	+ r_116	+
                     r_117	+ r_118	+ r_119	+ r_120	+ r_121	+  r_122	+ r_123	+  r_124	+ r_125	+  r_126	+
                     r_127	+  r_128	+ r_129	+ r_130	+  r_131	+ r_132	+ r_133	+  r_134	+ r_135	+  r_136	+
                     r_137	+  r_138	+ r_139	+ r_140	+  r_141	+ r_142	+ r_143	+  r_144	+ r_145	+  r_146	+
                     r_147	+  r_148	+ r_149	+ r_150	+  r_151	+ r_152	+ r_153	+  r_154	+ r_155	+  r_156	+
                     r_157	+  r_158	+ r_159	+ r_160	+  r_161	+ r_162	+ r_163	+  r_164	+ r_165	+  r_166	+
                     r_167	+  r_168	+ r_169	+ r_170	+  r_171	+ r_172	+ r_173	+ r_174	+ r_175	+ r_176	+
                     r_177	+  r_178	+ r_179	+ r_180	+  r_181	+ r_182	+ r_183	+ r_184	+  r_185	+ r_186	+
                     r_187	+  r_188	+ r_189,
                  data = par_set, 
                  index = c("Year"),
                  effect = "time",
                  model = "within")
        summary(pyl)
 
        #On full data 
         fyl <- plm(log(1+Intensity) ~ border+ r_1 + r_2 + r_3 +
                      r_4	+ r_5	+ r_6	+  r_7	+ r_8	+ r_9	+ r_10	+ r_11	+ r_12	+ r_13	+ r_14	+ r_15	+ r_16	+
                      r_17	+ r_18	+ r_19	+ r_20	+ r_21	+ r_22	+ r_23	+ r_24	+ r_25	+ r_26	+ r_27	+ r_28	+ r_29	+
                      r_30	+ r_31	+ r_32	+ r_33	+ r_34	+ r_35	+ r_36	+ r_37	+ r_38	+ r_39	+ r_40	+ r_41	+ r_42	+
                      r_43	+ r_44	+ r_45	+ r_46	+ r_47	+ r_48	+ r_49	+ r_50	+ r_51	+ r_52	+ r_53	+ r_54	+ r_55	+
                      r_56	+ r_57	+ r_58	+ r_59	+ r_60	+ r_61	+ r_62	+ r_63	+ r_64	+ r_65	+ r_66	+ r_67	+ r_68	+
                      r_69	+ r_70	+ r_71	+ r_72	+ r_73	+ r_74	+ r_75	+ r_76	+ r_77	+ r_78	+ r_79	+ r_80	+ r_81	+
                      r_82	+ r_83	+ r_84	+ r_85	+ r_86	+ r_87	+ r_88	+ r_89	+ r_90	+ r_91	+ r_92	+ r_93	+ r_94	+
                      r_95	+ r_96	+ r_97	+ r_98	+ r_99	+ r_100 + r_101	+ r_102	+ r_103	+ r_104	+ r_105	+ 
                      r_106	+
                      r_107	+ r_108	+ r_109	+ r_110	+ r_111	+ r_112	+ r_113	+ r_114	+ r_115	+ r_116	+
                      r_117	+ r_118	+ r_119	+ r_120	+ r_121	+  r_122	+ r_123	+  r_124	+ r_125	+  r_126	+
                      r_127	+  r_128	+ r_129	+ r_130	+  r_131	+ r_132	+ r_133	+  r_134	+ r_135	+  r_136	+
                      r_137	+  r_138	+ r_139	+ r_140	+  r_141	+ r_142	+ r_143	+  r_144	+ r_145	+  r_146	+
                      r_147	+  r_148	+ r_149	+ r_150	+  r_151	+ r_152	+ r_153	+  r_154	+ r_155	+  r_156	+
                      r_157	+  r_158	+ r_159	+ r_160	+  r_161	+ r_162	+ r_163	+  r_164	+ r_165	+  r_166	+
                      r_167	+  r_168	+ r_169	+ r_170	+  r_171	+ r_172	+ r_173	+ r_174	+ r_175	+ r_176	+
                      r_177	+  r_178	+ r_179	+ r_180	+  r_181	+ r_182	+ r_183	+ r_184	+  r_185	+ r_186	+
                      r_187	+  r_188	+ r_189, 
                   data = full_set,
                   index = c("Year"),
                   effect = "time",
                   model = "within")
         summary(fyl)

    #Reg Model II: Degree centrality ~ border 
        dc_py <- plm(dcb ~ dc, 
                   data = dc,
                   index = c("ID","Year"),
                   model = "within", 
                   effect = "twoways")
      summary(dc_py)
       