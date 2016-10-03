# Masters-Simina
Master's thesis code

############################################################################################################################
####################### FLOWERING ############################################################################################
##############################################################################################################################
##############################################################################################################################

#### convertin raw data to time series for fourier analysis

setwd("R:/rsrch/cb751/lab/Simina")

#######################################################################################################################

calc.four_norm <- function(F_df = Fourier_df) {
  F_df$spec_norm_a <- NA
  for (ID in unique(Fourier_df$ID)) {
    F_df1 <- F_df[F_df$ID == ID,]
    a <- 0
    for (i in 2: NROW(F_df1)){
      a <- a + ((F_df1$freq[i] - F_df1$freq[i-1]) * 
                  ((F_df1$spec[i] + F_df1$spec[i-1])/2))
    }
    F_df$spec_norm_a[F_df$ID == ID] <- F_df1$spec / a
  }
  return(F_df)
}


#Upload the libraries you need for the results table and plots 
library(plyr)
library(RColorBrewer)
library(ggplot2)
library(reshape)
library(grid)
library(circular)

########################### Prepare time series for analysis #################################################################

#create a new dataframe that will act as a universal month x year identifier

timeline1 <- data.frame(year=rep(1900:2015,each=12),
                        month=rep(c(1,2,3,4,5,6,7,8,9,10,11,12), NROW(1900:2015)))
timeline1$jmonth <- 1:NROW(timeline1)
timeline1$date <- strptime(paste("01", timeline1$month, timeline1$year, sep = "-"), format = "%d-%m-%Y")


#prepare the time series function for the analyses 
ts_fun <- function(d_l){
  ts(data = d_l[[1]], start = c(timeline1$year[timeline1$jmonth == min(d_l[[2]])],
                                + timeline1$month[timeline1$jmonth == min(d_l[[2]])]), freq = 12)
}

# Specify and upload the files you want to insert in your loop 

sites <- dir(paste0(getwd(), "/Data/Flowers/Sites/"), pattern = ".csv") #find all csv files in flowers folder; These functions produce a character vector of the names of files or directories in the named directory.

# Create a list that will contain all sites 

big_list=list()

for(ss in 1:length(sites)){
  
  data_ls <- list() #create a list to fill in later with all the species for one site 
  
  my_file=read.csv(paste0(getwd(), "/Data/Flowers/Sites/", sites[[ss]]))
  
  for(sp in levels(my_file$Species)){
    dat.sp <- my_file[my_file$Species == sp,] # Distinguish between all Species within all sites
    dat.sp <- dat.sp[order(dat.sp$Tree.ID, dat.sp$jmonth),]
    
    tree_ls <- list()
    
    for (i in unique(dat.sp$Tree.ID)){
      mon1 <-  dat.sp$jmonth[dat.sp$Tree.ID == i]  ### using numeric indexing we extract jmonth; for every Tree iD; raw data  
      mon <- min(mon1, na.rm=T):max(mon1, na.rm=T) #min till max: completes the numbers with the missing jmonth 
      dat <- rep(NA, NROW(mon)) # repllicates NA in the length of mon (95)
      dat1 <- dat.sp$value[dat.sp$Tree.ID == i]   # Finds flowering score for every tree ID
      used.mons <- mon %in% mon1 #what is in mon1 and not in mon 
      
      if(sum(used.mons) != NROW(dat1)){ #if used.mons is not equal to dat1; length wise 
        dups <- duplicated(mon1) # no duplicates in  mon1, the initial jmonth in our data
        mon1 <- mon1[!dups]# ! logical operator (negative). 
        dat1 <- dat1[!dups]
        used.mons <- mon %in% mon1
      }
      
      dat[used.mons] <- dat1 
      dat.nas <- (is.na(dat)) * unlist(lapply(rle(is.na(dat))$lengths, seq_len))  ## finds conequtive NAs in the value column 
      #print(sum(!is.na(dat)))
      
      if(sum(!is.na(dat))>2) {
        dat[is.na(dat)] <- approx(x= mon[!is.na(dat)], y = dat[!is.na(dat)], xout = mon[is.na(dat)] )$y  ##interpolates all NAs - linear
        
        
        ## If gaps are too long, cut dataset
        if(any(dat.nas>3)){ # if gaps are longer than 3 months then cut 
          data_list <- list()
          cuts <- which(dat.nas>3) # cuts where there is more than 3 NAs in consecutive 3 months 
          cuts <- c(cuts, NROW(dat.nas) + 1)  ## need to add to the cuts the last value in the vector - but plus one because the code below assumes that the cut vector identifies the first na, not the last non-na and the last value is required.
          bot <- 1
          for (j in 1:NROW(cuts)){
            dat_l <- dat[bot:(cuts[j] - 1)]   ### cuts original data using the positions of the first value that is not na and last one before a cut(cuts identify the FIRST na values, not the last non na)
            mon_l <- mon[bot:(cuts[j] - 1 )]
            dat_l <- dat_l[!is.na(dat_l)]   ### removes residual nas left at start of sequence if only 3 - approx does not extrapolate.
            mon_l <- mon_l[!is.na(dat_l)]
            if(NROW(dat_l) > 3){
              data_list[[j]] <- list(dat_l, mon_l)
            }
            bot <- cuts[j] + 1
            
          }
          
          data_list <- data_list[unlist(lapply(data_list, NROW))>0]
          tree.ts <- lapply(data_list, ts_fun) #apply the ts_function on data_list 
          tree_ls <- c(tree_ls, tree.ts)
          
        } else { 
          mon <- mon[!is.na(dat)] # mon takes all the values that are not NAs 
          dat <- dat[!is.na(dat)] # dat takes all the values that are not NAs 
          data_list <- list(dat, mon) # creates a list where 
          
          
          tree.ts <- ts_fun(data_list)
          tree_ls <- c(tree_ls, list(tree.ts))
        }
      }
    }
    tree_ls <- tree_ls[lapply(tree_ls, NROW)>60]  ### Strips any fragments of time series less than 60 months long (5 years) 
    tree_ls <- tree_ls[lapply(tree_ls, sum)>0]  ### Strips any trees that never flower
    # combine each species
    data_ls[[sp]]<-tree_ls
  }
  data_ls <- data_ls[lapply(data_ls, NROW)>9]### strips any species with fewer than 10 individuals that flowered
  big_list[[ss]]=data_ls
}

names(big_list) <- unlist(lapply(strsplit(sites, "_"), function(x) x[1]))

# Remove all the sites that have no data available to analyse (data too short, length: 0) 
# Sites that did not pass the above filters were removed: Sites that will not be taken into consideration for this analysis 

##########################################################################################################################

### Now we have a list of all tree species, each with a list of all individuals for each site in data_ls
### The following code works on all the individuals of a species, so we must loop through species.

FILEPATH="/R:/rsrch/cb751/lab/Simina/Data/Plot/"

for(ss in 1:length(sites)){
  
  # for(ss in 1:2){
  
  pdf(paste0(sites[ss], ".pdf"), height = 6, width = 8)
  
  for (sp in names(big_list[[ss]])){   ### extract tags from big_list for species
    tree_ls <- big_list[[ss]][[sp]]
    
    if (NROW(tree_ls)>9) {
      
      individuals<-length(tree_ls)
      IDs<-1:individuals
      
      data_df<-data.frame()
      # selecting complete dates
      
      #d=data.frame(Date=seq(from=strptime(paste(2004,9,01),"%Y %m %d"),to=strptime(paste(2012,7,01),"%Y %m %d"), by="month"))          
      
      for(i in 1:individuals){
        d=data.frame(Date=seq(from=strptime(paste(start(tree_ls[[i]])[1],
                                                  start(tree_ls[[i]])[2],01),"%Y %m %d"),
                              to=strptime(paste(end(tree_ls[[i]])[1],
                                                end(tree_ls[[i]])[2],01),"%Y %m %d"), by="month"))          
        d$Year <- c(format(d$Date, format="%Y")) # gets year
        d$Month <- c(format(d$Date, format="%m")) # gets month
        d$TS <- tree_ls[[i]][1:length(tree_ls[[i]])]
        d$ID <-  IDs[i]
        data_df <- rbind(data_df,d)
      }
      #################################### To make circular box plot - FIGURE 3 ######################################
      
      #Summarise phenology scores across sample for each month and year combination
      #(to give proportions for boxplots)
      
      data_df_sample_summary<-ddply(data_df,.(Year,Month),summarize,
                                    meanPhenoScore=mean(TS),
                                    propPhenophase=length(TS[TS>0])/length(IDs))
      
      #Summarise phenology scores across sample and years for each month (to give mean
      #score for boxplot fill colour)
      meanPhenoScore_fill<-ddply(data_df_sample_summary,.(Month),summarize,
                                 meanScore_month=mean(meanPhenoScore))
      
      meanPhenoScore_fill$levels<-cut(meanPhenoScore_fill$meanScore_month,
                                      breaks=c(0,0.5,1,1.5,2,2.5,3,3.5,4),
                                      labels=c(0.5,1,1.5,2,2.5,3,3.5,4))
      
      data_df_sample_summary$Month<-mapvalues(data_df_sample_summary$Month,from=c("01","02","03","04","05","06","07","08","09","10","11","12"), to=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))
      data_df_sample_summary$Month<-factor(data_df_sample_summary$Month,ordered=T,levels=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))
      
      # Set fill colours according to mean score for each month
      myColors <- data.frame(Colour=(brewer.pal(9,"OrRd")),levels=c(NA,0.5,1,1.5,2,2.5,3,3.5,4))
      meanPhenoScore_fill<-merge(meanPhenoScore_fill,myColors,by="levels",all.x=T)
      colours_fill<-as.character(meanPhenoScore_fill[order(meanPhenoScore_fill$Month),]$Colour)
      
      #Write circular boxplot graphic
      Phenology_circular_boxplot<-ggplot(data_df_sample_summary,aes(x=factor(Month),y=propPhenophase))+
        geom_boxplot(aes(fill=factor(Month)),outlier.colour="lightgrey")+
        scale_colour_manual(values=colours_fill)+  
        scale_fill_manual(values=colours_fill)+
        coord_polar(start=0)+
        theme_minimal(base_size=18)+
        ggtitle(paste0("Smoothed periodogramsby individual ", sp))+
        theme(legend.position="none",axis.title.x=element_blank(),axis.ticks=element_blank())+
        ylab("")
        
      write.csv(data_df_sample_summary,paste0("R:/rsrch/cb751/lab/Simina/test/",sites[ss], ".csv"),row.names = F) # summary statistics: all species per site 
      
      #Write colour scale bar as legend for boxplot fill
      my.colors = colorRampPalette(brewer.pal(9,"OrRd"))
      df <- melt(outer(seq(0,4,len=100), 1), varnames = c("X1", "X2"))
      df$fill<- as.factor(cut(df$value, breaks=100,labels=F))
      
      colourbar<-ggplot(df, aes(X2, value))+
        geom_tile(aes(fill=fill))+
        scale_fill_manual(values=my.colors(100),guide=F)+ 
        labs(x = "",  y = "") + 
        theme_minimal(base_size=18)+
        theme(axis.ticks = element_blank(), axis.text.x = element_blank(),axis.text.y = element_blank(),panel.grid.minor=element_blank(), panel.grid.major=element_blank())+
        coord_flip()
      
      
      ##Plot results: #timeseries+periodogram plot
      vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(5,5,heights=unit((c(4,3,4,3,3)),"null"))))
      print(Phenology_circular_boxplot, vp = vplayout(1:4,1:5))
      print(colourbar, vp = vplayout(5,2:4))
      
      ############################# Smoothed periodograms and duration plots ##############################
      
      ##Functions for periodogram analysis 
      spec_fun<-function(x) spectrum(x,spans=spans_comb[[which.min(abs(months_comb-length(x)))]] ,plot=F,demean=T,detrend=T) #spectrum function for normal smoother periodogram
      
      ##Spans
      #The following lists give appropriate successive spans for the corresponding
      #length of data wto give a smoothed periodogram with bandwidth in magnitude of
      #0.1.
      months_comb<-c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30)*12
      spans_comb<-list(3,3,3,3,c(3,3),c(3,5),c(3,5),c(5,5),c(5,5),c(5,7),c(5,7),c(7,7),c(7,9),c(7,9),c(7,9))
      
      ##Run fourier analysis using function "spectrum" on each list object (individual
      ##time series)
      #  Fourier_ls<-llply(data_ls,function(x) spectrum(x,spans=spans_comb[[which.min(abs(months_comb-length(x)))]],demean=T,detrend=T,plot=F))
      Fourier_ls<-llply(tree_ls,function(x) spectrum(x,spans=spans_comb[[which.min(abs(months_comb-length(x)))]],demean=T,detrend=T,plot=F))
      
      #Merge spectrum outputs into a dataframe
      Fourier_df<-data.frame()
      for (i in 1:individuals){
        d<-data.frame(freq=Fourier_ls[[i]]$freq/12)
        d$spec=Fourier_ls[[i]]$spec
        d$spec_norm=Fourier_ls[[i]]$spec/sum(Fourier_ls[[i]]$spec)
        d$ID=as.factor(IDs[i])
        Fourier_df<-rbind(Fourier_df,d)
      }
      Fourier_df <- calc.four_norm()
      
      #Summarise dominant peaks
      Fourier_df_sample_summary<-ddply(Fourier_df, .(ID), summarize,
                                       maxFreq=freq[which.max(spec)],
                                       maxSpec_norm=max(spec_norm_a))
      
      Dominant_wavelength<-1/median(Fourier_df_sample_summary$maxFreq)
      
      ##Plot periodograms for visual inspection
      #Periogorams on one graph
      Periodograms<-ggplot(data=Fourier_df,aes(x=freq,y=spec_norm_a,group=ID,colour=ID)) +
        geom_line()+
        theme_minimal(base_size=18)+
        ylab("standardised spectrum")+
        xlab("frequency (cycles per month")+
        ggtitle(paste0("Smoothed periodogramsby individual ", sp))+
        theme(legend.position="none")
      
      #Optional to plot periodograms as facet
      Periodograms_facet<-ggplot(data=Fourier_df,aes(x=freq/12,y=spec_norm_a,colour=ID)) +
        geom_line()+
        theme_minimal()+
        xlab("frequency (cycles per month)")+
        ylab("standardised spectrum")+
        theme(legend.position="none")+
        facet_wrap(~ID,ncol=2)
      #ggtitle("Smoothed periodograms of phenology cycles by individual")
      
      #Plot individual time series durations
      Timeseries_duration<-ggplot(data=data_df,aes(y=Date,x=ID))+
        geom_point(aes(colour=ID))+
        coord_flip()+
        theme_minimal(base_size=18)+
        xlab("")+
        ylab("year")+
        theme(legend.position="none",axis.title.y=element_blank(),axis.ticks = element_blank(), axis.text.y = element_blank())
      #ggtitle("Duration of each timeseries")
      
      
      ##Plot results: #timeseries+periodogram plot
      vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(1,3,heights=unit((c(4)),"null"))))
      print(Periodograms, vp = vplayout(1,2:3))
      print(Timeseries_duration, vp = vplayout(1,1))
    }
  }
  dev.off()
}


#####################################################################################
#### Periodogram analysis with confidence intervals to find frequency and phase #####
#####################################################################################

##Find frequency of dominant cycles and test against null hypothesis

#Functions
spec_fun<-function(x) spectrum(x,spans=spans_comb[[which.min(abs(months_comb-length(x)))]],
                               plot=F,demean=T,detrend=T) #spectrum function for normal smoother periodogram
spec_null_fun<-function(x) spectrum(x,spans=spans_super_smooth_comb[[which.min(abs(months_comb-length(x)))]] ,plot=F,demean=T,detrend=T) #spectrum function for null hypothesis spectrum

spec_fun_a <-function(x) {
  spec <- spectrum(x,spans=spans_comb[[which.min(abs(months_comb-length(x)))]],
                   plot=F,demean=T,detrend=T) #spectrum function for normal smoother periodogram
  a <- 0
  for (i in 2: NROW(spec$spec)){
    a <- a + ((spec$freq[i]/12 - spec$freq[i-1]/12) * 
                ((spec$spec[i] + spec$spec[i-1])/2))
  }
  a
}

confidence_fun<-function(p, ind_ls = big_list[[3]][[1]]){
  ts<-ind_ls[[p]]
  d<-data.frame(ID=p)
  d$length<-length(ts)
  d$freq_max<-(spec_fun(ts)$freq[which.max(spec_fun(ts)$spec)])/12
  d$spec_max<-max(spec_fun(ts)$spec)
  d$spec_max_norm<-d$spec_max/spec_fun_a(ts) #normalise spectrum so can compare between individuals
  d$df<-spec_fun(ts)$df
  d$bw_smooth<-spec_fun(ts)$bandwidth
  d$bw_super_smooth<-spec_null_fun(ts)$bandwidth
  d$lower_ci<-(d$df*d$spec_max)/(qchisq(c(0.975),d$df))
  d$null_spec_max<-spec_null_fun(ts)$spec[which(abs((spec_null_fun(ts)$freq)/12-d$freq_max)==min(abs((spec_null_fun(ts)$freq)/12-d$freq_max)))]                     
  d$sig_null<-d$lower_ci>d$null_spec_max #If the null spectrum is lower than the lower confidence interval, dominant peak can be considered "signficantly different" to null hypothesis
  data.frame(d)
}

##Spans
#The following lists give appropriate successive spans for the corresponding
#length of data to give a smoothed periodogram with bandwidth in magnitude of
#0.1 (spans_comb) and a null continuum smoother periodogram with bandwidth in
#magnitude of 1 (spans_super_smooth_comb)
months_comb<-c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30)*12
spans_comb<-list(3,3,3,3,c(3,3),c(3,5),c(3,5),c(5,5),c(5,5),c(5,7),c(5,7),c(7,7),c(7,9),c(7,9),c(7,9))
spans_super_smooth_comb<-list(c(5,7),c(11,11),c(15,17),c(19,21),c(25,27),c(29,31),c(35,37),c(39,41),c(43,45),c(47,49),c(55,57),c(59,61),c(65,67),c(73,75),c(75,79))#spans_smooth_fun<-function(x) ifelse(odd(round(sqrt(sqrt(x))))==TRUE,round(sqrt(sqrt(x))),round(sqrt(sqrt(x)))+1)

# Run loop to apply spectrum functions to each individual time series - currently using only tree species 1 (data_ls[[1]])
confidence_results<-ldply(1:length(big_list[[3]][[1]]),.fun=confidence_fun, ind_ls = big_list[[3]][[1]])  ### CMB

#Summary of frequency of dominant cycle (wavelength = 1/freq)
summary(confidence_results$freq_max)
Dominant_wavelength<-round(1/median(confidence_results[confidence_results$sig_null==T,]$freq_max))# you can set it as 12 months instead of taking the dominant
#Dominant_wavelength<-6 # this was also run to find out the accurate phase for the 6 months frequencies
## Find phase of significant dominant cycles and assess synchroncity between individuals

#Functions
spec_cross_fun<-function(x)spec.pgram(x,spans=spans_comb[[which.min(abs(months_comb-length(x)))]] ,demean=T,detrend=T,plot=F)

phase_fun<-function(p, ind_ls = big_list[[3]][[1]]){
  ts_cross<-ts.intersect(simulated_phase_ts,ind_ls[[p]]) #intersect them 
  d<-data.frame(ID=IDs[p])
  d$freq<-1/12*spec_cross_fun(ts_cross)$freq[which.max(spec_cross_fun(ts_cross)$spec[,1])]
  d$phase<-spec_cross_fun(ts_cross)$phase[which.max(spec_cross_fun(ts_cross)$spec[,1])]
  d$coh<-spec_cross_fun(ts_cross)$coh[which.max(spec_cross_fun(ts_cross)$spec[,1])]
  d
}

#Simulate a "template" time series using dominant wavelength identified above
#with phase 0 to act as guide in co-Fourier analysis
simulated_phase<-list(t=1:360,sd=0,w=2*pi*(1/Dominant_wavelength),R=4,p=0,ts=NA)
simulated_phase$ts<-simulated_phase$R*cos(simulated_phase$w*simulated_phase$t+simulated_phase$p)

simulated_phase_ts<-ts(simulated_phase$ts,start=c(1986,1),freq=12)

# Run loop to apply co-fourier spectrum function to each individual time series
# in combination with simulated guide time series to identify phase
phase_results<-ldply(1:length(big_list[[3]][[1]]),.fun=phase_fun, ind_ls = big_list[[3]][[1]])   ### only for first species in list at the moment CMB
results<-cbind(confidence_results,phase_results[,c(1,3,4)])

results <- data.frame(ID =NA, length =NA, freq_max =NA,   spec_max=NA, spec_max_norm=NA,
                      df=NA,   bw_smooth =NA, bw_super_smooth=NA,  lower_ci=NA,  
                      null_spec_max=NA,  sig_null=NA, freq = NA,   phase=NA, coh=NA, species = NA, site =NA)

results <- results[!is.na(results[,1]),]

## Now add a column for the species and another for the site.
#big_results=list()

for(ss in 1:length(sites)){
  # setwd("R:/rsrch/cb751/lab/Simina")
  # my_file=read.csv(paste0(getwd(), "/Data/Flowers/", sites[[ss]]))
  #my_file=read.csv(sites[ss])
  for (sp in 1:NROW(big_list[[ss]])){
    confidence_results<-ldply(1:length(big_list[[ss]][[sp]]),.fun=confidence_fun, ind_ls = big_list[[ss]][[sp]])  ### CMB
    phase_results<-ldply(1:length(big_list[[ss]][[sp]]),.fun=phase_fun, ind_ls = big_list[[ss]][[sp]])   ### only for first species in list at the moment CMB
    results.t <-cbind(confidence_results,phase_results[,c("freq","phase","coh")])
    results.t  <- results.t[!is.nan(results.t$spec_max_norm),]
    results.t$species <- names(big_list[[ss]])[sp]
    results.t$site <- names(big_list)[ss]
    results <- rbind(results, results.t)
    #  results$site <- ## from site loop...
  }
  #  big_results[[ss]]=results
}

results$cyclicity<-ifelse(1/results$freq_max>results$length/2,"No cyclicity","Cyclicity")
results$cycle_category<-ifelse(1/results$freq_max>results$length/2,"No cyclicity",
                               ifelse(1/results$freq_max>13,"Supra-annual",
                                      ifelse(1/results$freq_max>11,"Annual","Sub-annual")))

a<-(2*pi/12)
results$phase_month <-ifelse(results$phase > 0, results$phase_month <- (results$phase)/a,((results$phase)+(2*pi)) / a)
results$freq_in_months<-1/results$freq_max
results<-results[!results$cycle_category== "No cyclicity", ]#remove the ones that do not have any cyclycity - approx 600 were removed

### results file contains individual tree flowering frequency, fidelity for a particular cycle and phase. 
### accurate phase is the one for the 12 months trees only


