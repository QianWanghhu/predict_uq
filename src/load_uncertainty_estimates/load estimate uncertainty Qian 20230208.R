
## Estimating uncertainty in constituent annual load estimation based on interpolation and beale ratio methods
## Baihua Fu  ## Edited by Qian Wang 
## 9/12/2021  ## 08/02/2023
## This edited version changes the way of calculating load uncertainty assuming that there is a linear relationship between the uncertainty of L_k-1 and L_k

# Set work space----
### set word directory
setwd("E:\\BaiduSyncdisk\\Projects\\predict_uq\\src\\load_uncertainty_estimates")
### load libraries
library(dplyr)
library(ggplot2)
library(zoo)
library(gridExtra)
library(lubridate)

# Generate functions ---------------------------------------------------------------

### Calculate hydrology year function
as.year <- function(x) as.integer(as.yearmon(x) - 6/12)

### Linear interpolation function 
# Test the functions below
conclinear <- function(flowdf, concdf, minconc){
  flowzoo <- read.zoo(flowdf, header = TRUE)
  conczoo <- read.zoo(concdf, header = TRUE, aggregate = mean) #take average for same time step
  
  # merge flow and concentration data
  flowconczoo <- merge(flowzoo, conczoo)
  names(flowconczoo) <- c("flow", "conc")
  
  # assign tiedowns: begin with half of min value of concentration, and end of the min value.
  flowconczoo[1,"conc"] <- minconc/2
  flowconczoo$conc[nrow(flowconczoo)] <- minconc
  
  # interpolate concentration
  flowconczoo$conc <- na.approx(flowconczoo$conc, rule = 1)
  return(flowconczoo)
  
}


loadlinear <- function(flowconczoo, seconds){
  
  # hourly load, mg/L * m3/s = g/s
  flowconczoo$load <- flowconczoo$conc*flowconczoo$flow
  
 # aggregate to annual mean instant flow, concentration and load
  yearly <- aggregate(flowconczoo, as.year(time(flowconczoo)), mean)
  yearly <- as.data.frame(yearly)
  yearly <- yearly[-1,]  # removed the first first row (10 data points, i.e. 10 hours) due to unresolved issue on timezone
  
  # calculate annual total load via: multiply mean instant load by number of seconds in that year, then convert unit to t/yr
  yearly$load <- round(yearly$load*seconds * 0.000001, 3) #g/s * s/yr = g/yr. 10^-6g/yr = t/yr
  return(yearly)
}


loadunc <- function(flowconczoo, concunc, flowunc, concdf, minconc, covar=1){
  conczoo <- read.zoo(concdf, header = TRUE, aggregate = mean)
  # Add tiedowns to the observed concentration data so as to enable interpolation
  concdfnew <- data.frame(Datetime = index(conczoo), conc = as.data.frame(conczoo)$conc)
  concdfnew <- rbind(concdfnew, data.frame(Datetime=flowunc$Datetime[1], conc=minconc/2))
  concdfnew <- rbind(concdfnew, data.frame(Datetime=flowunc$Datetime[length(flowunc$Datetime)], conc=minconc))
  concdfnew <- concdfnew[order(concdfnew$Datetime),]
  
  df_temp <- as.data.frame(flowconczoo)
  flowcondf <- data.frame(Datetime=flowunc$Datetime, flow=df_temp$flow, conc=df_temp$conc)
  flowcondf$loadraw <- flowcondf$flow * flowcondf$conc
  time_conc_sample = concdfnew$Datetime
  sample_conc_length <- length(time_conc_sample)
  # assign initial values used for calculation
  flowcondf$slope = 1
  concunc$concuncsq <- concunc$concuncbound**2
  names(concdfnew) <- c('Datetime', 'conc')
  flowcondf$flowuncbound <- flowunc$flowuncbound # Assign values of flow uncertainty bounds
  flowcondf$concuncbound <- 0
  flowcondf$conc_delta_i <- 0
  flowcondf$conc_delta_i2 <- 0
  time_end <- time_conc_sample[sample_conc_length]
  flow_time_end_bool <- flowcondf$Datetime == time_end
  flowcondf[flowcondf$Datetime == time_conc_sample[1],"concuncbound"] <- max(concunc$concuncbound)
  flowcondf[flow_time_end_bool, "concuncbound"] <- max(concunc$concuncbound)
  flowcondf[flow_time_end_bool, "conc_delta_i"] <-
    concdfnew[concdfnew$Datetime == time_end, "conc"] * flowcondf[flow_time_end_bool, "concuncbound"]
  
  flowcondf[flow_time_end_bool, "conc_delta_i2"] <-
    concdfnew[concdfnew$Datetime == time_end, "conc"] * flowcondf[flow_time_end_bool, "concuncbound"]
  
  
  for (i in 1:(sample_conc_length - 1)){
    if (i> 1){
      
      flowcondf[flowcondf$Datetime == time_conc_sample[i], "concuncbound"] <- mean(concunc[concunc$Datetime==time_conc_sample[i], "concuncbound"])
    }
  }
  
  for (i in 1:(sample_conc_length - 1)){
    # If a data is obtained by interpolation, calculate the concentration uncertainty using the eq. in chapter 4.
    diff_time_sample <- as.integer(difftime(time_conc_sample[i + 1], time_conc_sample[i], units = 'hours'))
    
    if (diff_time_sample > 1){
      # browser()
      t1_bool <- flowcondf$Datetime >= time_conc_sample[i]
      t2_bool <- flowcondf$Datetime < time_conc_sample[i+1]
      # # Vectorize the calculation to speed up
      flowcondf[t1_bool & t2_bool, "slope"] <- c(seq(diff_time_sample, 1) / diff_time_sample)
      flowcondf[t1_bool & t2_bool, "conc_delta_i"] <-
        concdfnew[concdfnew$Datetime == time_conc_sample[i], "conc"] * flowcondf[flowcondf$Datetime == time_conc_sample[i], "concuncbound"]
      
      flowcondf[t1_bool & t2_bool, "conc_delta_i2"] <-
        concdfnew[concdfnew$Datetime == time_conc_sample[i+1], "conc"] * flowcondf[flowcondf$Datetime == time_conc_sample[i+1], "concuncbound"]
    }
    # END if
  }
  if (covar==1){
    flowcondf[, "concuncbound"] <- ((flowcondf[, "slope"]*flowcondf[, "conc_delta_i"])** 2 + ((1 - flowcondf[, "slope"])*flowcondf[, "conc_delta_i2"])** 2)
    flowcondf[, "concuncbound"] <- round(sqrt(flowcondf[, "concuncbound"]) / flowcondf[, "conc"], 3)
    # Note that flowcondf$concuncbound is not squared.
    # The code in the line below calculate the square of load uncertainty at each time step without co-variance.
    flowcondf$loadbase = flowcondf$loadraw
    flowcondf$loaduncbound = (flowcondf$concuncbound + flowcondf$flowuncbound) * flowcondf$loadbase
  } else if (covar==0){
    flowcondf[, "concuncbound"] <- ((flowcondf[, "slope"]*flowcondf[, "conc_delta_i"])** 2 + ((1 - flowcondf[, "slope"])*flowcondf[, "conc_delta_i2"])** 2) / (flowcondf[, "conc"]**2)
    # The code in the line below calculate the square of load uncertainty at each time step without covariance.
    # Note that flowcondf$concuncbound is already the squared values.
    flowcondf[, "concuncbound"] <- round(flowcondf[, "concuncbound"], 3)
    flowcondf$loadbase = flowcondf$loadraw**2
    flowcondf$loaduncbound = (flowcondf$concuncbound + flowcondf$flowuncbound ** 2) * flowcondf$loadbase
  }
  else{
    flowcondf[, "concuncbound"] <- ((flowcondf[, "slope"]*flowcondf[, "conc_delta_i"])** 2 + ((1 - flowcondf[, "slope"])*flowcondf[, "conc_delta_i2"])** 2) / (flowcondf[, "conc"]**2)
    # The code in the line below calculate the square of load uncertainty at each time step without covariance.
    # Note that flowcondf$concuncbound is already the squared values.
    flowcondf[, "concuncbound"] <- round(flowcondf[, "concuncbound"], 3)
    flowcondf$loadbase = flowcondf$loadraw
    flowcondf$loaduncbound = sqrt(flowcondf$concuncbound + flowcondf$flowuncbound ** 2) * flowcondf$loadbase
    
  }
  # End if
  return(flowcondf)
}



loaduncannual <- function(flowcondf){
  loadunczoo <- read.zoo(flowcondf, header = TRUE)
  
  # aggregate to annual mean instant flow, concentration and load 
  yearlyunc <- aggregate(loadunczoo, as.year(time(loadunczoo)), sum)
  yearlyunc <- as.data.frame(yearlyunc)
  yearlyunc <- yearlyunc[-1,] 
  yearlyunc$loaduncbound <- round(yearlyunc$loaduncbound / yearlyunc$loadbase, 3)

  return(yearlyunc)
}

# Use a function to call the linear interpolation related functions
calloadunc <- function(flowdf, concdf, minconc, seconds, flowunc, concunc, yearlyloadbase, covar, boundtype){
  flowconczoo <- conclinear(flowdf, concdf, minconc)
  # browser()
  flowcondf <- loadunc(flowconczoo, concunc, flowunc, concdf, minconc, covar)
  yearlyloadunc <- loaduncannual(flowcondf)
  if (boundtype == "upper"){
    yearlyloadunc$load <- yearlyloadbase$load * (1 + yearlyloadunc$loaduncbound)
  } else {
    yearlyloadunc$load <- yearlyloadbase$load * (1 - yearlyloadunc$loaduncbound)
    }
  
  return(yearlyloadunc)
}


### Random sampling functions for sample frequency analysis
# random sample n, repeat N times. concentration data
randomall <- function(concdf, n, N){
  names(concdf) <- c("Datetime", "conc")
  randomout <- bind_rows(replicate(N, concdf %>% sample_n(n), simplify=F), .id="Obs")
  return(randomout)
}

# random sample n, repeat N times for a given hydrology year. concentration data
randomyear <- function(concdf, nprop, N){
  names(concdf) <- c("Datetime", "conc")
  concdf$hyear <- as.year(concdf$Datetime)
  set.seed(200)
  randomout <- bind_rows(replicate(N, concdf %>%
                                     group_by(hyear) %>%
                                     mutate(yearsamples = length(conc)) %>%
                                     slice_sample(prop=nprop, weight_by=yearsamples),
                                   simplify = FALSE),
                         .id="Obs")
  return(randomout)
}


### Estimating combined uncertainty with choice of random sample methods, data, and load estimation methods.
## Qian to change the settings of arguments and inside loop
samplefreqfun <- function(flowdf, concdf, seconds,minconc, loadfun, biasfun, randomfun, 
                          randomseeds, randomseedprop, randomtimes, flowunc, concunc, covar, boundtype){
  # get random samples
  if (randomfun == "randomall"){
    concsamplefreq <- randomall(concdf, randomseeds, randomtimes)
  } else if (randomfun == "randomyear"){
    concsamplefreq <- randomyear(concdf, randomseedprop, randomtimes)
  }
  
  hyear <- as.year(flowdf$Datetime)
  
  conclist <- data.frame(matrix(0, randomtimes, length(unique(hyear))), row.names = seq(1, randomtimes))
  names(conclist) <- unique(hyear)
  
  if (loadfun == "linear"){
    for (i in (unique(concsamplefreq$Obs))){
      concsub <- subset(concsamplefreq, concsamplefreq$Obs == i)
      concsub <- as.data.frame(concsub[,c("Datetime", "conc")])
      # conclist[[i]] <- loadlinear(flowdf, concsub, seconds, minconc)
      # Calculate the new base load with re-sampled concentration
      # browser()
      flowconcbase_temp <- conclinear(flowdf, concsub, minconc)
      yearlyloadbase_temp <- loadlinear(flowconcbase_temp, seconds)
      load_temp <- calloadunc(flowdf, concsub, minconc, seconds, flowunc, concunc, yearlyloadbase_temp, covar, boundtype) # Use Qian's new settings
      conclist[i, ] <- load_temp$load
      # conclist[[i]]$hyear <- unique(hyear) ## TODO
    }
  } else if (loadfun == "beale"){
    for (i in (unique(concsamplefreq$Obs))){
      concsub <- subset(concsamplefreq, concsamplefreq$Obs == i)
      concsub <- as.data.frame(concsub[,c("Datetime", "conc")])
      conclist[[i]] <- loadbeale(flowdf, concsub, seconds, biasfun)
      conclist[[i]]$hyear <- unique(hyear) ## TODO
    }
  }
  
  # concsampleall <- bind_rows(conclist, .id = "samples")
  # browser()
  # concsampleminmax <- concsampleall %>%
  #   group_by(hyear, flow) %>%
  #   dplyr::summarise(across(.cols = c("load"),list (min = min, max = max), na.rm=TRUE))
  # 
  return(conclist)
  
}


# Import data -------------------------------------------------------------
### Import and format flow data
# import flow data
flow <- read.csv("126001A_flows.csv", skip = 1, header = F)

# tidy flow data to date, flow, qc format
colnames(flow) <- rep(c("Datetime", "flow", "QC","X"), 10)
flow <- rbind(flow[,c(1:3)], flow[,c(5:7)], flow[,c(9:11)], flow[,c(13:15)], flow[,c(17:19)], 
              flow[,c(21:23)], flow[,c(25:27)], flow[,c(29:31)], flow[,c(33:35)])

flow <- na.omit(flow) 

flow$Datetime <- strptime(flow$Datetime, format="%d/%m/%Y %H:%M", tz = "Australia/Brisbane")

# Qian's edit
flowunc <- flow
flowunc$flowuncbound <- 0
  
# create flow upper bounds and lower bounds data
flowqc <- flow
flow <- flow[,-3]

# create flow upper bounds and lower bounds data
flowuncup <- mutate(flowunc, flowuncbound = case_when(
  QC == 9 ~ 0.05,
  QC == 10 ~ 0.1,
  QC == 20 ~ 0.2,
  QC == 30 ~ 0.5,
  QC == 60 ~ 0.75,
  TRUE ~ flowuncbound))

flowunclow <- mutate(flowunc, flowuncbound = case_when(
  QC == 9 ~ 0.05,
  QC == 10 ~ 0.1,
  QC == 20 ~ 0.2,
  QC == 30 ~ 0.5,
  QC == 60 ~ 0.75,
  TRUE ~ flowuncbound))
# Form the flowunc* into two columns (Datatime, flowuncbound)
flowuncup <- flowuncup[,c(1, 4)]
flowunclow <- flowunclow[,c(1, 4)]
flowuncbest <- flowunc[,c(1, 4)]


### import and format concentration data
conc <- read.csv("126001A_din_concentrations_conditions.csv", header=TRUE, stringsAsFactors=FALSE) #TSS mg/L, NOx mg/L

# take only time, TSS, NOX, conditions
conc <- conc[,c(3:11)]
colnames(conc) <- c("Datetime", "NOx", "NH4", "DIN", "Flow_at_sample_time", "Equipment_contamination", "Sample_Temp", "Days_from_sample_to_analysis", "Days_to_filtering")

# round the concentration data to nearest hour
conc$Datetime <- strptime(conc$Datetime, format="%d/%m/%Y %H:%M", tz = "Australia/Brisbane")
conc$Datetime <- round(conc$Datetime, units = "hours")

# convert less than values to the corresponding threshold values
conc$TSS[conc$NH4=="<0.002"] = 0.002
conc$NOx[conc$NOx=="<0.001"] = 0.001

# convert NOx, NH4 and DIn data from character to numeric
indx <- sapply(conc[,c(2:4)], is.character)
conc[indx] <- lapply(conc[indx], function(x) as.numeric(x))

# convert blank data to NA
conc$Equipment_contamination[conc$Equipment_contamination == "" ] <- NA

# add lab condition
conc$Lab <- "All"

# create separate TSS and Nox concentration data
nox <- conc[,c(-3, -4)]
din <- conc[,c(-2, -3)]
dinraw <- din[, c(1, 2)]
nhraw <- conc[, c(1, 3)]

## convert days in category conditions
thresholds <- read.csv("time thresholds.csv")

# Inputs: chose a threshold type, here "Guideline"
# Obtain the threshold value for different conditions
toanalysis_optimal <- thresholds$Guideline[thresholds$Condition == "Days from sample to analysis" & thresholds$State == "Optimal"]
toanalysis_suboptimal <- thresholds$Guideline[thresholds$Condition == "Days from sample to analysis" & thresholds$State == "Suboptimal"]

tofilter_optimal_ma <- thresholds$Guideline[thresholds$Condition == "Days to filtering" & thresholds$State == "Optimal" & thresholds$Sampling == "MA"]
tofilter_suboptimal_ma <- thresholds$Guideline[thresholds$Condition == "Days to filtering" & thresholds$State == "Suboptimal" & thresholds$Sampling == "MA"]

tofilter_optimal_as <- thresholds$Guideline[thresholds$Condition == "Days to filtering" & thresholds$State == "Optimal" & thresholds$Sampling == "AS"]
tofilter_suboptimal_as <- thresholds$Guideline[thresholds$Condition == "Days to filtering" & thresholds$State == "Suboptimal" & thresholds$Sampling == "AS"]

nox$Days_from_sample_to_analysis = with(nox, ifelse(nox$Days_from_sample_to_analysis <= toanalysis_optimal, "Optimal", 
                                                    ifelse(nox$Days_from_sample_to_analysis <= toanalysis_suboptimal, "Suboptimal", "Bad")))


nox$Days_to_filtering = with(nox, ifelse((nox$Equipment_contamination == "MA" & nox$Days_to_filtering <= tofilter_optimal_ma) | (nox$Equipment_contamination == "AS" & nox$Days_to_filtering <= tofilter_optimal_as), "Optimal", 
                                         ifelse((nox$Equipment_contamination == "MA" & nox$Days_to_filtering <= tofilter_suboptimal_ma) | (nox$Equipment_contamination == "AS" & nox$Days_to_filtering <= tofilter_suboptimal_as) | (nox$Equipment_contamination == "contaminated" & nox$Days_to_filtering <= tofilter_suboptimal_as), "Suboptimal", "Bad")))

nox$comb <- paste(paste("Flow at sample time", nox$Flow_at_sample_time, sep=" "),
                  paste("Equipment contamination", nox$Equipment_contamination, sep=" "),
                  paste("Sample Temp", nox$Sample_Temp, sep=" "), 
                  paste("Days from sample to analysis", nox$Days_from_sample_to_analysis, sep=" "), 
                  paste("Days to filtering", nox$Days_to_filtering, sep=" "), 
                  paste("Lab", nox$Lab, sep=" "),
                  sep = ",")

### create concentration uncertainty library
bounds <- read.csv("condition uc bounds.csv")
conditions <- unique(bounds$Condition)

nalookup <- read.csv("na lookup.csv")
bounds <- plyr::rbind.fill(bounds, nalookup)

bounds$condstat <- paste(bounds$Condition, bounds$State, sep=" ")

bounds$State <- addNA(bounds$State)
bound2 <- reshape(bounds, idvar=c("Condition", "State", "condstat"), timevar = "Bound", direction="wide" )

expandfun <- function(df, bound){
  df1 <- subset(df, df$Condition == conditions[1])[,bound]
  df2 <- subset(df, df$Condition == conditions[2])[,bound]
  df3 <- subset(df, df$Condition == conditions[3])[,bound]
  df4 <- subset(df, df$Condition == conditions[4])[,bound]
  df5 <- subset(df, df$Condition == conditions[5])[,bound]
  df6 <- subset(df, df$Condition == conditions[6])[,bound]
  newdf <- expand.grid(df1, df2, df3, df4, df5, df6)
}

# This is where to change for Qian's work
# Combine the uncertain bounds in a multiplicative way.
combfun <- function(df){
  combound <- round(sqrt((df$Var1**2) + (df$Var2**2) + (df$Var3**2) + (df$Var4**2) + (df$Var6**2)), 2)
}

Calcol <- c("Median.Upper", "Median.Lower", "Absolute.extremes.Upper", "Absolute.extremes.Lower", "Refined.extremes.Upper", "Refined.extremes.Lower", "ID1.Upper", "ID1.Upper","ID1.Lower","ID2.Upper","ID2.Lower","ID3.Upper","ID3.Lower","ID4.Upper","ID4.Lower","ID5.Upper","ID5.Lower","ID6.Upper","ID6.Lower","ID7.Upper","ID7.Lower","ID8.Upper","ID8.Lower","ID9.Upper","ID9.Lower","ID10.Upper","ID10.Lower","ID11.Upper","ID11.Lower","ID12.Upper","ID12.Lower")

temp <- list()

# Calculate the uncertainty bounds for different Calcol category
for (i in Calcol){
  temp[[i]] <- expandfun(bound2, i)
  temp[[i]] <- combfun(temp[[i]])
}

comball <- do.call(cbind.data.frame, temp)

expcond <- expandfun(bound2, "condstat")
comball <- cbind(expcond, comball)

comball$comb <- paste(comball$Var1, comball$Var2, comball$Var4, comball$Var3, comball$Var5, "Lab All", sep = ",") #slight variation in the order of the variables, better code needed

# select the Median.Upper and Median.Lower to represent the uncertainty bounds
comball_sub <- comball[,c("comb", "Median.Upper", "Median.Lower")]
noxcomb <- merge(nox, comball_sub, by="comb", all= TRUE)
noxcomb <- noxcomb[!is.na(noxcomb$Datetime),] 

noxcomb <- noxcomb[,c("Datetime", "NOx", "Median.Upper", "Median.Lower")]
noxcomb <- noxcomb[order(noxcomb$Datetime),]
rownames(noxcomb) <- 1:nrow(noxcomb)

# Qian to change the code here for preparing the concdf
noxuncup <- noxcomb[,c(1,3)]
names(noxuncup)<-c('Datetime', 'concuncbound')
noxunclow <- noxcomb[,c(1,4)]
names(noxunclow)<-c('Datetime', 'concuncbound')
noxraw <- mutate(noxcomb, NOx=NOx)[,c(1:2)] ##remove extra column to keep just datetime and nox.
names(noxraw)<-c('Datetime', 'concuncbound')

## setting other input data
yearlyseconds <- c(31532400,31532400,31618800,31532400,31532400,31532400,31618800,31532400,31532400)

MinNH4 = 0.002
MinNOx = 0.001
MinDIN = MinNH4 + MinNOx


# ### scenarios data

# Different unc scenarios of flow and concentration
# flow_scenarios <- list(flow = flow, flowup = flowup, flowlow = flowlow)
# nox_scenarios <- list(nox = noxcomb, noxup = noxup, noxlow = noxlow)

# Estimate and plot combined uncertainty -----------------------------------------------------------

## generating combined uncertainty for upper bound
# noxmixlow <- samplefreqfun(flowdf = flow, concdf = noxraw, seconds = yearlyseconds, minconc = MinNOx,
#                            loadfun = "linear",  biasfun = 1, randomfun = "randomyear", randomseedprop = 1.0,
#                            randomseeds = 10, randomtimes = 2, flowunc = flowunclow, concunc = noxunclow,
#                            covar = FALSE, boundtype = "lower")

# write.csv(noxmixlow, 'output/noxmixlow_nocovar_100.csv')

# noxmixlowcovar <- samplefreqfun(flowdf = flow, concdf = noxraw, seconds = yearlyseconds, minconc = MinNOx,
#                                       loadfun = "linear",  biasfun = 1, randomfun = "randomyear", randomseedprop = 0.7,
#                                       randomseeds = 10, randomtimes = 1000, flowunc = flowunclow,
#                                       concunc = noxunclow, covar = TRUE, boundtype = "lower")
# write.csv(noxmixlowcovar, 'output/noxmixlow_covar_70.csv')

# noxmixup <- samplefreqfun(flowdf = flow, concdf = noxraw, seconds = yearlyseconds, minconc = MinNOx,
#                            loadfun = "linear",  biasfun = 1, randomfun = "randomyear", randomseedprop = 1.0,
#                            randomseeds = 10, randomtimes = 2, flowunc = flowuncup, concunc = noxuncup,
#                            covar = FALSE, boundtype = "upper")
# 
# write.csv(noxmixup, 'output/noxmixup_nocovar_100.csv')

# noxmixupcovar <- samplefreqfun(flowdf = flow, concdf = noxraw, seconds = yearlyseconds, minconc = MinNOx,
#                                 loadfun = "linear",  biasfun = 1, randomfun = "randomyear", randomseedprop = 0.7,
#                                 randomseeds = 10, randomtimes = 1000, flowunc = flowuncup,
#                                 concunc = noxuncup, covar = TRUE, boundtype = "upper")
# write.csv(noxmixupcovar, 'output/noxmixup_covar_70.csv')

# ## generating combined uncertainty for lower bound for din
dinmixlow <- samplefreqfun(flowdf = flow, concdf = dinraw, seconds = yearlyseconds, minconc = MinDIN,
                                      loadfun = "linear",  biasfun = 1, randomfun = "randomyear", randomseedprop = 1.0,
                                      randomseeds = 10, randomtimes = 2, flowunc = flowunclow,
                                      concunc = noxunclow, covar = 0.1, boundtype = "lower")

write.csv(dinmixlow, 'output/resample_freq100_new_assumption/dinmixlow_covar_100.csv')

# dinmixlowcovar <- samplefreqfun(flowdf = flow, concdf = dinraw, seconds = yearlyseconds, minconc = MinDIN,
#                                            loadfun = "linear",  biasfun = 1, randomfun = "randomyear", randomseedprop = 0.7,
#                                            randomseeds = 10, randomtimes = 1000, flowunc = flowunclow,
#                                            concunc = noxunclow, covar = TRUE, boundtype = "lower")
# write.csv(dinmixlowcovar, 'output/dinmixlow_covar_70.csv')

## generating combined uncertainty for upper bound for din
dinmixup <- samplefreqfun(flowdf = flow, concdf = dinraw, seconds = yearlyseconds, minconc = MinDIN,
                           loadfun = "linear",  biasfun = 1, randomfun = "randomyear", randomseedprop = 1.0,
                           randomseeds = 10, randomtimes = 2, flowunc = flowuncup,
                           concunc = noxuncup, covar = 0.1, boundtype = "upper")

write.csv(dinmixup, 'output/resample_freq100_new_assumption/dinmixup_covar_100.csv')

# dinmixupcovar <- samplefreqfun(flowdf = flow, concdf = dinraw, seconds = yearlyseconds, minconc = MinDIN,
#                                 loadfun = "linear",  biasfun = 1, randomfun = "randomyear", randomseedprop = 0.7,
#                                 randomseeds = 10, randomtimes = 1000, flowunc = flowuncup,
#                                 concunc = noxuncup, covar = TRUE, boundtype = "upper")
# write.csv(dinmixupcovar, 'output/dinmixup_covar_70.csv')

 
# noxmix <- loadlinear(flow, noxcomb, yearlyseconds, MinNOx)
# noxmix$noxuploadmax <- noxmixup$load_max
# noxmix$noxuploadmin <- noxmixup$load_min
# noxmix$noxlowloadmax <- noxmixlow$load_max
# noxmix$noxlowloadmin <- noxmixlow$load_min
# noxmix

## plotting combined uncertainty
# plotdata <- noxmix
# plotdata$Year <- as.integer(row.names(plotdata))
# plotdata$noxlowloadmean <- (plotdata$noxlowloadmax + plotdata$noxlowloadmin)/2
# plotdata$noxuploadmean <- (plotdata$noxuploadmax + plotdata$noxuploadmin)/2
# 
# noxtop <- ggplot(plotdata, aes(x=as.factor(Year)))+
#   geom_col(aes(y=load), fill = "wheat3")+
#   geom_errorbar(aes(ymin = noxlowloadmin, ymax = noxuploadmax), colour="black", width=0.3)+
#   xlab(NULL)+
#   ylab("NOx Load (t/yr)")+
#   theme_bw()
# 
# noxmid <- ggplot(plotdata, aes(x=as.factor(Year)))+
#   geom_point(aes(y=noxlowloadmin/load), colour = "steel blue")+
#   geom_point(aes(y=noxuploadmax/load), colour = "tomato")+
#   geom_hline(yintercept=1, linetype="dashed") +
#   # ylim (0,2.5)+
#   labs(y = "Load Ratio (Subsamples/full samples)", x="Hydrological year")+
#   scale_colour_brewer("Beale's Bias Correction", palette="Dark2")+
#   theme_bw()
# 
# grid.arrange(noxtop, noxmid, ncol=2, nrow=1)
