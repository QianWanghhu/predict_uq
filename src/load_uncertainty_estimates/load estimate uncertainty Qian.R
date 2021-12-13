## Estimating uncertainty in constituent annual load estimation based on interpolation and beale ratio methods
## Baihua Fu  ## Edited by Qian Wang 
## 9/12/2021  ## 13/12/2021


# Set work space----
### set word directory
setwd("E:\\cloudStor\\Projects\\predict_uq \\Projects\\src\\load_uncertainty_estimates")
### load libraries
library(dplyr)
library(ggplot2)
library(zoo)
library(gridExtra)


# Generate functions ---------------------------------------------------------------

### Calculate hydrology year function
as.year <- function(x) as.integer(as.yearmon(x) - 6/12)

### Linear interpolation function 
loadlinear <- function(flowdf, concdf, seconds, minconc){
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


### Beale ratio function 

loadbeale <- function(flowdf, concdf, seconds, biasfun){
  flowzoo <- read.zoo(flowdf, header = TRUE)

  conczoo <- read.zoo(concdf, header = TRUE, aggregate = mean) #take average for same time step
  
  # merge flow and concentration data
  flowconczoo <- merge(flowzoo, conczoo)
  names(flowconczoo) <- c("flow", "conc")
  
  # convert data from zoo to dataframe format
  flowconc <- data.frame(Datetime = index(flowconczoo), as.data.frame(flowconczoo))
  
  # create hydrology year column
  flowconc$hyear <- as.year(flowconc$Datetime)
  flowconc$hyear[c(1:10)] <- flowconc$hyear[11] ## wrap the first 10 hours to year 2009, due to unresolved issue to make time zone match.
  
  # calculate load for each hour time step
  flowconc$load <- flowconc$flow*flowconc$conc

  # calculate mean instant load Lm
  Lm <- flowconc %>%
    group_by(hyear) %>% 
    dplyr::summarise_at(c("load", "flow"), mean, na.rm=TRUE)
  
  # create identifier (A) to tag dates with/without concentration data
  flowconc$A <- 1
  flowconc$A[is.na(flowconc$conc)] <- NA

  # calculate AQ: 
  flowconc$AQ <- flowconc$flow*flowconc$A
  
  # calculate flow ratio adjustment factor Fa
  Fa <- flowconc %>%
    group_by(hyear) %>%
    dplyr::summarise_at(c("flow", "AQ"), mean, na.rm=TRUE) %>%
    mutate(Fa = flow/AQ) #mean of annual flow divided by the mean of flows with concentration samples
  
  # calculate ACQ^2
  flowconc$ACQ2 <- flowconc$A*flowconc$conc*flowconc$flow*flowconc$flow
  
  # calculate AQ^2
  flowconc$AQ2 <- flowconc$A*flowconc$flow*flowconc$flow

  # calculate ACQ
  flowconc$ACQ <- flowconc$A*flowconc$conc*flowconc$flow
  
  # calculate Q^2
  flowconc$Q2 <- flowconc$flow*flowconc$flow

  # calculate bias correction factors

  if(biasfun == 1){
    Bf <- flowconc %>%
      group_by(hyear) %>%
      dplyr::summarise(across(.cols = c("AQ", "ACQ", "Q2"), ~ mean(.x, na.rm = TRUE)), 
                       across(.cols = c("ACQ2", "AQ2", "A"),  ~ sum(.x, na.rm = TRUE)), 
                       across(.cols = c("flow"), ~ length(.x))) %>%
      mutate(Slq = (1/(A - 1))*(ACQ2 - A*AQ*ACQ)) %>%
      mutate(Sq2 = (1/(A - 1))*(AQ2 - A*Q2)) %>%
      mutate(Bf = (1+(1/A)*(Slq/(ACQ*AQ)))/(1+(1/A)*(Sq2/AQ2)))
  } else if (biasfun == 2){
    Bf <- flowconc %>%
      group_by(hyear) %>%
      dplyr::summarise(across(.cols = c("AQ", "ACQ", "Q2"), ~ mean(.x, na.rm = TRUE)), 
                       across(.cols = c("ACQ2", "AQ2", "A"),  ~ sum(.x, na.rm = TRUE)), 
                       across(.cols = c("flow"), ~ length(.x))) %>%
      mutate(Slqb2 = (1/(A - 1))*(ACQ - A*AQ*ACQ)) %>%
      mutate(Sq2b2 = (1/(A - 1))*(AQ - A*Q2)) %>%
      mutate(Bf = (1+(1/A)*(Slqb2/(ACQ*AQ)))/(1+(1/A)*(Sq2b2/AQ2)))
  } else if (biasfun == 3){
    Bf <- flowconc %>%
      group_by(hyear) %>%
      dplyr::summarise(across(.cols = c("AQ", "ACQ", "Q2"), ~ mean(.x, na.rm = TRUE)), 
                       across(.cols = c("ACQ2", "AQ2", "A"),  ~ sum(.x, na.rm = TRUE)), 
                       across(.cols = c("flow"), ~ length(.x))) %>%
      mutate(Slq = (1/(A - 1))*(ACQ2 - A*AQ*ACQ)) %>%
      mutate(Sq2 = (1/(A - 1))*(AQ2 - A*Q2)) %>%
      mutate(Bf = (1+((1-A/flow)/A)*(Slq/(ACQ*AQ)))/(1+((1-A/flow)/A)*(Sq2/AQ2)))
  }
  
  colnames(Bf)[colnames(Bf) == "flow"] <- "days"

  # join inputs for Beale ratio equations together
  Beale <- left_join(Lm, Fa, by=c("hyear", "flow"))
  Beale <- left_join(Beale, Bf, by=c("hyear", "AQ"))
  
  # calculate annual load
  Bealeload <- as.data.frame(Beale)
  Bealeload$bealeload = round(Bealeload$load*Bealeload$Fa*Bealeload$Bf*seconds * 0.000001, 3)
  Bealeload$Biasfactor = biasfun
  
  return(Bealeload)
  
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

samplefreqfun <- function(flowdf, concdf, seconds,minconc, loadfun, biasfun, randomfun, randomseeds, randomseedprop, randomtimes){
  # get random samples
  if (randomfun == "randomall"){
    concsamplefreq <- randomall(concdf, randomseeds, randomtimes)
  } else if (randomfun == "randomyear"){
    concsamplefreq <- randomyear(concdf, randomseedprop, randomtimes)
  }
  
  hyear <- as.year(flowdf$Datetime)
  
  conclist <- list()
  
  if (loadfun == "linear"){
    for (i in (unique(concsamplefreq$Obs))){
      concsub <- subset(concsamplefreq, concsamplefreq$Obs == i)
      concsub <- as.data.frame(concsub[,c("Datetime", "conc")])
      conclist[[i]] <- loadlinear(flowdf, concsub, seconds, minconc)
      conclist[[i]]$hyear <- unique(hyear) ## TODO
    }
  } else if (loadfun == "beale"){
    for (i in (unique(concsamplefreq$Obs))){
      concsub <- subset(concsamplefreq, concsamplefreq$Obs == i)
      concsub <- as.data.frame(concsub[,c("Datetime", "conc")])
      conclist[[i]] <- loadbeale(flowdf, concsub, seconds, biasfun)
      conclist[[i]]$hyear <- unique(hyear) ## TODO
    }
  }
  
  concsampleall <- bind_rows(conclist, .id = "samples")
  
  concsampleminmax <- concsampleall %>%
    group_by(hyear, flow) %>%
    dplyr::summarise(across(.cols = c("load"),list (min = min, max = max), na.rm=TRUE))
  
  return(concsampleminmax)
  
}




# Import data -------------------------------------------------------------
### Import and format flow data
# import flow data
flow <- read.csv("126001A_flows.csv", skip = 1, header = F)

# tidy flow data to date, flow, qc format
colnames(flow) <- rep(c("Datetime", "flow", "QC","X"), 10)
flow <- rbind(flow[,c(1:3)], flow[,c(5:7)], flow[,c(9:11)], flow[,c(13:15)], flow[,c(17:19)], 
              flow[,c(21:23)], flow[,c(25:27)], flow[,c(29:31)], flow[,c(33:35)], flow[,c(37:39)])

flow <- na.omit(flow) 

flow$Datetime <- strptime(flow$Datetime, format="%d/%m/%Y %H:%M", tz = "Australia/Brisbane")


# create flow upper bounds and lower bounds data
flowup <- mutate(flow, flow = case_when(
  QC == 9 ~ flow*1.05,
  QC == 10 ~ flow*1.1,
  QC == 20 ~ flow*1.2,
  QC == 30 ~ flow*1.5,
  QC == 60 ~ flow*1.75,
  TRUE ~ flow))

flowlow <- mutate(flow, flow = case_when(
  QC == 9 ~ flow*0.95,
  QC == 10 ~ flow*0.9,
  QC == 20 ~ flow*0.8,
  QC == 30 ~ flow*0.5,
  QC == 60 ~ flow*0.25,
  TRUE ~ flow))

flowqc <- flow
flow <- flow[,-3]
flowup <- flowup[,-3]
flowlow <- flowlow[,-3]


### import and format concentration data
conc <- read.csv("126001A_concentrations_conditions.csv", header=TRUE, stringsAsFactors=FALSE) #TSS mg/L, NOx mg/L

# take only time, TSS, NOX, conditoins
conc <- conc[,c(3:10)]
colnames(conc) <- c("Datetime", "TSS", "NOx", "Flow_at_sample_time", "Equipment_contamination", "Sample_Temp", "Days_from_sample_to_analysis", "Days_to_filtering")

# round the concentration data to nearest hour
conc$Datetime <- strptime(conc$Datetime, format="%d/%m/%Y %H:%M", tz = "Australia/Brisbane")
conc$Datetime <- round(conc$Datetime, units = "hours")

# convert less than values to the corresponding threshold values
conc$TSS[conc$TSS=="<1"] = 1
conc$NOx[conc$NOx=="<0.001"] = 0.001

# convert TSS and NOx data from character to numeric
indx <- sapply(conc[,c(2:3)], is.character)
conc[indx] <- lapply(conc[indx], function(x) as.numeric(x))

# convert blank data to NA
conc$Equipment_contamination[conc$Equipment_contamination == "" ] <- NA

# add lab condition
conc$Lab <- "All"

# create separate TSS and Nox concentration data
tss <- conc[,-3]
nox <- conc[,-2]

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

combfun <- function(df){
  combound <- round((1+df$Var1)*(1+df$Var2)*(1+df$Var3)*(1+df$Var4)*(1+df$Var6)-1, 2)
}

Calcol <- c("Median.Upper", "Median.Lower", "Absolute.extremes.Upper", "Absolute.extremes.Lower", "Refined.extremes.Upper", "Refined.extremes.Lower", "ID1.Upper", "ID1.Upper","ID1.Lower","ID2.Upper","ID2.Lower","ID3.Upper","ID3.Lower","ID4.Upper","ID4.Lower","ID5.Upper","ID5.Lower","ID6.Upper","ID6.Lower","ID7.Upper","ID7.Lower","ID8.Upper","ID8.Lower","ID9.Upper","ID9.Lower","ID10.Upper","ID10.Lower","ID11.Upper","ID11.Lower","ID12.Upper","ID12.Lower")

temp <- list()

for (i in Calcol){
  temp[[i]] <- expandfun(bound2, i)
  temp[[i]] <- combfun(temp[[i]])
}

comball <- do.call(cbind.data.frame, temp)

expcond <- expandfun(bound2, "condstat")
comball <- cbind(expcond, comball)

comball$comb <- paste(comball$Var1, comball$Var2, comball$Var4, comball$Var3, comball$Var5, "Lab All", sep = ",") #slight variation in the order of the variables, better code needed

comball_sub <- comball[,c("comb", "Median.Upper", "Median.Lower")]
noxcomb <- merge(nox, comball_sub, by="comb", all= TRUE)
noxcomb <- noxcomb[!is.na(noxcomb$Datetime),] 

noxcomb <- noxcomb[,c("Datetime", "NOx", "Median.Upper", "Median.Lower")]
noxcomb <- noxcomb[order(noxcomb$Datetime),]
rownames(noxcomb) <- 1:nrow(noxcomb)

noxup <- mutate(noxcomb, NOx = NOx*(1+Median.Upper))[,c(1,2)]

noxlow <- mutate(noxcomb, NOx = NOx*(1+Median.Lower))[,c(1,2)]

noxcomb <- noxcomb[,c(1:2)] ## remove extra column to keep just datetime and nox. 
noxup <- noxup[,c(1:2)]
noxlow <- noxlow[,c(1:2)]

### setting other input data
yearlyseconds <- c(31532400,31532400,31618800,31532400,31532400,31532400,31618800,31532400,31532400,31532400)

MinTSS = 1
MinNOx = 0.001

### scenarios data

flow_scenarios <- list(flow = flow, flowup = flowup, flowlow = flowlow)
nox_scenarios <- list(nox = noxcomb[,c(1,2)], noxup = noxup[,c(1,2)], noxlow = noxlow[,c(1,2)])


# Estimate and plot flow uncertainty ------------------------------------------------------------

# test impacts of flow, linear  
noxout <- list()
for (i in names(flow_scenarios)){
  noxout[[i]] <- loadlinear(flow_scenarios[[i]], noxcomb, yearlyseconds, MinNOx)
  noxout[[i]]$flowscenario <- i
  noxout[[i]]$hyear <- row.names(noxout[[i]])
}

# alternatively: test impacts of flow, beale 
noxout <- list()
for (i in names(flow_scenarios)){
  noxout[[i]] <- loadbeale(flow_scenarios[[i]], noxcomb, yearlyseconds, biasfun = 1)
  noxout[[i]]$flowscenario <- i
}


# plotting flow uncertainty
flowscennox <- as.data.frame(noxout)

flowscennox %>%
  group_by(flow.hyear) %>%
  summarise(flowupbound = flowup.load/flow.load -1, flowlowbound = 1-  flowlow.load/flow.load)

noxtop <- ggplot(flowscennox, aes(x=flow.hyear))+
  geom_col(aes(y=flow.load), fill = "wheat3")+
  geom_errorbar(aes(ymin = flowlow.load, ymax = flowup.load), colour="black", width=0.3)+
  xlab("Hydrological year")+
  ylab("NOx Load (t/yr)")+
  theme_bw()

noxmid <- ggplot(flowscennox, aes(x=flow.hyear))+
  geom_point(aes(y=flowup.load/flow.load), colour = "steel blue")+
  geom_point(aes(y=flowlow.load/flow.load), colour = "tomato")+
  geom_hline(yintercept=1, linetype="dashed") +
  ylim (0.5,1.5)+
  labs(y = "Load Ratio (Flow Bounds/Flow Obs)", x="Hydrological year")+
  scale_colour_brewer("Beale's Bias Correction", palette="Dark2")+
  theme_bw()

noxbot <- ggplot(flowscennox, aes(x=flow.flow))+
  geom_point(aes(y=flowup.load/flow.load), colour = "steel blue")+
  geom_point(aes(y=flowlow.load/flow.load), colour = "tomato")+
  geom_hline(yintercept=1, linetype="dashed") +
  ylim (0.5,1.5)+
  labs(y = "Load Ratio (Flow Bounds/Flow Obs)", x="Annual Mean Hourly Flow (Cumecs)")+
  scale_colour_brewer("Beale's Bias Correction", palette="Dark2")+
  theme_bw()

grid.arrange(noxtop,  noxmid,  noxbot, ncol=3, nrow=1)#, widths=c(1, 1), heights=c(1, 1,1))


# Estimate and plot concentration uncertainty ---------------------------------------------------

# test impacts of concentration, linear  
noxout <- list()

for (i in names(nox_scenarios)){
  noxout[[i]] <- loadlinear(flow, nox_scenarios[[i]], yearlyseconds, MinNOx)
  noxout[[i]]$concscenario <- i
  noxout[[i]]$hyear <- row.names(noxout[[i]])
}

# alternatively, test impacts of concentration, beale  
noxout <- list()

for (i in names(nox_scenarios)){
  noxout[[i]] <- loadbeale(flow, nox_scenarios[[i]], yearlyseconds, biasfun = 1)
  noxout[[i]]$concscenario <- i
}


# plotting concentration uncertainty
flowscennox <- as.data.frame(noxout)

noxtop <- ggplot(flowscennox, aes(x=nox.hyear))+
  geom_col(aes(y=nox.load), fill = "wheat3")+
  geom_errorbar(aes(ymin = noxlow.load, ymax = noxup.load), colour="black", width=0.3)+
  xlab("Hydrological year")+
  ylab("NOx Load (t/yr)")+
  theme_bw()

noxmid <- ggplot(flowscennox, aes(x=nox.hyear))+
  geom_point(aes(y=noxup.load/nox.load), colour = "steel blue")+
  geom_point(aes(y=noxlow.load/nox.load), colour = "tomato")+
  geom_hline(yintercept=1, linetype="dashed") +
  labs(y = "Load ratio ([C] bounds/[C] best-guess)", x="Hydrological year")+
  scale_colour_brewer("Beale's Bias Correction", palette="Dark2")+
  theme_bw()

noxbot <- ggplot(flowscennox, aes(x=nox.flow))+
  geom_point(aes(y=noxup.load/nox.load), colour = "steel blue")+
  geom_point(aes(y=noxlow.load/nox.load), colour = "tomato")+
  geom_hline(yintercept=1, linetype="dashed") +
  labs(y = "Load ratio ([C] bounds/[C] best-guess)", x="Annual Mean Hourly Flow (Cumecs)")+
  scale_colour_brewer("Beale's Bias Correction", palette="Dark2")+
  theme_bw()

grid.arrange(noxtop, noxmid, ncol=2, nrow=1)#, widths=c(1,1))


# Estimate and plot load estimation function uncertainty ---------------------------------------------------
## estimating load estimation function uncertainty

noxlinear <- loadlinear(flow, noxcomb, yearlyseconds, MinNOx)
noxbeale1 <- loadbeale(flow, noxcomb, yearlyseconds, biasfun = 1)
noxbeale2 <- loadbeale(flow, noxcomb, yearlyseconds, biasfun = 2)
noxbeale3 <- loadbeale(flow, noxcomb, yearlyseconds, biasfun = 3)

colnames(noxbeale2) <- colnames(noxbeale1)
noxbeale <- rbind(noxbeale1, noxbeale2)


# plotting estimation function uncertainty
noxlibe <- merge(noxlinear, noxbeale, by.x="row.names", by.y = "hyear")

noxtop <- ggplot(noxlibe, aes(x=as.factor(Row.names)))+
  geom_col(aes(y=load.x/2000), fill = "wheat3")+
  geom_point(aes(y = bealeload/1000, colour = as.factor(Biasfactor)), size=2)+
  xlab("Hydrological year")+
  ylab("NOx Load (kt/yr)")+
  scale_colour_brewer("Beale's Bias Correction", palette="Dark2")+
  theme_bw()+
  guides(colour=FALSE)+
  theme(axis.title.x = element_blank())

noxmid <- ggplot(noxlibe, aes(x=as.factor(Row.names)))+
  geom_point(aes(y = (bealeload/load.x), colour = as.factor(Biasfactor)),  size=2)+
  geom_hline(yintercept=1, linetype="dashed") +
  ylim (0,3)+
  labs(y = "NOx Load Ratio (Beale/Linear)", x="Hydrological year")+
  scale_colour_brewer("Beale's Bias Correction", palette="Dark2")+
  theme_bw()+
  theme(legend.position = c(0.8,0.8))+
  guides(colour=FALSE)


noxbot <- ggplot(noxlibe, aes(x=flow.x))+
  geom_point(aes(y = (bealeload/load.x), colour = as.factor(Biasfactor)),  size=2)+
  geom_hline(yintercept=1, linetype="dashed") +
  ylim (0,3)+
  labs(y = "NOx Load Ratio (Beale/Linear)", x="Annual Mean Hourly Flow (Cumecs)")+
  scale_colour_brewer("Beale's Bias Correction", palette="Dark2")+
  theme_bw()+
  theme(legend.position = c(0.8,0.8))+
  guides(colour=guide_legend(title="Bias Correction", ncol=2))

grid.arrange(noxtop,  noxmid, noxbot, ncol=3, nrow=1)


# Estimate and plot sample frequency uncertainty ----------------------------------------------------------

## estimating sample frequency uncertainty with linear estimation function
nox_base_linear_randomyear <- samplefreqfun(flowdf = flow, concdf = noxcomb, seconds = yearlyseconds, minconc = MinNOx, loadfun = "linear",  biasfun = 1, randomfun = "randomyear", randomseedprop = 0.7, randomseeds = 10, randomtimes = 5000)

noxlinear <- loadlinear(flow, noxcomb, yearlyseconds, MinNOx)

nox_base_linear_randomyear$conc <- noxlinear$conc
nox_base_linear_randomyear$load <- noxlinear$load

## plotting sample frequency uncertainty
noxrandom <- as.data.frame(nox_base_linear_randomyear)

noxtop <- ggplot(noxrandom, aes(x=as.factor(hyear)))+
  geom_col(aes(y=load), fill = "wheat3")+
  geom_errorbar(aes(ymin = load_min, ymax = load_max), colour="black", width=0.3)+
  xlab(NULL)+
  ylab("NOx Load (t/yr)")+
  theme_bw()

noxmid <- ggplot(noxrandom, aes(x=as.factor(hyear)))+
  geom_point(aes(y=load_min/load), colour = "steel blue")+
  geom_point(aes(y=load_max/load), colour = "tomato")+
  geom_hline(yintercept=1, linetype="dashed") +
  ylim (0,2.5)+
  labs(y = "Load Ratio (Flow Bounds/Flow Obs)", x="Hydrological year")+
  scale_colour_brewer("Beale's Bias Correction", palette="Dark2")+
  theme_bw()

noxbot <- ggplot(noxrandom, aes(x=flow))+
  geom_point(aes(y=load_min/load), colour = "steel blue")+
  geom_point(aes(y=load_max/load), colour = "tomato")+
  geom_hline(yintercept=1, linetype="dashed") +
  ylim (0,2.5)+
  labs(y = "Load Ratio (Flow Bounds/Flow Obs)", x="Annual Mean Hourly Flow (Cumecs)")+
  scale_colour_brewer("Beale's Bias Correction", palette="Dark2")+
  theme_bw()

grid.arrange(noxtop,  noxmid,  noxbot, ncol=3, nrow=1)


## estimating sample frequency uncertainty with beale ratio function
nox_base_beale_randomyear <- samplefreqfun(flowdf = flow, concdf = nox, seconds = yearlyseconds, minconc = MinNOx, loadfun = "beale",  biasfun = 1, randomfun = "randomyear", randomseedprop = 0.7, randomseeds = 10, randomtimes = 5000)

noxlinear <- loadlinear(flow, noxcomb, yearlyseconds, MinNOx)

nox_base_beale_randomyear$conc <- noxlinear$conc
nox_base_beale_randomyear$load <- noxlinear$load

## plotting frequency uncertainty
noxrandombeale <- as.data.frame(nox_base_beale_randomyear)

noxtop <- ggplot(noxrandombeale, aes(x=as.factor(hyear)))+
  geom_col(aes(y=load), fill = "wheat3")+
  geom_errorbar(aes(ymin = load_min, ymax = load_max), colour="black", width=0.3)+
  xlab(NULL)+
  ylab("NOx Load (t/yr)")+
  theme_bw()

noxmid <- ggplot(noxrandombeale, aes(x=as.factor(hyear)))+
  geom_point(aes(y=load_min/load), colour = "steel blue")+
  geom_point(aes(y=load_max/load), colour = "tomato")+
  geom_hline(yintercept=1, linetype="dashed") +
  ylim (0,2.5)+
  labs(y = "Load Ratio (Subsamples/full samples)", x="Hydrological year")+
  scale_colour_brewer("Beale's Bias Correction", palette="Dark2")+
  theme_bw()

noxbot <- ggplot(noxrandombeale, aes(x=flow))+
  geom_point(aes(y=load_min/load), colour = "steel blue")+
  geom_point(aes(y=load_max/load), colour = "tomato")+
  geom_hline(yintercept=1, linetype="dashed") +
  ylim (0,2.5)+
  labs(y = "Load Ratio (Subsamples/full samples)", x="Annual Mean Hourly Flow (Cumecs)")+
  scale_colour_brewer("Beale's Bias Correction", palette="Dark2")+
  theme_bw()


grid.arrange(tsstop, noxtop, tssmid, noxmid, ncol=2, nrow=2)


# Estimate and plot combined uncertainty -----------------------------------------------------------

## generating combined uncertainty

noxmixup <- samplefreqfun(flowdf = flowup, concdf = noxup, seconds = yearlyseconds, minconc = MinNOx, loadfun = "linear",  biasfun = 1, randomfun = "randomyear", randomseedprop = 0.7, randomseeds = 10, randomtimes = 5000)

noxmixlow <- samplefreqfun(flowdf = flowlow, concdf = noxlow, seconds = yearlyseconds, minconc = MinNOx, loadfun = "linear",  biasfun = 1, randomfun = "randomyear", randomseedprop = 0.7, randomseeds = 10, randomtimes = 5000)

noxmix <- loadlinear(flow, noxcomb, yearlyseconds, MinNOx)
noxmix$noxuploadmax <- noxmixup$load_max
noxmix$noxuploadmin <- noxmixup$load_min
noxmix$noxlowloadmax <- noxmixlow$load_max
noxmix$noxlowloadmin <- noxmixlow$load_min
noxmix

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



