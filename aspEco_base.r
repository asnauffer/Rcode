# reads ASP data from a csv file
# runs SnowMelt from EcoHydRology using T & P data

if(R.version[1]=="x86_64-pc-linux-gnu") {
  library("EcoHydRology", lib.loc="/net/home/asnauffer/R/x86_64-pc-linux-gnu-library/3.1")
  library("R.matlab", lib.loc="/net/home/asnauffer/R/x86_64-pc-linux-gnu-library/3.1")
  setwd("/net/home/asnauffer/Documents/PhD/Rcode")
} else {
  library("EcoHydRology", lib.loc="C:/Users/drew/Documents/R/win-library/3.1")
  library("R.matlab", lib.loc="C:/Users/drew/Documents/R/win-library/3.1")
  setwd("C:/Users/Drew/Documents/PhD/Rcode")
}

source('SnowMelt2L.R')
source('SnowAccum.R')

skipreload <- FALSE
if(skipreload){
  print("skipping data load & variable init")
} else {
  print("loading and initializing...")
  asplatlon <- data.frame(readMat("../matfiles/ASPalllatlon.mat"))
  stnname <- asplatlon$stnname
  stnlat <- asplatlon$stnlatlon.1  # lat used in SnowMelt
  stnlon <- asplatlon$stnlatlon.2
  lstn <- length(stnname)
  aspphysionum <- unlist(readMat("../matfiles/aspphysionum.mat"))
  aspstnsel <- unlist(readMat("../matfiles/aspEcostn.mat"))

  # init vectors
#   sm_asp <- vector("list",lstn)
#   rmse_asp <- vector("numeric",lstn)
  out_asp <- vector("list",1)
  out_plot <- vector("list",1)
#   exectime <- vector("numeric",lstn)
}

# set up for reruns
#rerun <- c(6,27,29,38,43,54,71)
#rerun2 <- c(2,3,57)
#rerun3 <- 11#16#38

stnrun <- data.frame(stnname=character(),
                     rmse=numeric(),
                     exectime=numeric(),
                     stringsAsFactors=FALSE)
srct <- 0
yrct <- 0

for(ireg in 2:3){#c(1:3,5)){
  aspstnnums <- which(aspstnsel==ireg)
  for(istn in aspstnnums[1]){
    print(paste("region",ireg,"stn:",istn,stnname[istn]))
    #for(istn in 1:lstn){
    #for(istn in rerun){
  
    # ASP
    
    fn <- paste("../data/ASP/",stnname[istn],".csv",sep="")
    #  print(paste(istn,fn))
    aspallna <- read.csv(fn,skip=8)
    # read and interpret date
    if(is.element(istn,c(2,3,38,57))){
      fmt <- "%m/%d/%Y"
    }else{
      fmt <- "%Y-%m-%d"
    }
    aspdatetest <- as.Date(aspallna$Date,format=fmt)
    aspall <- aspallna[is.finite(aspdatetest),]
    aspallswe <- aspall$Snow.Water.Equivalent
    aspallprec <- aspall$Precipitation
    aspalltmin <- aspall$Temp..Min.
    aspalltmax <- aspall$Temp..Max.
    aspalldate <- as.Date(aspall$Date,format=fmt)
    aspallY <- as.numeric(format(aspalldate,"%Y"))
    aspYunique <- unique(aspallY)
    
    # check for missing and invalid data
    lognatx <- is.na(aspalltmax)
    lognatn <- is.na(aspalltmin) 
    lognaswe <- is.na(aspall$Snow.Water.Equivalent) 
    logtxlttn <- aspalltmax < aspalltmin
    
    # check for increase in SWE without corresponding daily precipitation
    sweinctol <- 0.0
    swediff <- diff(aspallswe)
    logsweinc <- c(0,swediff) > sweinctol
    logsweinc[is.na(logsweinc)]=F # if NA, this means it didn't increase and can be set to F - Q:ok?
    lognap <- is.na(aspallprec)
    logzerop <- aspallprec==0
    lognop <- lognap | logzerop
    logswenop <- logsweinc & lognop

    # create discard vector (Q: can replace na p with 0?)
    logbadtp <- lognatx | lognatn | logtxlttn | lognap #| logswenop
    logdiscard <- logbadtp | logswenop
    
    for(iyr in aspYunique){
      dateseasoni = paste(iyr-1,'-10-01',sep='')
      dateseasonf = paste(iyr,'-07-01',sep='')
      logyr <- aspalldate>=dateseasoni & aspalldate<=dateseasonf
      logvalidmodel <- logyr & !logdiscard
      logvalidcompare <- logvalidmodel & !lognaswe
      #logvalidcompareyr <- logvalidcompare[logyr]
      logvalidcompareyr <- logvalidcompare[logvalidmodel]
      swediffnop <- swediff[logyr & logswenop] # case described by Dan
      sumswenop <- sum(swediffnop[swediffnop>0],na.rm=T) 
      swediffdisc <- swediff[logyr & logdiscard] # this is really what is needed
      sumswebadtp <- sum(swediffdisc[swediffdisc>0],na.rm=T) 
      validdates <- aspalldate[logvalidmodel]
      datediff <- diff(validdates)
      if(sum(logvalidmodel)==0){
        next
      }
      
      # run and time SnowMelt
      yrct <- yrct+1
      asp <- aspall[logvalidmodel,]
      aspdate <- as.Date(asp$Date,format=fmt)
      
      out_asp[yrct] <- paste(istn,try(
        exectime <- round(system.time({
          sm_asp <- SnowMelt(Date=aspdate, precip_mm=asp$Precipitation,
                                     Tmax_C=asp$Temp..Max., Tmin_C=asp$Temp..Min., lat_deg=stnlat[istn])
          sm_asp2L <- SnowMelt2L(Date=aspdate, precip_mm=asp$Precipitation,
                               Tmax_C=asp$Temp..Max., Tmin_C=asp$Temp..Min., lat_deg=stnlat[istn])
          sm_accum <- SnowAccum(Date=aspdate, precip_mm=asp$Precipitation,
                                 Tmax_C=asp$Temp..Max., Tmin_C=asp$Temp..Min., lat_deg=stnlat[istn]
                                logNoAccumWarm=TRUE)
          sm_accumwarm <- SnowAccum(Date=aspdate, precip_mm=asp$Precipitation,
                                Tmax_C=asp$Temp..Max., Tmin_C=asp$Temp..Min., lat_deg=stnlat[istn],
                                logNoAccumWarm=FALSE)
          
        }),3)[3]
      ))
      #                     },error=function(cond){out_asp[istn]=cond},warning=function(cond){})
      aspswe <- asp$Snow.Water.Equivalent
      smswe <- sm_asp$SnowWaterEq_mm
      smswe2L <- sm_asp2L$SnowWaterEq_mm
      smsweaccum <- sm_accum$SnowWaterEq_mm
      rmse_asp <- sqrt(mean((aspswe[logvalidcompareyr]-smswe[logvalidcompareyr])^2))
      mae_asp <- mean(abs(aspswe[logvalidcompareyr]-smswe[logvalidcompareyr]))
      
      print(paste(stnname[istn],iyr,round(rmse_asp,3),exectime))

      yri <- data.frame(stnname[istn],
                        dateyr=iyr,
                        numdates=sum(logyr),
                        numbadtp=sum(logyr & logbadtp),
                        numswebadtp=sum(logyr & logdiscard),
                        sumswebadtp=sumswebadtp,
                        numvalidmodel=sum(logvalidmodel),
                        numvalidcompare=sum(logvalidcompare),
                        maxdategap=max(datediff,na.rm=T),
                        exectime=exectime,
                        mae=mae_asp,
                        rmse=rmse_asp
      )
      yrfi <- data.frame(yri,
      #                   stnname[istn],
      #                  dateyr=iyr,
      #                  numdates=sum(logyr),
                        numnatn=sum(logyr & lognatn),
                        numnatx=sum(logyr & lognatx),
                        numtxlttn=sum(logyr & logtxlttn,na.rm=T),
                        numnap=sum(logyr & lognap),
      #                  numbadtp=sum(logyr & logbadtp),
                        numswenop=sum(logyr & logswenop),
      #                  numswebadtp=sum(logyr & logdiscard),
                        sumswenop=sumswenop,
      #                  sumswebadtp=sumswebadtp,
      #                  numvalidmodel=sum(logvalidmodel),
                        numnaswe=sum(logyr & lognaswe)
      #                  numvalidcompare=sum(logvalidmodel & !lognaswe)
      )
      if(yrct==1){
        yearrun <- yri
        yearrunfull <- yrfi
      }
      else {
        yearrun[yrct,] <- yri
        yearrunfull[yrct,] <- yrfi
      }
      if(is.na(sumswebadtp)){
        print('sumswebadtp = NA')
      }
      #plot asp and snowmelt curves
      if(sumswebadtp<100){
        linew <- 3
        out_plot[yrct] <- paste(istn,try({
          #  pdf(paste("../plots/aspEco/",stnname[istn],".pdf",sep=""),width=6*3,height=6*2,pointsize=24)
          jpeg(paste("../plots/aspEco/aspEcoswebad_rmse/",sprintf("%02d",sumswebadtp),"_",round(mae_asp),"_",stnname[istn],"_",iyr,".jpg",sep=""),
               width=480*3,height=480*2,pointsize=24,quality=100)
          plot(aspdate,aspswe,col="black","l",xlab="",ylab="SWE (mm)",lwd=linew,ylim=c(0,max(aspswe,smswe,smswe2L,smsweaccum,na.rm=T)))
          lines(aspdate,smswe,col="red",lwd=linew)
          lines(aspdate,smswe2L,col="blue",lwd=linew)
          lines(aspdate,smsweaccum,col="green",lwd=linew)
          title(paste("Station",stnname[istn],iyr))#,"Exec time =",exectime[istn],"sec"))
          legend("topright",c("ASP measured","EcoH modeled","EcoH 2L modeled","Precip snow measured"),col=c("black","red","blue","green"),lwd=linew,bty="n")
          dev.off()
        }))
      }
    }
    #    next

    # check for gaps in time series > 1 day
#     aspdatediff <- as.numeric(diff(aspdate)) 
#     aspDF <- data.frame(aspdate,c(NA,aspdatediff)) # add NA to beginning of diff vector to make same length
#     aspDFgt1 <- aspDF[aspdatediff > 1,]
    
#     # run and time SnowMelt
#     asp <- aspall[is.finite(aspall$Precipitation) &
#                     is.finite(aspall$Temp..Max.) &
#                     is.finite(aspall$Temp..Min.) &
#                     is.finite(aspall$Snow.Water.Equivalent) &
#                     (aspall$Temp..Max. >= aspall$Temp..Min.),]
#     asp <- aspall[!logdiscard,]
#     aspdate <- as.Date(asp$Date,format=fmt)
# 
#     out_asp[istn] <- paste(istn,try(
#       exectime[istn] <- round(system.time(
#         sm_asp[[istn]] <- SnowMelt(Date=aspdate, precip_mm=asp$Precipitation,
#                         Tmax_C=asp$Temp..Max., Tmin_C=asp$Temp..Min., lat_deg=stnlat[istn])
#       ),3)[3]
#     ))
#     #                     },error=function(cond){out_asp[istn]=cond},warning=function(cond){})
  
#     aspswe <- asp$Snow.Water.Equivalent
#     smswe <- sm_asp[[istn]]$SnowWaterEq_mm
#     rmse_asp[istn] <- sqrt(mean((aspswe-smswe)^2,na.rm=T))
#     print(paste(stnname[istn],round(rmse_asp[istn],3),exectime[istn]))

#     srct=srct+1
#     stnrun[srct,1] = as.character(stnname[istn])
#     stnrun[srct,2:3] = data.frame(round(rmse_asp[istn],3),exectime[istn])
    
  }
}

