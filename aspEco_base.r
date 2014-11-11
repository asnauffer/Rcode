# reads ASP data from a csv file
# runs SnowMelt from EcoHydRology using T & P data

if(R.version[1]=="x86_64-pc-linux-gnu") {
  library("EcoHydRology", lib.loc="/net/home/asnauffer/R/x86_64-pc-linux-gnu-library/3.1")
  library("R.matlab", lib.loc="/net/home/asnauffer/R/x86_64-pc-linux-gnu-library/3.1")
  setwd("/net/home/asnauffer/PhD/R")
} else {
  library("EcoHydRology", lib.loc="C:/Users/drew/Documents/R/win-library/3.1")
  library("R.matlab", lib.loc="C:/Users/drew/Documents/R/win-library/3.1")
  setwd("C:/Users/Drew/Documents/PhD/R")
}

skipreload <- FALSE
if(skipreload){
  print("skipping data load & variable init")
} else {
  print("loading and initializing...")
  # load tmax tmin and precip values for all stns
#   tmax <- data.frame(readMat("../matfiles/MERRAtmaxatASP.mat"))-273.16 # convert Kelvin to Celsius
#   tmin <- data.frame(readMat("../matfiles/MERRAtminatASP.mat"))-273.16 # convert Kelvin to Celsius
#   prec <- data.frame(readMat("../matfiles/MERRAprecatASP.mat"))*3600*24 # convert kg m^2 s^1 to mm d^-1
#   gpdate <- readMat("../matfiles/MERRAgdate.mat") 
  #tmax  <- data.frame(readMat("Drew/PhD/Research/matfiles/ERAtmaxatASP.mat"))-273.16 # convert Kelvin to Celsius
  #tmin <- data.frame(readMat("Drew/PhD/Research/matfiles/ERAtminatASP.mat"))-273.16 # convert Kelvin to Celsius
  #prec <- data.frame(readMat("Drew/PhD/Research/matfiles/ERAprecatASP.mat"))*3600*24 # convert kg m^2 s^1 to mm d^-1
  #gpdate <- readMat("Drew/PhD/Research/matfiles/ERAgdate.mat") 
  #gdate <- data.frame(gpdate)[,1] # cast to factor
  #gpdatefactor <- data.frame(gpdate)[,1] # cast to factor
#   gdate <- as.Date(gpdate$MERRAgdate)
#   gdatemo <- as.numeric(format(as.Date(gdate),"%m"))
#   gdateda <- as.numeric(format(as.Date(gdate),"%d"))
#   gdateyr <- as.numeric(format(as.Date(gdate),"%Y"))
#   gdateYmd <- format(as.Date(gdate),"%Y-%m-%d")
  asplatlon <- data.frame(readMat("../matfiles/ASPalllatlon.mat"))
  stnname <- asplatlon$stnname
  stnlat <- asplatlon$stnlatlon.1  # lat used in SnowMelt
  stnlon <- asplatlon$stnlatlon.2
#   d2 <- stnlat^2 + stnlon^2
#   stnsort <- sort(d2,decreasing=TRUE,index.return=TRUE)
  lstn <- length(stnname)
  aspphysionum <- unlist(readMat("../matfiles/aspphysionum.mat"))
  aspstnsel <- unlist(readMat("../matfiles/aspEcostn.mat"))

  # read climate norm file for ASPs only
  # fncn <- "../data/ASPlatlonelev_ClimateBC_Normal_1981_2010AM.csv"
  # cnorm <- read.csv(fncn)
  # cnorm[cnorm==-9999] <- NA
  
  # init vectors
  sm_asp <- vector("list",lstn)
#   sm_era <- vector("list",lstn)
#   sm_prec_corr <- vector("list",lstn)
#   sm_temp_corr <- vector("list",lstn)
#   sm_ttp_corr <- vector("list",lstn)
  rmse_asp <- vector("numeric",lstn)
#   rmse_gp <- vector("numeric",lstn)
#   rmse_gp_pcorr <- vector("numeric",lstn)
#   rmse_gp_tcorr <- vector("numeric",lstn)
#   rmse_gp_ttp <- vector("numeric",lstn)
  out_asp <- vector("list",lstn)
  out_plot <- vector("list",lstn)
  exectime <- vector("numeric",lstn)
#   out_era <- vector("list",lstn)
#   out_erap <- vector("list",lstn)
#   out_erat <- vector("list",lstn)
#   out_erattp <- vector("list",lstn)
}

# set up for reruns
#rerun <- c(6,27,29,38)
#rerun <- c(6,27,29,38,43,54,71)
#rerun2 <- c(2,3,57)
#rerun3 <- 11#16#38

stnrun <- data.frame(stnname=character(),
                     rmse=numeric(),
                     exectime=numeric(),
                     stringsAsFactors=FALSE)
srct <- 0
yearrun <- data.frame(stnname=character(),
                      rmse=numeric(),
                      exectime=numeric(),
                      stringsAsFactors=FALSE)
yrct <- 0

for(ireg in c(1:3,5)){
  aspstnnums <- which(aspstnsel==ireg)
  for(istn in aspstnnums){
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
    aspallY <- format(aspalldate,"%Y")
    aspYunique <- unique(aspallY)
    #print(aspYunique)
    
    # check for missing and invalid data
    lognatx <- is.na(aspalltmax)
    lognatn <- is.na(aspalltmin) 
    lognaswe <- is.na(aspall$Snow.Water.Equivalent) 
    loginvt <- aspalltmax < aspalltmin
    
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
    logbadtp <- lognatx | lognatn | loginvt | lognap #| logswenop
    logdiscard <- logbadtp | logswenop
    
    for(iyr in aspYunique){
      dateseasoni = paste(as.numeric(iyr)-1,'-10-01',sep='')
      dateseasonf = paste(iyr,'-07-01',sep='')
      logyr <- aspalldate>=dateseasoni & aspalldate<=dateseasonf
      sumswenop <- sum(swediff[logyr & logswenop]) # case described by Dan
      sumswebadtp <- sum(swediff[logyr & logdiscard]) # this is really what is needed
      validdatesi <- aspalldate[logyr & !logdiscard]
      if(length(validdatesi)>0){
        print(paste(dateseasoni,'to',dateseasonf))
      }
      yrct <- yrct+1
      yearrun[yrct,1] <- as.character(stnname[istn])
#       yearrun <- data.frame(numdates=sum(logyr),
#                             numvalidtmin=sum(logyr&!lognatn),
#                             numvalidtmax=sum(logyr&!lognatx),
#                             numtxgttn=sum(logyr&!loginvt[is.finite(loginvt)]),
#                             numvalidprec=sum(logyr&!lognap),
#                             numvalidswe=sum(logyr&lognaswe)
#                             )
      yearrun <- data.frame(numdates=sum(logyr),
                            numnatn=sum(logyr&lognatn),
                            numnatx=sum(logyr&lognatx),
                            numinvt=sum(logyr&loginvt[is.finite(loginvt)]),
                            numnap=sum(logyr&lognap),
                            numswenop=sum(logyr&logswenop),
                            numnaswe=sum(logyr&lognaswe)
      )
    }
#    next
# 
#     logswenapi <- logical(length=nrow(aspall))
#     logswenopi <- logical(length=nrow(aspall))
#     logswenoti <- logical(length=nrow(aspall))
#     for(iobs in c(2:length(aspsweall))){
#       # check for valid SWE before and on day
#       if(is.finite(aspsweall[iobs]) & is.finite(aspsweall[iobs-1])){
#         # check for SWE difference
#         if(aspsweall[iobs]-aspsweall[iobs-1] > sweinctol){
#           # check for no precips or temps
#           if(is.na(aspprecall[iobs])){
#             #print(paste(c(istn,aspalldate[iobs],aspprecall[iobs],aspsweall[iobs],aspsweall[iobs-1])))
#             logswenapi[iobs] <- TRUE
#           }
#           else if((aspprecall[iobs]==0)) {
#             #print(paste(c(istn,aspalldate[iobs],aspprecall[iobs],aspsweall[iobs],aspsweall[iobs-1])))
#             logswenopi[iobs] <- TRUE
#           }
#           # check for no temps
#           if(is.na(asptmaxall[iobs]) | is.na(asptminall[iobs])) {
#             logswenoti[iobs] <- TRUE
#           }
#         }
#       }
#     }

    # check for gaps in time series > 1 day
#     aspdatediff <- as.numeric(diff(aspdate)) 
#     aspDF <- data.frame(aspdate,c(NA,aspdatediff)) # add NA to beginning of diff vector to make same length
#     aspDFgt1 <- aspDF[aspdatediff > 1,]
    
    # run and time SnowMelt
    #    MRdate <- data.frame(dateasp)
    asp <- aspall[is.finite(aspall$Precipitation) &
                    is.finite(aspall$Temp..Max.) &
                    is.finite(aspall$Temp..Min.) &
                    is.finite(aspall$Snow.Water.Equivalent) &
                    (aspall$Temp..Max. >= aspall$Temp..Min.),]
    asp <- aspall[!logdiscard,]
    aspdate <- as.Date(asp$Date,format=fmt)

    out_asp[istn] <- paste(istn,try(
      exectime[istn] <- round(system.time(
        sm_asp[[istn]] <- SnowMelt(Date=aspdate, precip_mm=asp$Precipitation,
                        Tmax_C=asp$Temp..Max., Tmin_C=asp$Temp..Min., lat_deg=stnlat[istn])
      ),3)[3]
    ))
    #                     },error=function(cond){out_asp[istn]=cond},warning=function(cond){})
  
    aspswe <- asp$Snow.Water.Equivalent
    smasp <- sm_asp[[istn]]$SnowWaterEq_mm
    rmse_asp[istn] <- sqrt(mean((aspswe-smasp)^2,na.rm=T))
    print(paste(stnname[istn],round(rmse_asp[istn],3),exectime[istn]))
    srct=srct+1
    stnrun[srct,1] = as.character(stnname[istn])
    stnrun[srct,2:3] = data.frame(round(rmse_asp[istn],3),exectime[istn])

    #plot asp and snowmelt curves
    linew <- 3
  out_plot[istn] <- paste(istn,try({
  #  pdf(paste("../plots/aspEco/",stnname[istn],".pdf",sep=""),width=6*3,height=6*2,pointsize=24)
    jpeg(paste("../plots/aspEco/aspEcoPchk/",aspphysionum[istn],"_",nrow(asp),"_",stnname[istn],".jpg",sep=""),
         width=480*3,height=480*2,pointsize=24,quality=100)
    plot(aspdate,aspswe,col="black","l",xlab="",ylab="SWE (mm)",lwd=linew,ylim=c(0,max(aspswe,smasp)))
    lines(aspdate,smasp,col="red",lwd=linew)
    title(paste("Station",stnname[istn]))#,"Exec time =",exectime[istn],"sec"))
    legend("topright",c("ASP measured","EcoH modeled"),col=c("black","red"),lwd=linew,bty="n")
    dev.off()
  }))
  
  #writeMat("../matfiles/exectime.mat",exectime=exectime)
  
  }
}

#     # ERA
#   
#     tmax_momean <- numeric(12)
#     tmax_modiff <- numeric(12)
#     tmax_corr <- tmax
#     tmax_corr[,] <- NA
#     tmax_climnorm <- cnorm[istn,6:17]
#     tmin_momean <- numeric(12)
#     tmin_modiff <- numeric(12)
#     tmin_corr <- tmin
#     tmin_corr[,] <- NA
#     tmin_climnorm <- cnorm[istn,18:29]
#     prec_momean <- numeric(12)
#     prec_moratio <- numeric(12)
#     prec_corr <- prec
#     prec_corr[,] <- NA
#     prec_climnorm <- cnorm[istn,42:53]
#     
#     for(imo in 1:12) {
#       # correct precip using ratio
#       logmoi <- gdatemo==imo
#       numyrsi <- length(unique(gdateyr[logmoi]))
#       prec_momean[imo] <- sum(prec[logmoi,istn])/numyrsi
#       prec_moratio[imo] <- as.numeric(prec_momean[imo]/prec_climnorm[imo])
#       prec_corr[logmoi,istn] <- prec[logmoi,istn]/prec_moratio[imo]
#       
#       # correct tmax & tmin additively
#       tmax_momean[imo] <- mean(tmax[logmoi,istn])
#       tmax_modiff[imo] <- as.numeric(tmax_momean[imo]-tmax_climnorm[imo])
#       tmax_corr[logmoi,istn] <- tmax[logmoi,istn]-tmax_modiff[imo]
#       tmin_momean[imo] <- mean(tmin[logmoi,istn])
#       tmin_modiff[imo] <- as.numeric(tmin_momean[imo]-tmin_climnorm[imo])
#       tmin_corr[logmoi,istn] <- tmin[logmoi,istn]-tmin_modiff[imo]
#       
#         #     logTFd <- tmin[logmoi,istn] > tmax[logmoi,istn] 
#         #     logTcorrFd <- tmin_corr[logmoi,istn] > tmax_corr[logmoi,istn] 
#         #     logTxcorrFd <- tmin[logmoi,istn] > tmax_corr[logmoi,istn] 
#         #     logTncorrFd <- tmin_corr[logmoi,istn] > tmax[logmoi,istn] 
#         #     
#         #     print(c(istn,imo,sum(logTFd),sum(logTcorrFd),sum(logTxcorrFd),sum(logTncorrFd),-9999))
#         #     
#         #     if(sum(logTcorrFd)>0){
#         #       tmini <- tmin[logmoi,istn]
#         #       tmaxi <- tmax[logmoi,istn]
#         #       tminic <- tmin_corr[logmoi,istn]
#         #       tmaxic <- tmax_corr[logmoi,istn]
#         #       tnd <- tminic[logTcorrFd]-tmini[logTcorrFd]
#         #       tnc <- tminic[logTcorrFd]
#         #       txc <- tmaxic[logTcorrFd]
#         #       txd <- tmaxic[logTcorrFd]-tmaxi[logTcorrFd]
#         #       #    print(data.frame(tmini[logTcorrFd],tminic[logTcorrFd],tmaxi[logTcorrFd],tmaxic[logTcorrFd]))
#         #       print(data.frame(istn,imo,tnd,tnc,txc,txd))
#         #     }
#   
#     }
#     
#     # remove lines with corrected Tmin>Tmax (EcoH SnowMelt cannot handle)
#     logTcorrbad <- tmin_corr[,istn] > tmax_corr[,istn] 
#     gdate_valid <- gdate[!logTcorrbad]
#     gdateYmd_valid <- gdateYmd[!logTcorrbad]
#     prec_valid <- prec[!logTcorrbad,istn]
#     tmax_valid <- tmax[!logTcorrbad,istn]
#     tmin_valid <- tmin[!logTcorrbad,istn]
#     prec_corr_valid <- prec_corr[!logTcorrbad,istn]
#     tmax_corr_valid <- tmax_corr[!logTcorrbad,istn]
#     tmin_corr_valid <- tmin_corr[!logTcorrbad,istn]
#     
#     # run EcoHydRology SnowMelt model
#     print(paste('running snowmelt at stn',istn))
#   out_era[istn] <- paste(istn,try(
#     sm_era[[istn]] <- SnowMelt(Date=gdate_valid, precip_mm=prec_valid,
#               Tmax_C=tmax_valid, Tmin_C=tmin_valid, lat_deg=stnlat[istn])
#     ))
#   out_erap[istn] <- paste(istn,try(
#     sm_prec_corr[[istn]] <- SnowMelt(Date=gdate_valid, precip_mm=prec_corr_valid,
#               Tmax_C=tmax_valid, Tmin_C=tmin_valid, lat_deg=stnlat[istn])
#     ))
#   out_erat[istn] <- paste(istn,try(
#     sm_temp_corr[[istn]] <- SnowMelt(Date=gdate_valid, precip_mm=prec_valid,
#               Tmax_C=tmax_corr_valid, Tmin_C=tmin_corr_valid, lat_deg=stnlat[istn])
#   ))
#   out_erattp[istn] <- paste(istn,try(
#     sm_ttp_corr[[istn]] <- SnowMelt(Date=gdate_valid, precip_mm=prec_corr_valid,
#               Tmax_C=tmax_corr_valid, Tmin_C=tmin_corr_valid, lat_deg=stnlat[istn])
#   ))
#   
#     # extract calculated SWE
#     smera <- sm_era[[istn]]$SnowWaterEq_mm
#     smerap <- sm_prec_corr[[istn]]$SnowWaterEq_mm
#     smerat <- sm_temp_corr[[istn]]$SnowWaterEq_mm
#     smerattp <- sm_ttp_corr[[istn]]$SnowWaterEq_mm
#   
#     # find RMSE wrt ASP measured values
#     logaspdate <- is.element(aspdate,gdate_valid)
#     logeradate <- is.element(gdate_valid,aspdate)
#     errasp <- aspswe[logaspdate]-smasp[logaspdate]
#     rmse_asp[istn] <- sqrt(mean((errasp[is.finite(errasp)])^2))
#     errgp <- aspswe[logaspdate]-smera[logeradate]
#     rmse_gp[istn] <- sqrt(mean((errgp[is.finite(errgp)])^2))
#     errgpp <- aspswe[logaspdate]-smerap[logeradate]
#     rmse_gp_pcorr[istn] <- sqrt(mean((errgpp[is.finite(errgpp)])^2))
#     errgpt <- aspswe[logaspdate]-smerat[logeradate]
#     rmse_gp_tcorr[istn] <- sqrt(mean((errgpt[is.finite(errgpt)])^2))
#     errgpttp <- aspswe[logaspdate]-smerattp[logeradate]
#     rmse_gp_ttp[istn] <- sqrt(mean((errgpttp[is.finite(errgpttp)])^2))
# 
#     print(c(istn,rmse_asp[istn],rmse_gp[istn],rmse_gp_pcorr[istn],rmse_gp_tcorr[istn],rmse_gp_ttp[istn]))
#   } #end stn for

# # write RMSEs sorted by asp stn NW -> SE
# rmse <- data.frame(rmse_asp,rmse_gp,rmse_gp_pcorr,rmse_gp_tcorr,rmse_gp_ttp)[stnsort$ix,]
# aspstnID <- stnname[stnsort$ix]
# writeMat("../matfiles/MERRArmses.mat",MERRArmse=rmse,stnsorti=stnsort$ix,aspstnID=aspstnID)
# 
# # write SWE & SnowModel outputs
# writeMat("../matfiles/ASPswe.mat",aspswe=aspswe)
# writeMat("../matfiles/ASPsm.mat",aspsm=smasp)
# writeMat("../matfiles/aspdate.mat",aspdate=aspdate)
# writeMat("../matfiles/MERRAsm.mat",MERRAsm=smera)
# writeMat("../matfiles/MERRAsmp.mat",MERRAsm=smerap)
# writeMat("../matfiles/MERRAsmt.mat",MERRAsm=smerat)
# writeMat("../matfiles/MERRAsmttp.mat",MERRAsm=smerattp)
# writeMat("../matfiles/gpdate.mat",gpdate=gdate_valid)
# 
#writeMat("../matfiles/MRMERRA.mat",MERRAsm=smera,MERRAsmttp=smerattp)
#smasp=smasp,aspswe=aspswe,MRdate=MRdate,MRPasp=asp$Precipitation,MRTxasp=asp$Temp..Max.,MRTnasp=asp$Temp..Min.)
# conout <- "Drew/PhD/Research/matfiles/KRerasmTTPadj.mat"
# writeMat(conout,gpdate=gdateYmd_valid,MRsmera=MRsmera,MRsmerap=MRsmerap,MRsmerat=MRsmerat,MRsmerattp=MRsmerattp)
#,MRPera=prec[,istn],MRTxera=tmax[,istn],MRTnera=tmin[,istn])
#MRsmerad=MRsmerad,MRsmerax=MRsmerax,MRsmeran=MRsmeran,

