
get_dims = function(Data){
  NAA = Data@Misc$DataOM@Number[[1]]
  na = as.integer(dim(NAA)[2])
  nt = as.integer(dim(NAA)[3])
  nr = as.integer(dim(NAA)[4])
  nf = as.integer(dim(Data@Effort@Value)[2])
  nl = as.integer(length(Data@Misc$DataOM@OM@Stock[[1]]@Length@Classes))
  data.frame(nt=nt,na=na,nf=nf,nr=nr,nl=nl)
}

writedebug=function(debugfile, debugts, Data, Sim_Ind, debugDatafile){
  if(!is.na(debugDatafile))saveRDS(Data,debugDatafile)

  if(!is.na(debugfile)){

    if(max(Data@Years)==debugts){
      DebugData = list(Data = Data, Sim_Ind = Sim_Ind)
      saveRDS(DebugData, debugfile)
      cat(paste0("Debug data in time step ",debugts, " outputted to: ",debugfile, "\n"))
      stop()
    }
  }

}

mahiMPcheck = function(docheck, Data, TAC, Effort, Smin, Smax, BL, verbose=F){
  if(docheck){
    d = get_dims(Data); nt = d$nt; na = d$na; nf = d$nf; nr = d$nr

    if(!is.na(Effort[1])){
      if(!is.matrix(Effort)){
        if(!(length(Effort)%in%c(1,nf))) stop(paste0("Effort must be a vector of length 1 or ", nf))
      }else{
        if(!(nrow(Effort)==nf & ncol(Effort)==nr))stop(paste0("Effort matrices should be ",nf, " (fleets) rows by ", nr, " (areas) columns"))
      }
    }

    if(!is.na(TAC[1])){
      if(!is.matrix(TAC)){
        if(!(length(TAC)%in%c(1,nf))) stop(paste0("TAC must be a vector of length 1 or ", nf))
      }else{
        if(!(nrow(TAC)==nf & ncol(TAC)==nr))stop(paste0("TAC matrices should be ",nf, " (fleets) rows by ", nr, " (areas) columns"))
      }
    }

    if(!is.na(Smin[1]) & !length(Smin)%in%c(1,nf)) stop(paste0("Smin must be either a single number or a vector ",nr, " long"))
    if(!is.na(Smax[1]) & !length(Smax)%in%c(1,nf)) stop(paste0("Smax must be either a single number or a vector ",nr, " long"))
    if(!is.na(BL[1]) & length(BL)!=nf) stop(paste0("BL must be a vector ",nf, " long"))
    #cat("mahiMPcheck CLEARED \n")
  }
}

do_default = function(rel_TAE, TAC, rel_TAC, TAE, HCR_use){
  if(!is.na(TAC[1]) | !is.na(rel_TAC[1]) | !is.na(TAE[1]) | HCR_use) rel_TAE = NA
  rel_TAE
}

do_TAE = function(Data, TAE, rel_TAE, season, ad, nyrs=3, ns=4){
  ad@EffType = 'Abs'
  if(!is.na(TAE[1])){  # Absolute effort (can be total, by fleet, by fleet x area)
   ad@Effort = TAE
  }
  if(!is.na(rel_TAE[1])) {
    takets = match(Data@YearLH, Data@Years) - (nyrs * ns -1):0
    Effort = Data@Misc$DataOM@Effort[1,takets,]
    nf = dim(Effort)[2]
    muEffort = apply(array(Effort,c(ns,nyrs,nf)),c(1,3),mean)

    #refyrs = match(Data@YearLH, Data@Years)-3:0
    #ts = refyrs[season]  #  pick from last four historical seasons
    histEff = muEffort[season,] #Data@Misc$DataOM@Effort[1,ts,] #Data@Effort@Value[ts,]
    newEffort <- histEff*rel_TAE

    Dist = Data@Misc$DataOM@Distribution[1,takets,,]
    nr = dim(Dist)[3]
    muDist = apply(array(Dist,c(ns,nyrs,nf,nr)),c(1,3,4),mean)
    #EffortFleetArea <- newEffort* Data@Misc$DataOM@Distribution[1,ts,,]
    EffortFleetArea <- newEffort * muDist[season,,]
    ad@Effort = EffortFleetArea
  }
  ad
}

do_TAC = function(Data, TAC, rel_TAC, season, ad){
  unlimitedTAC = 1E10
  ad@TACType = "Removals"
  if(any(!is.na(TAC))){ # works with both total TAC and TAC by fleet
    sTAC =  aperm(rTAC,c(2,1,3))
    TACfrac = sTAC/apply(sTAC,1,sum)
    tempTAC = apply(TAC * TACfrac[,season,],1,sum) # fleet specific TAC summed over areas
    tempTAC[is.na(tempTAC)] = unlimitedTAC
    ad@TAC
  }
  if(!is.na(rel_TAC[1])){
    ad@TAC = rel_TAC * apply(rTAC[season,,],1,sum) # relTAC can be total, by fleet, by fleet x season
  }
  ad
}

do_PRM = function(Data, PRM, ad){
  LC = Data@Misc$DataOM@OM@Stock$Dolphinfish@Length@Classes
  ds = get_dims(Data)
  LC = Data@Misc$DataOM@OM@Stock$Dolphinfish@Length@Classes
  mal = rep(NA,ds$nl)
  dmlist = list()
  for(ff in 1:ds$nf){
    mal[]=PRM[[ff]]
    dmlist[[ff]] = DiscardMortality(MeanAtLength = mal,Classes = LC)#array(PRM[[ff]],c(nLC,dm$nt,dm$nr))
  }
  ad@DiscardMortality = dmlist
  ad
}

do_size = function(Data, SLmin, SLmax, ad){
  ds = get_dims(Data)
  LC = Data@Misc$DataOM@OM@Stock$Dolphinfish@Length@Classes
  retmat = matrix(1, ds$nl, ds$nr)
  retlist = list()
  for(ff in 1:ds$nf){
    retvec = c(0,1)[as.integer(LC > SLmin[ff] & LC < SLmax[ff])+1]
    retmat[]=retvec
    retlist[[ff]] = Retention(MeanAtLength = retmat , Classes = LC)
  }
  ad@Retention = retlist
  ad
}





