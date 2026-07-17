

weightstats = function(x, dms, MPnames){

  CAA = PropCAA = TotLen = array(0, c(dms$nsim, dms$na, dms$ny, dms$nf, dms$nMP))
  dimnames(CAA)=dimnames(PropCAA)=dimnames(TotLen)=list(paste("sim",1:dms$nsim), paste("age",1:dms$na), Time_steps, Fleets, MPnames)
  LAA = x@OM@Stock$Dolphinfish@Length@MeanAtAge # Sim, age, all ts
  CAAp = apply(x@LandingsAtAge$Dolphinfish,c(1,2,3,4,6),sum)  # Sim, age, proj ts, fleet, MP
  CAAh = apply(x@Hist@LandingsAtAge$Dolphinfish,1:4,sum)    # Sim, age, hist ts, fleet
  hind = MSEtool:::TEG(c(dms$nsim, dms$na, dms$nt, dms$nf, dms$nMP))
  CAA[hind] = CAAh[hind[,1:4]]
  pind = MSEtool:::TEG(c(dms$nsim, dms$na, dms$np, dms$nf, dms$nMP))
  pind[,3] = pind[,3]+dms$nt
  CAA[pind] = CAAp
  Cind = MSEtool:::TEG(dim(CAA))
  TotLen = CAA * LAA[Cind[,1:3]]
  PropCAA = CAA / apply(CAA, c(1,3,4,5),sum,na.rm=T)[Cind[,c(1,3,4,5)]]
  MuLen = apply(TotLen,c(1,3,4,5),sum,na.rm=T) / apply(CAA,c(1,3,4,5),sum,na.rm=T)
  list(MuLen=MuLen, PropCAA=PropCAA)
}

getVB = function(x, dms, MPnames, doret=T){

  NAA = TotWt = array(0, c(dms$nsim, dms$na, dms$ny, dms$nf, dms$nMP))
  dimnames(NAA)=dimnames(TotWt)=list(paste("sim",1:dms$nsim), paste("age",1:dms$na), Time_steps, Fleets, MPnames)
  WAA = x@OM@Stock$Dolphinfish@Weight@MeanAtAge # Sim, age, all ts

  NAAp = apply(x@Number$Dolphinfish,c(1,2,3,5),sum)  # Sim, age, proj ts, MP
  NAAh = apply(x@Hist@Number$Dolphinfish,1:3,sum)    # Sim, age, hist ts
  hind = MSEtool:::TEG(c(dms$nsim, dms$na, dms$nt, dms$nf, dms$nMP))
  NAA[hind] = NAAh[hind[,1:3]]
  pind = nind = MSEtool:::TEG(c(dms$nsim, dms$na, dms$np, dms$nf, dms$nMP))
  nind[,3] = pind[,3]+dms$nt
  NAA[nind] = NAAp[pind[,c(1,2,3,5)]]
  Nind = MSEtool:::TEG(dim(NAA))
  TotWt = NAA * WAA[Nind[,1:3]]

  sel = array(unlist(lapply(x@OM@Fleet$Dolphinfish, function(y)y@Selectivity@MeanAtAge[,,,1])),c(dms$nsim,dms$na,dms$ny,dms$nf,dms$nMP)) # area 1 is fine
  retsel = sel

  if(doret){
    Retention  =  x@Misc$Retention
    for(mp in 1:dms$nMP){
      for(ff in 1:dms$nf){
        mpmatch = match(MPnames[mp],names(Retention))
        ret =Retention[[mpmatch]][[1]][[ff]]
        if(!is.null(ret)){
          retsel[,,dms$nt + 1:dms$np,ff,mp] = retsel[,,dms$nt + 1:dms$np,ff,mp] * ret$MeanAtAge[,,,1] # area 1 is fine
        }
      }
    }
  }

  apply(TotWt * retsel, c(1,3,4,5),sum)

}


mahiResults = function(MSElist, stock = 1, ns = 4, doMSY = T){

  if(class(MSElist)!='list')MSElist=list(MSElist)
  ylab = unique(floor(Time_steps))
  nylab = length(ylab)
  cat("Processing MSE Results \n")
  OMnames = names(MSElist)

  dms = getMSEdims(MSElist)                                             # Dimensions
  cat(paste0("Dimensions: ",paste(paste0(names(dms)," = ",dms),collapse=", ")," \n"))

  SSB = B = array(NA, c(dms$nsim, dms$ny, dms$nMP, dms$nMSE))           # sim, allyear, MP, OM,
  dimnames(SSB)=list(paste("sim",1:dms$nsim),Time_steps,names(MSElist[[1]]@MPs),OMnames)
  Landings = array(NA, c(dms$nsim, dms$ny, dms$nf, dms$nMP, dms$nMSE))  # sim, allyear, fleet, MP, OM
  dimnames(Landings) = list(paste("sim",1:dms$nsim), Time_steps, Fleets, names(MSElist[[1]]@MPs), OMnames)

  cat("Extracting Spawning Stock Biomass (SSB) \n")
  # --- SSB  ----------------------------------------------------------------------------------------
  SSBa = lapply(MSElist, function(x)x@SBiomass[,stock,,])    # projected
  SSB[,dms$nt+1:dms$np,,] = array(unlist(SSBa),c(dms$nsim, dms$np, dms$nMP, dms$nMSE))
  SSBh = lapply(MSElist,function(x)x@Hist@SBiomass[,stock,]) # historical
  for(mm in 1:dms$nMP) SSB[,1:dms$nt,mm,] = array(unlist(SSBh),c(dms$nsim, dms$nt, dms$nMSE))
  SSB = aperm(SSB,c(4,1,2,3))

  cat("Extracting Biomass (B) \n")
    # --- B -------------------------------------------------------------------------------------------
  Ba = lapply(MSElist, function(x)x@Biomass[,stock,,])    # projected
  B[,dms$nt+1:dms$np,,] = array(unlist(Ba),c(dms$nsim, dms$np, dms$nMP, dms$nMSE))
  Bh = lapply(MSElist,function(x)x@Hist@Biomass[,stock,]) # historical
  for(mm in 1:dms$nMP) B[,1:dms$nt,mm,] = array(unlist(Bh),c(dms$nsim, dms$nt, dms$nMSE))
  B = aperm(B,c(4,1,2,3))

  cat("Extracting Landings (Landings) \n")
  # --- Landings ------------------------------------------------------------------------------------
  Landa = lapply(MSElist, function(x)x@Landings[,stock,,,])    # projected
  Landings[,dms$nt+1:dms$np,,,] = array(unlist(Landa),c(dms$nsim, dms$np, dms$nf, dms$nMP, dms$nMSE))
  Landh = lapply(MSElist,function(x)x@Hist@Landings[,stock,,]) # historical
  for(mm in 1:dms$nMP) Landings[,1:dms$nt,,mm,] = array(unlist(Landh),c(dms$nsim, dms$nt, dms$nf, dms$nMSE))
  Landings_byfleet_s = aperm(Landings,c(5,1,2,3,4)) # OM first
  Landings_byfleet = apply(array(Landings_byfleet_s,c(dms$nMSE,dms$nsim,ns,dms$ny/ns,dms$nf,dms$nMP)),c(1,2,4,5,6),sum)
  dimnames(Landings_byfleet) = list(OMnames,paste("sim",1:dms$nsim),unique(floor(Time_steps)),Fleets,names(MSElist[[1]]@MPs))
  Landings = apply(Landings, c(5,1,2,4),sum) # summed over fleets, OM first

  cat("Calculating Catch Diff (Cdif) \n")
  # --- Catch diff ----------------------------------------------------------------------------------
  t1= 1:(dim(Landings_byfleet)[3]-1)
  t2= 2:dim(Landings_byfleet)[3]
  CDif = abs((Landings_byfleet[,,t2,,,drop=F]- Landings_byfleet[,,t1,,,drop=F])/ Landings_byfleet[,,t1,,,drop=F])
  dimnames(CDif) = list(OMnames,paste("sim",1:dms$nsim),ylab[2:nylab],Fleets,names(MSElist[[1]]@MPs))

  cat("Calculating Proportion of Landings in each age class Biomass (PropCAA) \n")
  # --- Size of catch -------------------------------------------------------------------------------
  szlist = lapply(MSElist, weightstats, dms=dms, MPnames = names(MSElist[[1]]@MPs)) # MuLen [sim, ts, fleet, mp] and PropCAA [sim, age, ts, fleet mp]
  MuLen = array(NA, c(dms$nsim, dms$ny, dms$nf, dms$nMP, dms$nMSE))
  #dimnames(MuLen) = list(paste("sim",1:dms$nsim), Time_steps, Fleets, names(MSElist[[1]]@MPs), paste0("OM_",1:dms$nMSE))
  PropCAA = array(NA, c(dms$nsim, dms$na, dms$ny, dms$nf, dms$nMP, dms$nMSE))
  dimnames(PropCAA) = list(paste("sim",1:dms$nsim), paste("age",1:dms$na), Time_steps, Fleets, names(MSElist[[1]]@MPs), OMnames)

  for(mm in 1:dms$nMSE){
    MuLen[,,,,mm] = szlist[[mm]]$MuLen
    PropCAA[,,,,,mm] = szlist[[mm]]$PropCAA
  }

  MuLen = aperm(MuLen,c(5,1,2,3,4))
  PropCAA=aperm(PropCAA, c(6,1,2,3,4,5))
  MuLen = apply(array(MuLen,c(dms$nMSE,dms$nsim,ns,dms$ny/ns,dms$nf,dms$nMP)),c(1,2,4,5,6),sum)
  dimnames(MuLen) = dimnames(Landings_byfleet)
  yind = match(MSElist[[1]]@OM@CurrentYear,unique(floor(Time_steps)))

  refMuLen = array(MuLen[,,yind,,],c(dms$nMSE,dms$nsim,dms$nf,dms$nMP))
  MLind = MSEtool:::TEG(dim(MuLen))
  MuLen[MLind] = MuLen[MLind]/refMuLen[MLind[,c(1,2,4,5)]]

  cat("Calculating Vulnerable Biomass as a retained catch rate predictor (VBr) \n")
  # --- Retained Catch Rate ---------------------------------------------------------------------------------
  VBrlist = lapply(MSElist, getVB, dms=dms, MPnames = names(MSElist[[1]]@MPs)) # MuLen [sim, ts, fleet, mp] and PropCAA [sim, age, ts, fleet mp]
  VBr = array(unlist(VBrlist),c(dms$nsim, dms$ny, dms$nf,dms$nMP,dms$nMSE))
  #dimnames(VB) = list(paste("sim",1:dms$nsim),Time_steps, Fleets, names(MSElist[[1]]@MPs), paste0("OM_",1:dms$nMSE))
  VBr = aperm(VBr,c(5,1,2,3,4))
  VBr = apply(array(VBr,c(dms$nMSE,dms$nsim,ns,dms$ny/ns,dms$nf,dms$nMP)),c(1,2,4,5,6),sum)
  dimnames(VBr) = dimnames(Landings_byfleet)
  yind = match(MSElist[[1]]@OM@CurrentYear,unique(floor(Time_steps)))
  refVBr = array(VBr[,,yind,,],c(dms$nMSE,dms$nsim,dms$nf,dms$nMP))
  VBrind = MSEtool:::TEG(dim(VBr))
  VBr[VBrind] = VBr[VBrind]/refVBr[VBrind[,c(1,2,4,5)]]
  dimnames(VBr) = dimnames(Landings_byfleet)

  cat("Calculating Vulnerable Biomass as a catch rate predictor (VB) \n")
  # --- Retained Catch Rate ---------------------------------------------------------------------------------
  VBlist = lapply(MSElist, getVB, dms=dms, MPnames = names(MSElist[[1]]@MPs), doret=F) # MuLen [sim, ts, fleet, mp] and PropCAA [sim, age, ts, fleet mp]
  VB = array(unlist(VBlist),c(dms$nsim, dms$ny, dms$nf,dms$nMP,dms$nMSE))
  #dimnames(VB) = list(paste("sim",1:dms$nsim),Time_steps, Fleets, names(MSElist[[1]]@MPs), paste0("OM_",1:dms$nMSE))
  VB = aperm(VB,c(5,1,2,3,4))
  VB = apply(array(VB,c(dms$nMSE,dms$nsim,ns,dms$ny/ns,dms$nf,dms$nMP)),c(1,2,4,5,6),sum)
  dimnames(VB) = dimnames(Landings_byfleet)
  yind = match(MSElist[[1]]@OM@CurrentYear,unique(floor(Time_steps)))
  refVB = array(VB[,,yind,,],c(dms$nMSE,dms$nsim,dms$nf,dms$nMP))
  VBind = MSEtool:::TEG(dim(VB))
  VB[VBind] = VB[VBind]/refVB[VBind[,c(1,2,4,5)]]
  dimnames(VB) = dimnames(Landings_byfleet)

  # --- MSY calcs ----------------------------------------------------------------------------------
  Brel = Frel = NULL

  if(doMSY){ # do MSY calcs using Martell and Walters approach
    msy_out = lapply(MSElist, calcmahi_MSY) # calculate MSY quantities
    cat("Calculating MSY reference points (Brel, Frel) \n")

    # --- Calculate U / UMSY (last historical year) aka Frel ---------------------------------------
    UMSY = sapply(msy_out,function(x)(x$UB))
    #Landings_allf = apply(Landings_byfleet_s, c(1,2,3,5),sum)
    Landings_y = apply(array(Landings,c(dms$nMSE,dms$nsim,ns,dms$ny/ns,dms$nMP)),c(1,2,4,5),sum)
    Biomass_y = apply(array(B,c(dms$nMSE,dms$nsim,ns,dms$ny/ns,dms$nMP)),c(1,2,4,5),sum)
    dimnames(Landings_y) = dimnames(Biomass_y) =
      list(OMnames,paste("sim",1:dms$nsim),unique(floor(Time_steps)),names(MSElist[[1]]@MPs))
    U =  Landings_y / Biomass_y
    Uind = MSEtool:::TEG(dim(U))
    Frel = array(U / UMSY[Uind[,c(2,1)]],dim(U))
    dimnames(Frel) = dimnames(U)    # matplot(Frel[1,,,1],type="l"); grid()
    #Frel = aperm(Frel, c(4,1,2,3))

    # --- Calculate SSB / SSBMSY (last historical year) aka Brel -----------------------------------
    SSBMSY = sapply(msy_out,function(x)x$eqSBMSY)
    Sind = MSEtool:::TEG(dim(SSB))
    Brels = array(SSB / SSBMSY[Sind[,c(2,1)]],dim(SSB))
    Brel = apply(array(Brels,c(dms$nMSE,dms$nsim,ns,dms$ny/ns,dms$nMP)),c(1,2,4,5),mean)
    dimnames(Brel) = dimnames(Frel)  # matplot(Brel[2,,,1],type="l") # matplot(Landings[1,,],type="l")


  }
  cat("Returning mahiResults object \n")
  MR = list(Landings = Landings, Landings_byfleet = Landings_byfleet, CDif = CDif, SSB = SSB,
       Brel = Brel, Frel = Frel, VB = VB, VBr = VBr, MuLen = MuLen, PropCAA = PropCAA, dms = dms,
       CurrentYear = MSElist[[1]]@OM@CurrentYear, MPs = MSElist[[1]]@MPs, OMs = OMnames, mahiMSE = packageVersion('mahiMSE'))

  class(MR) = "mahiResults"
  MR

}


getMSEdims = function(MSElist){
  nsim = nSim(MSElist[[1]])
  nt = nYear(MSElist[[1]]@OM)
  np = pYear(MSElist[[1]]@OM)
  ny = nt+np
  na = as.integer(nAge(MSElist[[1]]@OM))
  nf = nFleet(MSElist[[1]]@OM)
  nl = length(MSElist[[1]]@OM@Stock[[1]]@Length@Classes)
  nr = nArea(MSElist[[1]]@OM)
  nMSE = length(MSElist)
  nMP = length(MSElist[[1]]@MPs)
  data.frame(nsim=nsim, nt=nt, np=np, ny=ny, na=na, nf=nf, nl=nl, nr=nr, nMSE = nMSE, nMP = nMP)
}

Combine = function(Reslist, type="OM"){

  newRes = Reslist[[1]]
  slots = names(Reslist[[1]])
  modslots = slots[!slots%in%c("dms","CurrentYear","MPs","OMs","mahiMSE")]

  for(slot in modslots){
    cat(".")
    slotlist = lapply(Reslist,function(x,slot)x[[slot]],slot=slot)
    ndim = length(dim(slotlist[[1]]))
    if(type == "OM") newRes[[slot]] = abind(slotlist, along = 1)    # OM is always first dimension
    if(type == "MP") newRes[[slot]] = abind(slotlist, along = ndim) # MP is always last dimension
  }
  cat(" \n")
  if(type == "OM") newRes$dms$nMSE = sum(sapply(Reslist,function(x)x$dms$nMSE))
  if(type == "MP"){
    newRes$dms$nMP = sum(sapply(Reslist,function(x)x$dms$nMP))
    newRes$MPs = unlist(lapply(Reslist,function(x)x$MPs))
  }
  cat("new dimensions are: \n")
  print(newRes$dms)
  newRes
}

