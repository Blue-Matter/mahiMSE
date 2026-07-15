

mahiMetrics = function(MR,
  Measures = c("Brel","Frel"),
  Mnames = c("SSB/SSBMSY", "F/FMSY"),
  RefLevs = list(c(0.5,1),1),
  RefDir = c("L","H"),
  Yrs = list(c(2026,2040),c(2026, 2030),c(2031, 2035),c(2036, 2040)),
  ynames= c("all","S","M","L"),
  Quantiles = c(0.05,0.25, 0.5, 0.75, 0.95)){

  # MR = mahiResults(myMSE); Measures = c("Brel","Frel");   Mnames = c("SSB/SSBMSY", "F/FMSY"); RefLevs = list(c(0.5,1),1); RefDir = c("L","H"); Yrs = list(c(2026,2040),c(2026, 2030),c(2031, 2035),c(2036, 2040));  ynames= c("all","S","M","L"); Quantiles = c(0.05,0.25, 0.5, 0.75, 0.95)
  ns = 4
  nm = length(Measures)
  ny = length(Yrs)
  dms = MR$dms
  MPs = dimnames(MR$Landings)[[4]]
  MM = list()

  # --- MSY quants -----------------------------------------------------------------
  for(mm in 1:nm){
    val = MR[[Measures[mm]]]
    dims = dim(val)
    mpdim = length(dims) # mp is always last dimension
    ydim = match(dms$ny/ns, dims) # find the year dimension (usually 3)
    ys = as.numeric(dimnames(val)[[ydim]]) # get the years
    MM[[Measures[mm]]] = list()
    for(yy in 1:ny){

      tind = ys>=Yrs[[yy]][1] & ys<=Yrs[[yy]][2]
      qval = apply(val[,,tind,,drop=F],mpdim,quantile, Quantiles)
      mu = apply(val[,,tind,,drop=F],mpdim,mean)
      sd = apply(val[,,tind,,drop=F],mpdim,sd)
      cv = sd/mu

      if(!is.na(RefLevs[[mm]][1])){
        RL = NULL
        for(rl in 1:length(RefLevs[[mm]])){
          RLt = apply(val[,,tind,,drop=F]<RefLevs[[mm]][rl],mpdim,mean)
          if(RefDir[mm] =="H") RLt = 1-RLt
          RL=rbind(RL,RLt)
        }
        opp = c("<",">")[as.integer(RefDir[mm]=="H")+1]
        row.names(RL) = paste0("P(",Measures[mm],opp,RefLevs[[mm]],")")
      }

      mets = rbind(mu,sd,cv,qval,RL)
      colnames(mets) = MPs
      mets = as.data.frame(t(mets))
      ylab = ynames[yy] #paste0(Yrs[[yy]][1],"-", Yrs[[yy]][2])
      MM[[Measures[mm]]][[ylab]] = mets
    }
  }


  # PGK
  Brel = MR$Brel
  Frel = MR$Frel
  dims = dim(Brel)
  mpdim = length(dims) # mp is always last dimension
  ydim = match(dms$ny/ns, dims) # find the year dimension (usually 3)
  ys = as.numeric(dimnames(Brel)[[ydim]]) # get the years
  MM[["PGK"]] = list()
  for(yy in 1:ny){

    tind = ys>=Yrs[[yy]][1] & ys<=Yrs[[yy]][2]
    mets = matrix(apply(Brel[,,tind,,drop=F]>1 & Frel[,,tind,,drop=F]<1,mpdim,mean),nrow=1)
    mets = as.data.frame(mets)
    colnames(mets) = MPs
    rownames(mets) = "mu"
    ylab = ynames[yy] #paste0(Yrs[[yy]][1],"-", Yrs[[yy]][2])
    MM[["PGK"]][[ylab]] = t(mets)
  }


  # --- Fleet specific metrics -------------------------------------------------

  slots =  c("Landings_byfleet", "VB",         "VBr",            "CDif" ,     "MuLen")
  labels = c("Landings",         "Catch Rate", "Ret Catch Rate", "Catch Dif", "Mean length")

  nfm = length(slots)
  for(fm in 1:nfm){
    for(ft in 1:5){

      val = array(MR[[slots[fm]]][,,,ft,],c(dms$nMSE,dms$nsim,dms$ny/ns,dms$nMP))
      dims = dim(val)
      mpdim = length(dims) # mp is always last dimension
      ydim = match(dms$ny/ns, dims) # find the year dimension (usually 3)
      #ys = as.numeric(dimnames(val)[[ydim]]) # get the years
      mlab = paste0(Fleets[ft]," ",labels[fm])
      MM[[mlab]] = list()

      for(yy in 1:ny){

        tind = ys>=Yrs[[yy]][1] & ys<=Yrs[[yy]][2]
        qval = apply(val[,,tind,,drop=F],mpdim,quantile, Quantiles,na.rm=T)
        mu = apply(val[,,tind,,drop=F],mpdim,mean)
        sd = apply(val[,,tind,,drop=F],mpdim,sd)
        cv = sd/mu
        mets = rbind(mu,sd,cv,qval)
        colnames(mets) = MPs
        mets = as.data.frame(t(mets))

        ylab = ynames[yy] #ylab = paste0(Yrs[[yy]][1],"-", Yrs[[yy]][2])
        MM[[mlab]][[ylab]] = mets
      }

    }
  }

  MM
}
