

mahiSlick = function(MR,
  Measures = c("Brel","Frel"),
  Mnames = c("SSB/SSBMSY", "F/FMSY"),
  Yrs = list(c(2026,2040),c(2026, 2030),c(2031, 2035),c(2036, 2040)),
  ynames= c("all","S","M","L")){

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
      stat = apply(val[,,tind,,drop=F],c(1,2,4),mean)
      ylab = ynames[yy]
      MM[[Measures[mm]]][[ylab]] = stat
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
    stat = apply(Brel[,,tind,,drop=F]>1 & Frel[,,tind,,drop=F]<1,c(1,2,4),mean)
    ylab = ynames[yy]
    MM[["PGK"]][[ylab]] = stat
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
      mlab = paste0(Fleets[ft]," ",labels[fm])
      MM[[mlab]] = list()

      for(yy in 1:ny){

        tind = ys>=Yrs[[yy]][1] & ys<=Yrs[[yy]][2]
        stat = apply(val[,,tind,,drop=F],c(1,2,4),mean,na.rm=T)
        ylab = ynames[yy] #ylab = paste0(Yrs[[yy]][1],"-", Yrs[[yy]][2])
        MM[[mlab]][[ylab]] = stat
      }

    }
  }

  MM
}









