MR_sub = function(MR, OMs = NULL, MPs = NULL){
  MR_new = MR
  OMnames = MR$OMs
  if(!is.null(OMs)){
    OMind = match(OMs, OMnames)
    MR_new$OMs = MR_new$OMs[OMind]
  }
  MPnames = names(MR$MPs)
  if(!is.null(MPs)){
    MPind = match(MPs,MPnames)
    MR_new$MPs = MR_new$MPs[MPind]
  }
  nslot = length(MR)
  snames = names(MR)
  for(ss in 1:nslot){
    if(class(MR_new[[snames[ss]]])[1]=="array"){
      dims = dim(MR_new[[snames[ss]]])
      nd = length(dims)
      if(!is.null(MPs)){
        if(nd == 4) MR_new[[snames[ss]]]= MR_new[[snames[ss]]][,,,MPind,drop=F]
        if(nd == 5) MR_new[[snames[ss]]]= MR_new[[snames[ss]]][,,,,MPind,drop=F]
        if(nd == 6) MR_new[[snames[ss]]]= MR_new[[snames[ss]]][,,,,,MPind,drop=F]
        MR_new$dms$nMP = length(MPind)
      }
      if(!is.null(OMs)){
        if(nd == 4) MR_new[[snames[ss]]]= MR_new[[snames[ss]]][OMind,,,,drop=F]
        if(nd == 5) MR_new[[snames[ss]]]= MR_new[[snames[ss]]][OMind,,,,,drop=F]
        if(nd == 6) MR_new[[snames[ss]]]= MR_new[[snames[ss]]][OMind,,,,,,drop=F]
        MR_new$dms$nMSE = length(OMind)
      }
    }
  }
  MR_new
}
