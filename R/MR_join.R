MR_join = function(MR1, MR2, type = "MP"){
  MRnew = MR1
  nslot = length(MR1)
  snames = names(MR1)
  for(ss in 1:nslot){
    dims = dim(MRnew[[snames[ss]]])
    if(class(MRnew[[snames[ss]]])=="array"){
      if(type == "MP"){
        along = length(dims)
        MRnew$dms$nMP = MR1$dms$nMP + MR2$dms$nMP
      }else if(type =="OM"){
        along = 1 # OM dim is always 1
        MRnew$dms$nMSE = MR1$dms$nMSE + MR2$dms$nMSE
      }else if(type =="Sim"){
        along = 2
        MRnew$dms$nsim = MR1$dms$nsim + MR2$dms$nsim
      }
      MRnew[[snames[ss]]] = abind::abind(MR1[[snames[ss]]],MR2[[snames[ss]]], along=along)
    }

  }
  if(type=="MP")MRnew$MPs = c(MR1$MPs,MR2$MPs)
  MRnew
}
