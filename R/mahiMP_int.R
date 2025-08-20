

mahiMPcheck = function(DataList, Effort,  TAC,  Smin, Smax, BL,
                       F_qfr, S_fa , verbose=F){
  nr <- length(DataList[[1]])   # Second level is fleets as areas (same dimensions among stocks)
  #if(!(dim(F_qfr)[2] == dim(S_fa[1]) & nrow(F_fr) == dim(TAC_qfr)[2])) stop("Matrix arguments F_qfr, S_fa and TAC_qfr must have the same number of fleets (rows)")
  nf = dim(F_qfr)[2]

  if(!is.na(Effort[1])){
    if(!is.matrix(Effort)){
      if(!(length(Effort)%in%c(1,nr,nf))) stop(paste0("Effort must be a vector of length 1, ",nr, " or ", nf))
    }else{
      if(!(nrow(Effort)==nf & ncol(Effort)==nr))stop(paste0("Effort matrices should be ",nf, " (fleets) rows by ", nr, " (areas) columns"))
    }
  }

  if(!is.na(TAC[1])){
    if(!is.matrix(TAC)){
      if(!(length(TAC)%in%c(1,nr,nf))) stop(paste0("TAC must be a vector of length 1, ",nr, " or ", nf))
    }else{
      if(!(nrow(TAC)==nf & ncol(TAC)==nr))stop(paste0("TAC matrices should be ",nf, " (fleets) rows by ", nr, " (areas) columns"))
    }
  }

  if(!is.na(Smin[1]) & !length(Smin)%in%c(1,nf)) stop(paste0("Smin must be either a single number or a vector ",nr, " long"))
  if(!is.na(Smax[1]) & !length(Smax)%in%c(1,nf)) stop(paste0("Smax must be either a single number or a vector ",nr, " long"))
  if(!is.na(BL[1]) & length(BL)!=nf) stop(paste0("BL_area must be a vector ",nf, " long"))
  #cat("mahiMPcheck CLEARED \n")
}


map_fleet_dyn = function(Effort=NA, BL =NA, Smin=NA, Smax=NA, F_qfr, S_fa, BGmu, BGcv, qq, Eff_q_mult, hB_ar, B_ar, Len_age, LenCV){  # converts effort instructions by fleet (n=8) to fleet-as-area (n=7)

  F_fr = F_qfr[qq,,]
  nf = dim(F_qfr)[2]; nr = dim(F_qfr)[3]; na = ncol(S_fa)
  F_fra = array(NA,c(nf,nr,na))

  if(!is.matrix(Effort)){         # not effort by fleet and area

    FRA = as.matrix(expand.grid(1:nf,1:nr,1:na)); FR = FRA[,1:2]; FA = FRA[,c(1,3)]; Fo = FRA[,1]; Ro = FRA[,2]
    if(length(Effort)==nf){       # by fleet
      F_fra[FRA] = F_fr[FR] * S_fa[FA] * Effort[Fo] * Eff_q_mult[FR]
    }else if(length(Effort)==nr){ # by area
      F_fra[FRA] = F_fr[FR] * S_fa[FA] * Effort[Ro] * Eff_q_mult[FR]
    }

  }else{                          # effort by fleet and area
    F_fra[FRA] = F_fr[FR] * S_fa[FA] * Effort[FR] * Eff_q_mult[FR]
  }

  F_fr_nu = apply(F_fra,1:2,max)
  E_nu = apply(F_fr_nu,2,sum)/apply(F_fr,2,sum)
  F_ra = apply(F_fra,2:3,sum)
  S_ra = F_ra/apply(F_ra,1,max)
  R_ra = array(1,c(nr,na)) # retain all

  if(!is.na(BL[1])){
    BL_Fmod = getBL_Fmod(BL, qq, BGmu, BGcv, F_fra, hB_ar, B_ar)
    R_fra = array(NA,dim(F_fra))
    RF = FRA[,2:1]
    R_fra[FRA] = F_fra[FRA] * BL_Fmod[RF]
    R_ra = apply(R_fra,2:3,sum) / F_ra
  }

  if(!is.na(Smin[1])|!is.na(Smax[1])){
    R_SL_mod = getSL_Rmod(Len_age, LenCV, Smin, Smax, qq, F_fra)
    R_ra = R_ra * R_SL_mod # product of retention and size limit effect (baglimit affects 'height', Smin Smax affect shape)
  }

  list(S_ra = S_ra, E_nu = E_nu, R_ra = R_ra)
}

# calculate release rate from mean, CV and
calc_BL_RR = function(mu, cv, bl){
  if(!is.na(mu)){
    CR = 1:200
    sd = cv*mu
    var = sd*sd
    size = (mu^2)/(var-mu) # get size parameter
    p = mu/var
    # p2 = size/(size+mu)        #
    # var2 = size*(1-p)/(p*p) # check
    dens = dnbinom(CR,size,mu=mu)
    cond = CR>bl
    nfish = dens*CR
    return(sum(nfish[cond])/sum(nfish))
  }else{
    return(NaN)
  }
}

        #getBL_Fmod(BL, qq, BGmu, BGcv, F_qfr[qq,,], S_fa, hB_ar, B_ar)
getBL_Fmod=function(BL, qq, BGmu, BGcv, F_fra, hB_ar, B_ar){

  nf = nrow(S_fa); na = ncol(S_fa); nr = ncol(B_ar)
  hVB_fra = VB_fra = array(NA,c(nf,nr,na))
  FRA = as.matrix(expand.grid(1:nf,1:nr,1:na)); AR = FRA[,3:2]; FR = FRA[,1:2]; FA = FRA[,c(1,3)]
  hVB_fra[FRA] = hB_ar[AR] * F_fra[FRA]
  hVB_rf = apply(hVB_fra,2:1,sum)
  calib = t(BGmu[qq,,])/hVB_rf
  VB_fra[FRA] = B_ar[AR] * F_fra[FRA]
  VB_rf = apply(VB_fra,2:1,sum)
  mu_rf = VB_rf * calib
  ind_rf = as.matrix(expand.grid(1:nr,1:nf))
  cv_rf = t(BGcv[qq,,])
  RR = sapply(1:(nr*nf),function(x,ind_rf, mu_rf, cv_rf, BL){
    calc_BL_RR(mu_rf[ind_rf[x,,drop=F]], cv_rf[ind_rf[x,,drop=F]], BL[ind_rf[x,2]])},
    ind_rf=ind_rf, mu_rf = mu_rf, cv_rf = cv_rf, BL = BL)
  Fmod = 1-array(RR,c(nr,nf))
  Fmod[is.na(Fmod)] = 1
  Fmod

}

Blank_rec = function(DataList){
  Rec=new('list')
  for(ss in 1:length(DataList)){
    Rec[[ss]]=new('list')
    for(ff in 1:length(DataList[[ss]])){
      Rec[[ss]][[ff]] = new('Rec')
    }
  }
  Rec
}


get_BL_bio = function(x, DataList, qq){
  Bio = DataList[[1]][[1]]@Misc$StockPars$Biomass[x,,129:148,] # last five years
  apply(Bio[,c(1,5,9,13,17)+qq-1,],c(1,3),mean) # mean biomass by age and area for this quarter over the last 5 years.
}


Rec_Effort = function(Rec, x, DataList, Effort, BL, Smin, Smax, F_qfr, S_fa, BGmu, BGcv, qq){

  nr = dim(F_qfr)[3];  nf = dim(F_qfr)[2]
  if(length(Effort)==1)Effort = rep(Effort, nf) # expand to by fleet
  Effort[Effort==0]=1E-5                        # can't be zero, just very small
  Eff_q_mult = F_qfr[qq,,]/F_qfr[4,,]           # the relative increase just due to seasonal variation
  Eff_q_mult[is.na(Eff_q_mult)] = 0             # has to be a real number
  hB_ar = get_BL_bio(x, DataList, qq)           # get mean biomass by age over last 5 years for this quarter
  B_ar = getBio(x, DataList, F_qfr, S_fa)[qq,,]

  cy = length(DataList[[1]][[1]]@Year)
  Len_age = DataList[[1]][[1]]@Misc$StockPars$Len_age[x,,cy]
  na = length(Len_age)
  LenCV = DataList[[1]][[1]]@Misc$StockPars$LenCV[x]

  newFdyn = map_fleet_dyn(Effort, BL, Smin, Smax, F_qfr, S_fa, BGmu, BGcv, qq, Eff_q_mult, hB_ar, B_ar, Len_age, LenCV) # returns both new effort (apical F relative) and new selectivity
  EffortbyArea = newFdyn$E_nu

  for(ss in 1:length(DataList)){
    for(rr in 1:nr){                                    # fleets as areas in the OM
      Rec[[ss]][[rr]]@Effort = EffortbyArea[rr]
      Rec[[ss]][[rr]]@Misc$V_age = newFdyn$S_ra[rr,] # Effort_fleet: selectivity age change due to fleet change
      Rec[[ss]][[rr]]@Misc$R_age = newFdyn$R_ra[rr,] # Effort_fleet: selectivity age change due to fleet change
    }
  }

  Rec
}

map_TAC = function(TAC, BL = NA, Smin = NA, Smax = NA, F_qfr, S_fa, B_qar, qq, hB_ar, B_ar, Len_age, LenCV){

  nq = dim(F_qfr)[1]; nr = dim(F_qfr)[3];  nf = dim(F_qfr)[2]; na = dim(S_fa)[2]
 # F_fr = F_qfr[qq,,]
  cat_pred = VB = S_fra = array(NA,c(nq,nf,nr,na))
  QFRA = as.matrix(expand.grid(1:nq,1:nf,1:nr,1:na)); QFR = QFRA[,1:3]; FR = QFRA[,2:3]
  QAR = QFRA[,c(1,4,3)]; FA= QFRA[,c(2,4)]; AR = QFRA[,4:3]
  cat_pred[QFRA] = B_qar[QAR] * (1-exp(-F_qfr[QFR])) * S_fa[FA]
  cat_qfr = apply(cat_pred, 1:3, sum)

  nQFR = as.matrix(expand.grid(1:nq,1:nf,1:nr)); nFR = nQFR[,2:3]; nF = nQFR[,2]; nR = nQFR[,3]


  if(!is.matrix(TAC)){ # if you have fleet or area specific TAC

    if(length(TAC)==1){

      TAC_qfr = TAC * cat_qfr/sum(cat_qfr)

    }else if(length(TAC)==nf){ # fleet specific

      cat_f = apply(cat_qfr,2,sum)
      cat_frac_byf = array(cat_qfr[QFR] / cat_f[nF],c(nq,nf,nr))
      TAC_qfr = array(TAC[nF] * cat_frac_byf[nQFR],c(nq,nf,nr))

    }else if(length(TAC)==nr){ # area specific

      cat_r = apply(cat_qfr,3,sum)
      cat_frac_byr = array(cat_qfr[QFR] / cat_r[nR],c(nq,nf,nr))
      TAC_qfr = array(TAC[nR] * cat_frac_byr[nQFR],c(nq,nf,nr))

    }

  }else{               # TAC_fr is specified

    cat_fr = apply(cat_qfr,2:3,sum)
    cat_frac_byfr = array(cat_qfr[QFR] / cat_fr[nFR],c(nq,nf,nr))
    cat_frac_byfr[is.na(cat_frac_byfr)]=0
    TAC_qfr = array(TAC[nFR] * cat_frac_byfr[nQFR],c(nq,nf,nr))

  }

  TAC_r = apply(TAC_qfr[qq,,],2,sum) # tac prescription by area.

  ### Selectivity weighted by fleet U
  RFA = as.matrix(expand.grid(1:nr,1:nf,1:na)); RF = RFA[,1:2]; FA = RFA[,2:3]
  B_r = apply(B_qar[qq,,],2,sum) # biomass in each area
  u_rf = t(TAC_qfr[qq,,])/B_r
  s_rfa = array(NA,c(nr,nf,na))
  s_rfa[RFA] = u_rf[RF] * S_fa[FA]
  s_ra_ns = apply(s_rfa,c(1,3),sum)
  S_ra = s_ra_ns/apply(s_ra_ns,1,max) # U weighted (by fleet) selectivity by area
  R_ra = array(1,c(nr,na)) # retain all

  if(!is.na(BL[1])){
    BL_Fmod = getBL_Fmod(BL, qq, BGmu, BGcv, aperm(s_rfa,c(2,1,3)), hB_ar, B_ar)
    R_fra = array(NA,c(nf,nr,na))
    RF = RFA[,1:2]; FRA = RFA[,c(2,1,3)]
    ### you are fixing this!
    R_fra[FRA] = s_rfa[RFA] * BL_Fmod[RF]
    R_ra = apply(R_fra,2:3,sum) / apply(s_rfa,c(1,3),sum)
  }

  if(!is.na(Smin[1])|!is.na(Smax[1])){
    R_SL_mod = getSL_Rmod(Len_age, LenCV, Smin, Smax, qq, aperm(s_rfa,c(2,1,3)))
    R_ra = R_ra * R_SL_mod # product of retention and size limit effect (baglimit affects 'height', Smin Smax affect shape)
  }

  list(TAC_r = TAC_r, S_ra = S_ra, R_ra = R_ra)
}

Rec_TAC = function(Rec, x, DataList, TAC, BL, Smin, Smax, C_qfr, F_qfr, S_fa, BGmu, BGcv, qq){

  nr = dim(F_qfr)[3];  nf = dim(F_qfr)[2]; na = dim(S_fa)[2]
  TAC[TAC==0]=1E-5
  TAC_q_frac = (apply(C_qfr,1,sum)/sum(C_qfr))[qq]
  TAC_rq_frac = apply(C_qfr,c(3,1),sum)/apply(C_qfr,3,sum)
  TAC_r_frac = apply(C_qfr[qq,,],2,sum) / sum(C_qfr[qq,,])
  B_qar = getBio(x, DataList, F_qfr, S_fa)
  hB_ar = get_BL_bio(x, DataList, qq)           # get mean biomass by age over last 5 years for this quarter
  B_ar = getBio(x, DataList, F_qfr, S_fa)[qq,,]

  cy = length(DataList[[1]][[1]]@Year)
  Len_age = DataList[[1]][[1]]@Misc$StockPars$Len_age[x,,cy]
  na = length(Len_age)
  LenCV = DataList[[1]][[1]]@Misc$StockPars$LenCV[x]

  map = map_TAC(TAC, BL, Smin, Smax, F_qfr, S_fa, B_qar, qq, hB_ar, B_ar, Len_age, LenCV)

  for(ss in 1:length(DataList)){
    for(rr in 1:length(DataList[[ss]])){
      Rec[[ss]][[rr]]@TAC= map$TAC_r[rr]
      Rec[[ss]][[rr]]@Misc$V_age = map$S_ra[rr,] # TAC_fleet: selectivity age change due to fleet change
      Rec[[ss]][[rr]]@Misc$R_age = map$R_ra[rr,]  # retention calculated based on reduction of harvest rate
    }
  }
  Rec

}



Rec_Fdisc=function(Rec, Fdisc){
  nr = length(Rec[[1]])
  if(length(Fdisc)==1)Fdisc = rep(Fdisc,nr)
  for(ss in 1:length(Rec)){
    for(rr in 1:length(Rec[[1]])){
      Rec[[ss]][[rr]]@Fdisc = Fdisc[rr]
    }
  }
  Rec
}

# To distribute TACs sensibly you need the vulnerable biomass from the most recent relevant quarter
# getBio gets the biomass by age and area for the most recent quarter (actually the next quarter in the model)

getBio = function(x, DataList, F_qfr, S_fa, VB=F){
  nr = dim(F_qfr)[3];  nq = dim(F_qfr)[1]; na = dim(S_fa)[2]
  Bio = array(NA,c(nq,na,nr))
  ny = dim(DataList[[1]][[1]]@Misc$StockPars$Biomass)[3]
  cy = length(DataList[[1]][[1]]@Year)
  last4ind = rep(1:4,1000)[1:cy]
  last4B = last4ind[length(last4ind)-3:0]
  yind = max((1:cy)[last4ind==4])-3:0
  #yind = (1:cy)[cy-3:0]
  for(qq in 1:nq){
    yy = yind[qq]
    #qind = last4B[qq]
    if(!VB){ # if not vulnerable biomass
      if(yy <= ny) Bio[qq,,] = DataList[[1]][[1]]@Misc$StockPars$Biomass[x,,yy,]        # annoying - need to find it in historical or projected arrays
      if(yy >ny) Bio[qq,,] = DataList[[1]][[1]]@Misc$StockPars$Biomass_P[x,,yy-ny,]  # annoying - need to find it in historical or projected arrays
    }else{
      if(yy <= ny) Bio[qq,,] = DataList[[1]][[1]]@Misc$StockPars$VBiomass[x,,yy,]        # annoying - need to find it in historical or projected arrays
      if(yy >ny) Bio[qq,,] = DataList[[1]][[1]]@Misc$StockPars$VBiomass_P[x,,yy-ny,]  # annoying - need to find it in historical or projected arrays
    }
  }
  Bio
}

# This calculates the implied retention due to a min and or max size limit by fleet
getSL_Rmod=function(Len_age, LenCV, Smin, Smax, qq, F_fra){
  dimF = dim(F_fra); nf = dimF[1]; nr = dimF[2]; na = dimF[3]
  nuF = array(NA,dim(F_fra))
  for(ff in 1:nf){
    Fmod = recalc_sel_sizelim(Len_age, LenCV, Smin[ff], Smax[ff])
    for(rr in 1:nr){
      nuF[ff,rr,] = F_fra[ff,rr,] * Fmod
    }
  }
  F_ra = apply(F_fra,2:3,sum)
  nuF_ra = apply(nuF,2:3,sum)
  nuF_ra/F_ra  # predicts the retention change due to Smin and Smax
}


recalc_sel_sizelim=function(Len_age,LenCV, Smin=500, Smax=1100, nl = 100, maxlen_mult = 1.4){
  lens = seq(0,maxlen_mult*max(Len_age),length.out=nl)
  na = length(Len_age)
  dens = array(0,c(na,nl))
  SLeffect = rep(1,nl)
  SLeffect[lens<Smin] = 0
  SLeffect[lens>Smax] = 0
  for(aa in 1:na)    dens[aa,] = dlnorm(lens,log(Len_age[aa]),LenCV)
  tdens = dens/apply(dens,1,sum)
  for(aa in 1:na) tdens[aa,] = tdens[aa,]*SLeffect
  apply(tdens,1,sum)
}
