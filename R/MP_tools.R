

Fcur_Reclist = function(DataList, Effort){
  Rec=new('list')
  for(ss in 1:length(DataList)){
    Rec[[ss]]=new('list')
    nr = length(DataList[[ss]])
    for(ff in 1:nr){   # fleets as areas in the OM
      Rec[[ss]][[ff]] = new('Rec')
      Rec[[ss]][[ff]]@Effort = Effort

    }
  }
  Rec
}


fakeMP = function(x, DataList, repyr = 2026,...){
  Rec = Fcur_Reclist(DataList, 1)
  yr = max(DataList[[1]][[1]]@Year)
  if(yr >= repyr){
    qq = yr-2022-3
    if(qq==5)stop("Recorded four quarters")
    saveRDS(DataList,paste0("C:/GitHub/DolphinMSE/MP_design/example_mp_data_",qq,".rds"))
  }
  Rec
}
class(fakeMP) = "MMP"


get_F_4_proj = function(fyrs = 4){ # last four seasons

  # Get data to assign catches to Fleet, year and season
  annual.landings = readRDS("data/annual.landings.rds") # total annual landings by large area - fleet (for checking purposes)
  Fcatch = readRDS("data/Catches_s_y_F.rds")            # catch by season, yr, fleet
  Spat_Frac = readRDS("data/Spat_Cat_Frac.rds")         # fraction of exploitation by fleet and region
  dist = readRDS('data/VASTdist.rds')                   # a matrix of historical distributions (n time step x narea)
  index = read.csv('data/VASTindex.csv')$ALL
  Fcat = make_Catch_2(Fcatch, Spat_Frac, dist, index, annual.landings) # yr, f, r  make_catch_2 does spatial (within large area) distribution of catches according to the VAST index and a constant F

  mov =  readRDS('OMs/mov_matrices/mov_50.rds')
  fits = readRDS('fits/fits_mini.rds')  # readRDS('fits/fits_2.rds')
 # MOM_list = lapply(fits, RCM2MOM)
  #MOM=MOM_list[[1]]
  MOM = RCM2MOM(fits[[1]])
  Hist00 = SimulateMOM(MOM) # spatially aggregated reconstruction
  C00 = getC(Hist00)/1000 # cbind(C00[1,],apply(fit@data@Chist,1,sum)/1000, apply(Fcat,1,sum)/1000)
  F00 = getF(Hist00)
  B00 = getSB(Hist00)/1000
  nsim = MOM@nsim
  SMOM = addmov(MOM, mov[1:nsim,,,,])
  anams = sapply(dimnames(mov)[[3]],function(x)strsplit(x,"from_")[[1]][2])

  # 1 Initial historical reconstruction
  Hist0 = SimulateMOM(SMOM) # spatially disaggregated reconstruction
  sels = all_recent_sel(Hist0)

  # 2 The implied spatial F by fleet summed over areas (fleets as areas)
  Nmat = Hist0[[1]][[1]]@AtAge$Number                 # 1 Numbers
  dims = dim(Nmat)
  nsim = dims[1]; na = dims[2]; ny = dims[3]; nr = dims[4] # sim, age(quarter), year(quarter), region
  nf = dim(Fcat)[2] # nf=length(Hist[[1]]);
  Fi = getFinds(Hist0)                                # 3 Apical F last estimated
  Wt = Hist0[[1]][[1]]@AtAge$Weight                   # 5 Weight at age

  gridy = expand.grid(1:nsim, 1:ny, 1:nf, 1:nr)      # stuff to solve
  nrun = nsim*ny*nf*nr                               # numbers of runs

  setup()
  sfExport('solveF_int')
  Fs = sfSapply(1:nrun, solveF, gridy=gridy, Fcat=Fcat, Nmat=Nmat, sels=sels, Fi=Fi, Wt=Wt)
  Fmat = array(Fs,c(nsim,ny,nf,nr))
  Famat = array(NA,c(nsim,ny,nf,na,nr))

  Find = as.matrix(expand.grid(1:nsim,1:ny,1:nf,1:na,1:nr))
  Famat[Find] = Fmat[Find[,c(1,2,3,5)]] * sels[Find[,c(1,4,3)]]

  rind = ny-(fyrs:1)+1
  Fmu = apply(Famat[,rind,,,],2:5,mean) # y, f, a, r    # mean by fleet
  Fapl = apply(Fmu,c(1,2,4),max)        # y, f, r       # apical by fleet

  Finfo = list(F_qfr = Fapl, S_fa = apply(sels,3:2,mean), C_qfr = Fcat[rind,,])
  Finfo
}


getq = function(DataList){
  dd = DataList[[1]][[1]]                 # 'Years' are actually quarters so the labelling is messed up
  rep(1:4,1000)[max(dd@Year)-dd@LHYear+1] # This is the quarter following the data of this quarter - so one time step ahead and the correct index for TAC and Effort distribution
}

makeDMPVars = function(){
  Effort = 1; TAC = NA;
  Smin = NA; Smax = NA; PRM = 0.9; PRM_area = NA; BL = NA;
  F_fr = matrix(c(
    0.000, 0.036, 0.036, 0.047, 0.040, 0.000, 0.049,
    0.000, 2.703, 2.716, 0.000, 0.000, 0.000, 0.000,
    0.000, 0.000, 0.000, 1.329, 1.199, 0.000, 0.000,
    0.000, 1.402, 1.413, 0.000, 0.000, 0.000, 0.000,
    0.000, 0.000, 0.000, 0.167, 0.141, 0.000, 0.000,
    0.463, 0.000, 0.000, 0.000, 0.000, 0.509, 0.434,
    0.239, 0.175, 0.177, 0.224, 0.191, 0.262, 0.226,
    0.342, 0.000, 0.000, 0.000, 0.000, 0.375, 0.322),byrow=T,nrow=8)
  S_fa = matrix(c(
    0.000, 0.002, 0.100, 0.358, 0.596, 0.752, 0.846, 0.902, 0.936, 0.958, 0.972, 0.982, 0.988, 0.993, 0.996, 0.998, 1,
    0.001, 0.293, 0.910, 0.985, 0.996, 0.999, 0.999, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1,
    0.000, 0.271, 0.957, 0.997, 0.999, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1,
    0.000, 0.539, 0.976, 0.996, 0.999, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1,
    0.000, 0.235, 0.952, 0.997, 0.999, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1,
    0.000, 0.002, 0.100, 0.358, 0.596, 0.752, 0.846, 0.902, 0.936, 0.958, 0.972, 0.982, 0.988, 0.993, 0.996, 0.998, 1,
    0.000, 0.002, 0.100, 0.358, 0.596, 0.752, 0.846, 0.902, 0.936, 0.958, 0.972, 0.982, 0.988, 0.993, 0.996, 0.998, 1,
    0.000, 0.002, 0.100, 0.358, 0.596, 0.752, 0.846, 0.902, 0.936, 0.958, 0.972, 0.982, 0.988, 0.993, 0.996, 0.998, 1),
    byrow=T,nrow=8)
  TAC_fr = matrix(c(
    0.00,   134.038,   209.487,   3870.435,   2150.09,     0.00,  4347.824,
    0.00, 33474.924, 50146.395,      0.000,      0.00,     0.00,     0.000,
    0.00,     0.000,     0.000, 129832.700, 107274.31,     0.00,     0.000,
    0.00,  5896.919,  8832.711,      0.000,      0.00,     0.00,     0.000,
    0.00,     0.000,     0.000,  26976.448,  21507.89,     0.00,     0.000,
    32020.84,    0.000,     0.000,      0.000,      0.00, 93566.29, 34453.182,
    18234.03,   670.788,  1079.381,  19723.945,  10178.47, 52209.62, 19261.668,
    24343.37,     0.000,     0.000,      0.000,      0.00, 70757.28, 26257.258),
    byrow=T,nrow=8)
  x=1

  list2env(list(Effort=Effort,  TAC=TAC, Smin=Smin, Smax = Smax,
                PRM=PRM, BL=BL,
                F_fr=F_fr, S_fa=S_fa, TAC_fr=TAC_fr, x=x), envir = .GlobalEnv)
}


get_seasonal_F<-function(x,DataList,lhr){
  qq=rep(1:4,100)[max(DataList[[1]][[1]]@Year)-DataList[[1]][[1]]@LHYear+1]
  ny = ncol(DataList[[1]][[1]]@Misc$FleetPars$Find)
  last_year = ny-(4-qq)
  sapply(DataList[[1]],function(X)X@Misc$FleetPars$Find[x,last_year])
}
