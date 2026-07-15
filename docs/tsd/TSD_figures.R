

# === TSD Figures ==============================================================

library(MSEtool)
library(mahiMSE)

setwd("C:/GitHub/DolphinMSE")
source("code/OM_new_spatial/dol.input.data.R")
fit = readRDS("fits/microfit_March_7_2026.rds")
fleetnames = c("USCom","RecN","RecS","HireN","HireS","Intl","Disc","UnRep")

C_yfr = C_yfr2 = readRDS('data/MK_catch_sep_2025/C_yfr.rds')
C_yfr[,7:8,] = C_yfr[,7:8,] * 0.5                      # level 1 is half catch
C_yfr_list = list(C_yfr, C_yfr2); names(C_yfr_list) = c("Cat_0.5", "Cat_1")



# ---- F1 Natural Mortality rate --------------------------------------
jpeg('C:/GitHub/mahiMSE/docs/tsd/img/Data/M_stoch.jpg', res=300, height=3.5, width=5, units="in")

  par(mai=c(0.9,0.9,0.05,0.05))
  hist(MSEtool::trlnorm(100000,0.25,0.1),xlim=c(0,0.7), border="white",col="lightblue", main="",xlab="Natural Mortality (M) (annual)",ylab="Freq."); grid()
  hist(MSEtool::trlnorm(100000,0.5,0.1),add=T, col="grey",border='white')

dev.off()



# ---- F2 Steepness ---------------------------------------------------

jpeg('C:/GitHub/mahiMSE/docs/tsd/img/Data/h_stoch.jpg', res=300, height=3.5, width=5, units="in")

  par(mai=c(0.9,0.9,0.05,0.05))

  val1 = 0.7; val2 = 0.95; cv = 0.025
  h1 = rbeta(1E5,MSEtool::alphaconv(val1,cv),MSEtool::betaconv(val1,cv))
  h2 = rbeta(1E5,MSEtool::alphaconv(val2,cv),MSEtool::betaconv(val2,cv))

  hist(h1,xlim=c(0.6,1), border="white",col="lightblue", main="",xlab="Bev.-Holt Steepness (h)",ylab="Freq."); grid()
  hist(h2,add=T, col="grey",border='white')

dev.off()


# ---- F4 Recruitment -------------------------------------------------


jpeg('C:/GitHub/mahiMSE/docs/tsd/img/Data/Rec_factor.jpg', res=300, height=5.5, width=7.5, units="in")

  load(file = "C:/GitHub/DolphinMSE/OMs/Reference_set/OM_1.RData")
  esty = 17:148
  recs = OM_1@Stock[[1]]@SRR@RecDevHist[1,esty,drop=F]
  tlab = hTime_steps[esty]
  nsim = dim(recs)[1]
  ny = dim(recs)[2]
  nrs = 16
  yind = ny-nrs:1+1
  yh = 1:(ny-nrs)
  mur = mean(recs[,yind])
  muh = mean(recs[,yh])
  rat = mur/muh

  par(mai=c(0.9,0.9,0.05,0.05))
  plot(tlab, recs[1,], col='blue',pch=19, xlab="Year",ylab="Rec. Dev.");grid()
  abline(v=tlab[ny-nrs],col="grey",lty=2,lwd=2)
  points(tlab[yind], recs[1,yind],col="red",pch=19)
  abline(h=c(mur,muh),col=c("red","blue"),lty=2,lwd=2)

dev.off()


# ---- Weight at length --------------------------------------------------------

# --- Recreational size data ---

  fInd = 2:5
ages = 0:OM@maxage
MRIP = read.csv("data/Copy of SA_DOLPHIN lengths.csv")
Kcols = names(MRIP) %in% c("YEAR","month","season","CATCH","RELEASE","HARVEST","tot_len","wgt_ab1")
La = OM@Linf[1] * (1-exp(-OM@K[1]*(ages-OM@t0[1])))
Wa = OM@a*La^OM@b

Modes = c(rep("Private/Rental",2),rep("Charter",2))
Sstates = c("FLORIDA","GEORGIA","SOUTH CAROLINA","NORTH CAROLINA")
Nstates =  c("VIRGINIA", "DELAWARE","MARYLAND","NEW JERSEY","NEW YORK","CONNECTICUT","RHODE ISLAND","MASSACHUSETTS")
States = list(Nstates, Sstates, Nstates, Sstates)

doCAL = function(sdat,OM,size_bins,syr = 1986){
  ns = length(size_bins)
  nt = OM@nyears
  ts = rep(NA,nrow(sdat))
  fCAL = array(0,c(nt,ns-1))
  for(i in 1:nrow(sdat)){
    if(sdat$month[i] == 12){
      yind = sdat$YEAR[i]-syr+1
    }else{
      yind = sdat$YEAR[i]-syr
    }
    sind = match(sdat$season[i],c("Winter","Spring","Summer","Autumn"))
    ts[i] = yind*4+sind-1
    lind = (1:(ns-1))[sdat$tot_len[i] <= size_bins[2:ns] & sdat$tot_len[i] > size_bins[1:(ns-1)]]
    fCAL[ts[i],lind] =  fCAL[ts[i],lind]+1
  }
  fCAL
}

jpeg('man/img/Data/Rec_WL.jpg', res=300, height=8.5, width=5.5, units="in")

  par(mfrow=c(4,2),mai=c(0.4,0.4,0.01,0.01),omi=c(0.3,0.3,0.01,0.01))
  for(i in 1:length(Modes)){
    sdat = MRIP[MRIP$MODE_FX == Modes[i] & MRIP$ST %in% States[[i]] & MRIP$YEAR >= syr & MRIP$YEAR <= eyr,Kcols]
    sdat = sdat[sdat$tot_len > 0 & sdat$HARVEST ==1 & sdat$tot_len < 1500,]
    plot(sdat$tot_len,sdat$wgt_ab1,col='grey',xlab = "Observed length (mm)",ylab ="Observed weight (kg)");grid()
    legend('topleft',legend=Fnames[i+1],bty="n")
    predW = OM@a*sdat$tot_len^OM@b
    keeprat = sdat$wgt_ab1/predW
    keep = keeprat > (2/3) & keeprat < (3/2)
    sdat = sdat[keep,]
    points(sdat$tot_len,sdat$wgt_ab1,col="#0000ff90",pch=19); lines(La,Wa,col="red")
    CAL[,,fInd[i]] = doCAL(sdat,OM,size_bins,syr=syr)
    hist(sdat$tot_len,main="",breaks = size_bins[1:15])

  }
  mtext("Length (mm)",1,line=0.3,outer=T)
  mtext("Weight (lbs)",2,line=0.3,outer=T)

dev.off()

jpeg('C:/GitHub/mahiMSE/docs/tsd/img/Data/Growth_stoch.jpg', res=300, height=4.5, width=8, units="in")

  nsim = 32*40
  nsyr = 4
  mus = log(c(1289, 1.27 / nsyr))
  corl = -0.6
  cv_Linf = 25.95 / exp(mus[1])
  cv_K = 0.08/1.27
  sds = c(cv_Linf^2,cv_K^2)^0.5
  covar = corl * prod(sds)
  vcv = matrix(c(cv_Linf^2, covar, covar, cv_K^2), nrow = 2, byrow=T)
  samps = rmvnorm(nsim,mus,vcv)
  Linf = exp(samps[,1])
  K = exp(samps[,2])
  plot(Linf, K,xlab ="Asymptotic Length (mm) ('L-infinity')", ylab = "Maximum Growth Rate (K)")
  abline(v = exp(mus[1]), col="red", h = exp(mus[2]), lty=2,lwd=2)
  ages = 1:16 - 0.5
  nages = length(ages)
  agemat = array(rep(ages,each=nsim),c(nsim,nages))
  La = Linf * (1 - exp(-(K)*agemat))
  mLa = exp(mus[1]) * (1 - exp(-(exp(mus[2])*ages)))

  par(mfrow=c(1,2),mai = c(0.9,0.9,0.2,0.05))
  plot(Linf, K,xlab ="Asymptotic Length (mm) ('L-infinity')", ylab = "Maximum Growth Rate (K)", col="dark grey"); grid()
  points(exp(mus[1]),exp(mus[2]),pch=19,cex=1.2, col="red")
  matplot(ages, t(La[1:500,]), xlab="Age Class (quarter)", ylab = "Length(mm)",col="#99999910",type="l",lwd=2,lty = 1); grid(); abline(v=1:2,col="grey",lty=3)
  abline(h = c(458,560),col="blue")
  lines(ages,mLa,col="red",lwd=2)

dev.off()

# ---- Maturity ----------------------------------------------------------------


jpeg('C:/GitHub/mahiMSE/docs/tsd/img/Data/Mat.jpg', res=300, height=3.8, width=7.5, units="in")

  par(mai=c(0.6,0.3,0.05,0.15),omi=c(0.02,0.45,0.01,0.01))
  ages = 0:6
  mats = fit@OM@cpars$Mat_age[1,ages+1,1]
  plot(0:6,mats,pch=19,col="blue");grid()
  lines(0:6,mats,col="blue")
  mtext("Age class",1,line=2.0)

  #lens = fit@OM@Linf[1] * (1-exp(-fit@OM@K[1]* (ages - fit@OM@t0[1])))
  #plot(lens,mats,pch=19,col="blue");grid()
  #lines(lens,mats,col="blue")
  #mtext("Length at start of age class (mm)",1,line=2.0)

  mtext("Spawning fraction",2,line=0.8,outer=T)

dev.off()

# ---- Annual Catches ----------------------------------------------------------

# worm
colnames(annual.landings) = fleetnames
cols = c("black","red","red","green","green","blue", "orange","purple")
lty = c(1,1,3,1,3,1,1,3)
omlab = c("Catch level 1: 50% Disc and UnRep", "Catch level 2: 100% Disc and UnRep")

jpeg('man/img/Data/Annual_landings.jpg', res=300, height=6, width=12, units="in")
  par(mfrow=c(1,2),mai=c(0.35,0.85,0.45,0.05),omi=c(0.4,0.05,0.05,0.05))
  for(om in 1:2){
    ql = apply(C_yfr_list[[om]],1:2,sum)/1000
    al = apply(array(ql,c(4,148/4,8)),2:3,sum)
    matplot(year_start:year_end,al,type="l",
            col=cols,lty=lty,lwd=2, ylab="Landings (mt)",
            xlab =""); grid()
    mtext(omlab[om],line=0.3,font=2)

  }
  mtext("Year",1,line=0.3,outer=T)
  legend('topright',fleetnames,col=cols,lty=lty,text.col=cols,bty="n")
dev.off()



# bar

cols = c("black","red","red","green","green","blue", "orange","purple")
lty = c(1,1,3,1,3,1,1,3)


jpeg('tsd/img/Data/Landings_bar_2022.jpg', res=300, height=8, width=15, units="in")

  cols = c("black","red","red","green","green","blue", "orange","purple")

  tind = 145:148
  CT = C_yfr_list[[2]]
  par(mfrow=c(1,2),mai=c(0.85,0.85,0.45,0.85),omi=c(0.05,0.05,0.05,0.05))

  CTf = apply(CT[tind,,],2,sum)/1E6
  barplot(CTf, ylab="Catch in 2022 (kt)",xlab="Fleet",col=cols,ylim=c(0,max(CTf)*1.1));grid()
  CTfp = CTf / sum(CTf)*100
  rat = max(CTfp)/max(CTf)
  percs = pretty(seq(0,max(CTfp),length.out=8))
  abline(h = percs/rat,col="light blue")
  barplot(CTf, add=T, col=cols)
  axis(4,at=percs/rat,percs)
  mtext("Percentage of total catches",4,line=1.8)
  text(((1:8)-0.45)*1.2,CTfp/rat+0.1,paste0(round(CTfp,1),"%"),cex=0.9)


  CTr = apply(CT[tind,,],3,sum)/1E6
  barplot(CTr, ylab="Catch in 2022 (kt)",xlab="Area",ylim=c(0,max(CTr)*1.1));grid()
  CTrp = CTr / sum(CTr)*100
  rat = max(CTrp)/max(CTr)
  percs = pretty(seq(0,max(CTrp),length.out=8))
  abline(h = percs/rat,col="light blue")
  barplot(CTr, add=T)
  axis(4,at=percs/rat,percs)
  mtext("Percentage of total catches",4,line=1.8)
  text(((1:5)-0.45)*1.2,CTrp/rat+0.15,paste0(round(CTrp,1),"%"),cex=0.9)


dev.off()


# Equilibrium Catches
jpeg('C:/GitHub/mahiMSE/docs/tsd/img/Data/EqCat.jpg', res=300, height=8, width=9, units="in")

  Ct = apply(C_yfr_list[[2]],1:2,sum)/1000
  incl = 1:14
  padding = 1
  tpad = hTime_steps[1] - (1:padding)/4
  cpad = array(NA,c(padding,dim(Ct)[2]))
  dat = rbind(cpad,Ct[incl,])
  lab = c(tpad,hTime_steps[incl])

  par(mai=c(0.85,0.85,0.45,0.85),omi=c(0.05,0.05,0.05,0.05))

  cols = c("black","red","red","green","green","blue", "orange","purple")
  lty = c(1,1,3,1,3,1,1,3)
  matplot(lab,dat,type="l",lty=lty,col=cols,xlab = "Initial Years", ylab="Catch (t)",ylim=c(0,2000));grid()

  yno = 2
  matplot(hTime_steps[1:(yno*4)],Ct[1:(yno*4),], type="p",pch=19,col=cols,add=T)

  muC = apply(Ct[1:(yno*4),],2,mean)
  xmarkers = hTime_steps[c(1,yno*4)]+c(-0.125,0.125)
  abline(v=xmarkers,lwd=2,col="#99999990",lty=2)

  for(ff in 1:length(muC)) lines(c(hTime_steps[1],hTime_steps[yno*4]),rep(muC[ff],2),col=cols[ff],lty=lty[ff],lwd=2)
  lnam = paste0(fleetnames," (",round(muC,1),"t)")
  legend('topright',lnam,col=cols,lty=lty,text.col=cols,bty="n")

dev.off()


# ---- Catch at length  --------------------------------------------------------

CAL = fit@data@CAL
bins = fit@data@length_bin
fleets = 1:5; nf=length(fleets)
byfleet = apply(CAL[,,fleets],2:3,sum)
binmax = 1350
keep = bins<=binmax

jpeg('man/img/Data/CAL.jpg', res=300, height=5, width=7, units="in")

  par(mfrow =c(floor(nf/2),ceiling(nf/2)),mai=c(0.3,0.3,0.25,0.05),omi=c(0.25,0.25,0.02,0.02))
  for(i in 1:nf){
    barplot(byfleet[keep,i],names.arg =bins[keep],col="white",border=NA); grid()
    barplot(byfleet[keep,i],names.arg =bins[keep],col="slategrey",border=NA,add=T)
    mtext(fleetnames[i],line=0.25,cex=0.9)
  }
  mtext("Length (mm)",1,outer=T,line=0.25)
  mtext("Frequency",2,outer=T,line=0.25)

dev.off()



# ---- Total index -------------------------------------------------------------

#xs = rep(1986:2022,each=4)+c(0,0.25,0.5, 0.75)
#all = as.vector(d$ALL)
#all = apply(dist,1,sum)

oVAST <- read.csv('data/New_VAST/VAST_index_10.21.2025.csv')
tot = oVAST[oVAST$STRATA == 'Total',]
ts =tot$Year
all = tot$value

jpeg('tsd/img/Data/Total_index.jpg', res=300, height=4.5, width=7, units="in")
  par(mai=c(0.8,0.8,0.05,0.05))
  plot(xs,all,type='l',col="grey",ylim=c(0,max(all,na.rm=T)*1.025),xlab="Year",ylab="VAST total index");grid()
  cols = c("red","green","blue","purple")
  points(xs,all,pch=19,col=cols)
  legend('topright',c("Wint","Spr","Sum","Aut"),text.col=cols,bty="n")
dev.off()


# ----- Spatial indices --------------------------------------------------------

unique(oVAST$STRATA)
Areas = c("CAR+FLK", "NCA",  "NCFL", "NED", "NNC+VBM")

jpeg('tsd/img/Data/Area_index.jpg', res=300, height=7.5, width=7, units="in")
  par(mfrow=c(3,2),mai=c(0.4,0.3,0.25,0.05),omi=c(0.35,0.35,0.01,0.01))

  for(rr in 1:length(Areas)){
    rv = oVAST[oVAST$STRATA == Areas[rr],]
    ts =rv$Year;    all = rv$value / 1E6
    plot(xs,all,type='l',col="grey",ylim=c(0,max(all,na.rm=T)*1.025),xlab="",ylab="");grid()
    mtext(Areas[rr],line=0.2)
    cols = c("red","green","blue","purple")
    points(xs,all,pch=19,col=cols)
  }
  plot(1,1,col='white',axes=F,xlab="",ylab="",main="")
  legend('center',c("Wint","Spr","Sum","Aut"),text.col=cols)
  mtext("VAST Spatial Index (m)",2,line=0.9,outer=T)
  mtext("Year",1,line=0.4,outer=T)
dev.off()


jpeg('tsd/img/Data/Area_index_not_scaled.jpg', res=300, height=7.5, width=7, units="in")
  par(mfrow=c(3,2),mai=c(0.4,0.3,0.25,0.05),omi=c(0.35,0.35,0.01,0.01))
  ymax = 30
  for(rr in 1:length(Areas)){
    rv = oVAST[oVAST$STRATA == Areas[rr],]
    ts =rv$Year;    all = rv$value / 1E6
    plot(xs,all,type='l',col="grey",ylim=c(0,ymax),xlab="",ylab="");grid()
    mtext(Areas[rr],line=0.2)
    cols = c("red","green","blue","purple")
    points(xs,all,pch=19,col=cols)
  }
  plot(1,1,col='white',axes=F,xlab="",ylab="",main="")
  legend('center',c("Wint","Spr","Sum","Aut"),text.col=cols)
  mtext("VAST Spatial Index (m)",2,line=0.9,outer=T)
  mtext("Year",1,line=0.4,outer=T)
dev.off()



# ---- Data Simulation --------------------------------------------

library(mahiMSE)
library(mahiRefSet)

getCVdist = function(x)x@Landings@CV
statlist = list()
for(oo in 1:32){
  OM = get(paste0('OM_',oo))
  statlist[[oo]] = sapply(OM@Obs$Dolphinfish,getCVdist)
}
Data = Example_Data_0

jpeg('C:/GitHub/mahiMSE/docs/tsd/img/Data/PLL_sim.jpg', res=300, height=3, width=12, units="in")

  Sim_Ind = gen_ind_2(Data, Index_names="PLL", Index_cv=NA, Index_ac=NA,         # Generate abundance indices [2]
                      Index_lag=1, Index_seed=1, Index_beta=1, plot=T, return_stats = T)

  MSEtool:::Calc_Stats(Sim_Ind$resid[[1]])
dev.off()

# do the in-MP stat calcs for OM_1
myMP = function(Data){
  # Data = Example_data
  #  Data = readRDS("C:/temp/Data1.rds")
  season = get_season(Data)
  lhy = match(Data@YearLH, Data@Years) # last historical time step
  FDeadArea = Data@Misc$DataOM@FDeadArea$Dolphinfish[1,,lhy,,]
  na = dim(FDeadArea)[1]; nf = dim(FDeadArea)[2]; nr = dim(FDeadArea)[3]
  N = Data@Misc$DataOM@Number$Dolphinfish[1,,,]
  nt = dim(N)[2]  # Total time steps (historical and projected)
  lhobs = length(Data@Years) # last observation
  subIdat = Indices[Indices$Name %in% "PLL",]
  INames = "PLL"
  rowind = match("PLL", subIdat$Name)
  INo = subIdat$Number[rowind]
  ISeason = subIdat$Season[rowind]
  IArea = subIdat$Area[rowind]
  IType = subIdat$Type[rowind]
  nI = length(INames)
  ii=1
  Ind_yrs = Indices$Year[Indices$Name == "PLL"]

  ft = 1
  sel = Data@Misc$DataOM@OM@Fleet[[1]][[ft]]@Selectivity@MeanAtAge[1,,,] # age, ts, area
  wta = Data@Misc$DataOM@OM@Stock$Dolphinfish@Weight@MeanAtAge[1,,]      # age, ts
  I_area_vec = rep(0,nr); I_area_vec[IArea[[ii]]] = 1
  VBA = N * array(sel, c(na, nt, nr)) * array(rep(I_area_vec, each = na * nt), c(na, nt, nr)) * array(wta, c(na, nt, nr))
  var = apply(VBA,2, sum) / 1E6

  Ind_yrs_s = Ind_yrs + (ISeason[ii]/4)-1/4
  keep = match(Ind_yrs_s,Time_steps)
  keep = keep[keep<=lhy] # only historical years for stat calcs
  hvar = var[keep]

  Ind0 = Indices$Index[Indices$Name == "PLL"]
  Ind = Ind0[Ind_yrs <= Data@YearLH]
  plot(Ind/mean(Ind), hvar/mean(hvar))
  res = MSEtool:::CalcIndexResiduals(Ind/mean(Ind),matrix(hvar/mean(hvar),nrow=1))$LogResiduals[1,]
  df = MSEtool:::Calc_Stats(res)

  saveRDS(Data,"C:/temp/Data1.rds")
  write.table(df,file="C:/Users/tcar_/Dropbox/temp/Dolphin/mahiPLLstats.csv",sep=",",append=T,quote = FALSE, col.names=FALSE,row.names=FALSE)
  if(max(Data@Years)>2022)stop()

}
class(myMP) = "mp"

hist = Simulate(smallOM)#, DoMSYRefs = FALSE, DoRefLandings = FALSE, DoRefRemovals = FALSE, DoGenerateData = FALSE)
test = Project(hist,"myMP")


# mahiMP

jpeg('C:/GitHub/mahiMSE/docs/tsd/img/MP/TL_reb_dep.jpg', res=300, height=8, width=8, units="in")

  par(mfrow=c(2,2),mai=c(0.6,0.7,0.4,0.01))
  mahiMP(Example_Data_dep,plot=T)
  mahiMP(Example_Data_reb,plot=T)

dev.off()

jpeg('C:/GitHub/mahiMSE/docs/tsd/img/MP/TL_reb_dep_model.jpg', res=300, height=8, width=8, units="in")

  par(mfrow=c(2,2),mai=c(0.6,0.7,0.4,0.01))
  mahiMP(Example_Data,plot=T, TL_empirical = T)
  mahiMP(Example_Data,plot=T, TL_empirical = F)

dev.off()


# MPs

source("C:/GitHub/DolphinMSE/code/Analyses/CER_for_PR/CER_source.R")

#overwrite SQ with constant catch HCR for plotting
SQ = CER
formals(SQ)$HCR_ICP = list(c(0, 10))           # Control points are PLL 1.5 and 2.5
formals(SQ)$HCR_LCP = list(c(1, 1))       # TAC is between 50% and 150% of 2022 catches
class(SQ) = 'mp'


jpeg('C:/GitHub/mahiMSE/docs/tsd/img/MP/Emp_MP.jpg', res=300, height=6.5, width=9, units="in")

  par(mfrow=c(2,3),mai=c(0.5,0.7,0.3,0.05))
  HCR_info(SQ,ymax = 4, xmax = 6)
  HCR_info(CER, leg=F, ymax = 4, xmax = 6)
  HCR_info(CER_A, leg=F, ymax = 4, xmax = 6)
  HCR_info(CER_B, leg=F, ymax = 4, xmax = 6)
  HCR_info(CER_C, leg=F, ymax = 4, xmax = 6)
  HCR_info(CER_S, leg=F, ymax = 4, xmax = 6)

dev.off()


#jpeg('C:/GitHub/mahiMSE/docs/tsd/img/MP/CER_reb_dep.jpg', res=300, height=6.5, width=6.5, units="in")

  par(mfrow=c(2,2),mai=c(0.5,0.7,0.3,0.05))
  CER(Example_Data_dep,plot=T,onlyHCR=F)
  CER(Example_Data_reb,plot=T,onlyHCR=F)

#dev.off()

# ---- Slick stuff ---------------------

library(Slick)
slick = readRDS("C:/GitHub/mahiMSE/docs/tsd/data/Slick/mahiSlick_1.rds")
App(slick=slick)


plotTimeseries(slick, includeQuants = F)
plotTimeseries(slick, includeQuants = T)

#plotTimeseries(slick, 2, includeQuants = F)
#plotTimeseries(slick, 3, includeQuants = T)



# === END =========================================================================




