# Data = Example_Data; HCR_use = T; HCR_ICP = list(c(0.5, 1.5), c(0.1, 0.6)); HCR_LCP = list(c(0.5,1.5), c(0.5, 1.5)); HCR_calib = 145:148;
# HCR_fleet = list(c(1,6)); HCR_sens = c(1,0.5); HCR_lever = c("TAC", "TAE");
# HCR_up_max = c(0.2,0.1); HCR_down_max = c(0.2,0.2); HCR_input = c(1,2); plot=T; onlyHCR = F
do_HCR = function(Data, HCR_use, HCR_ICP, HCR_LCP, HCR_fleet, HCR_sens, HCR_lever,
                  HCR_up_max, HCR_down_max, HCR_input,  HCR_smooth, Sim_Ind, plot, onlyHCR, ad){

  # Notes:
  # Relative TAC level (relative to 2022) goes in in Data@Misc$rel_lev_vals and out in ad@Misc$rel_lev_vals (for determining acceptible rate of change)

  d = get_dims(Data); nt = d$nt; na = d$na; nf = d$nf; nr = d$nr
  season = get_season(Data)

  if(any(HCR_use)){

    # set up blank TAC and TAE vectors by fleet, or alternatively borrow these if they have already been specified in some way
    if("TAC" %in% HCR_lever){
      TAC = rep(NA,d$nf)
      if(!is.null(ad@TAC))TAC=ad@TAC
    }

    if("TAE" %in% HCR_lever){
      TAE = rep(NA,d$nf)
      if(!is.null(ad@Effort))TAE=ad@Effort
    }

    nHCR = length(HCR_fleet)
    #oldAdvice = Data@Advice # advice from previous year including Data@Misc@rel_lev_vals
    old_rel_lev_vals = Data@Misc$rel_lev_vals
    if(is.null(old_rel_lev_vals))old_rel_lev_vals = rep(1, nHCR) # if first advice year
    new_rel_lev_vals = rep(NA,length(old_rel_lev_vals))
    for(hh in 1:nHCR){ if(HCR_use[hh]){

      Index = Sim_Ind[[HCR_input[hh]]]

      # fleet index  - which fleets are having their advice updated
      if(is.na(HCR_fleet[[hh]][1])){
        fleets = rep(T,nf)
      }else{
        fleets = rep(F, nf)
        fleets[HCR_fleet[[hh]]] = T
      }

      # HCR lever
      if(HCR_lever[hh] == "TAE"){
        hind = match(Data@YearLH, Data@Years)-3:0
        hEffort = Data@Effort@Value[hind,]
        reflev=sum(hEffort[season,fleets])
      }else if(HCR_lever[hh] =="TAC"){
        reflev = sum(rTAC[season,fleets,], na.rm=T)
      }

      if(season == 1){ # Do control calcs for season 1 - updates are annual so can skip this for seasons 2-4
        if(!is.na(HCR_smooth[hh])) Index = Index_smooth(Index, HCR_smooth[hh], plot = plot, plotname = paste0("HCR Index Smooth #",hh)) # polynomial smoother parameterized as 'effective number of parameters' - a fraction of the total length of the time series
        xCP = HCR_ICP[[hh]]
        yCP = reflev * HCR_LCP[[hh]]
        curx = Index[length(Index)]                       # last historical annual index lagged by specified year
        rel_lev_val = old_rel_lev_vals[hh]
        levlab = paste0(HCR_lever[hh],": ", paste(Fleets[fleets],collapse=","))
        new_rel_lev_vals[hh] = do_HockStick(curx, xCP, yCP, reflev, rel_lev_val, HCR_up_max_hh = HCR_up_max[hh], HCR_down_max_hh = HCR_down_max[hh], HCR_sens_hh = HCR_sens[hh], plot=plot, onlyHCR = onlyHCR, hh, levlab=levlab, Data=Data)
      }else{ # same as season 1 adjustment
        new_rel_lev_vals[hh] = old_rel_lev_vals[hh]
      }

      newlev = new_rel_lev_vals[hh] * reflev

      fleetnos=(1:nf)[fleets]

      if(HCR_lever[hh]=="TAE"){
        hind = match(Data@YearLH,Data@Years)-3:0
        hEffort = Data@Effort@Value[hind,]
        sTAE=hEffort[season,]
        fleetfracs = sTAE[fleets]/reflev
        TAE[fleetnos] = fleetfracs * newlev
      }else if(HCR_lever[hh] =="TAC"){
        sTAC=rTAC[season,,]
        fleetfracs = apply(sTAC[fleets,,drop=F]/reflev,1,sum)
        TAC[fleetnos] = fleetfracs * newlev
      }

    }} #  HCR use loop / HCR loop

    if("TAC" %in% HCR_lever){
      TAC[is.na(TAC)] = 1E9 # 1 mil tons (unrestricted)
      ad@TAC = TAC
    }
    if("TAE" %in% HCR_lever)ad@Effort = TAE

    #if(season < 4)  ad@Misc$rel_lev_vals[hh] = old_rel_lev_vals            # in-season, still referenced from previous year
    #if(season == 4) ad@Misc$rel_lev_vals[hh] = new_rel_lev_val             # record update for next year
    ad@Misc$rel_lev_vals = new_rel_lev_vals

  } # use HCR condition

   ad
}

do_HockStick = function(curx, xCP, yCP, reflev, rel_lev_val, HCR_up_max_hh, HCR_down_max_hh, HCR_sens_hh, plot=F, onlyHCR=F, hh, levlab="Management Lever",Data){

  if(curx<xCP[1]){
    val = yCP[1]
  } else if(curx>xCP[2]){
    val = yCP[2]
  } else {
    val = yCP[1] + (yCP[2]-yCP[1])* (curx-xCP[1])/(xCP[2]-xCP[1])
  }
  ref_update = val/reflev                          # compared to historical ref level
  prop_change = ref_update/rel_lev_val             # size of update relative to last change
  sens_change = exp(log(prop_change)*HCR_sens_hh) # apply sensitivity control parameter
  update = update2 = sens_change
  if((update-1)>HCR_up_max_hh) update2 = 1 + HCR_up_max_hh
  if((update-1)< (-HCR_down_max_hh)) update2 = 1 - HCR_down_max_hh
  delta = update2 * rel_lev_val

  if(plot|onlyHCR){
    levlab2 = ""
    if(nchar(levlab)>24){
      sp =  strsplit(levlab,",")[[1]]
      levlab2 = paste0(sp[floor(length(sp)/2):length(sp)],collapse=", ")
      levlab = paste0(sp[1:(floor(length(sp)/2)-1)],collapse=", ")
    }
    plot(c(0,max(xCP,curx)*1.2),c(0,max(yCP,reflev)*1.1),col="white",xlab="Index",ylab="");grid()
    mtext(levlab,2,line=2.8,cex=0.7)
    mtext(levlab2,2,line=1.9,cex=0.7)
    lines(c(0,xCP[1]),rep(yCP[1],2),lwd=4); lines(xCP,yCP,lwd=4); lines(c(xCP[2],xCP[2]*1000),rep(yCP[2],2),lwd=4)
    abline(v=curx,col="blue",lwd=2)
    abline(h=reflev,col="orange",lwd=2)
    abline(h=val,col="blue",lwd=2)
    abline(h=rel_lev_val*reflev,col="grey",lwd=2)
    abline(h=update*reflev*rel_lev_val,col="green",lty=2,lwd=1.5)
    abline(h=update2*reflev*rel_lev_val,col="red", lty=2, lwd=2)

      legend('topleft',legend=c(
        paste0("Reference (",round(reflev,2),")"),
        paste0("Previous (",round(rel_lev_val,3),")"),
        paste0("With control rule (",round(ref_update,3),")"),
        paste0("With sensitivity (",round(update*rel_lev_val,3),")"),
        paste0("With change constraints (",round(update2*rel_lev_val,3),")")),
        cex=0.7,text.col=c("orange","grey","blue","green","red"),bty="n")
      xs = max(xCP,curx)*1.2 * c(0.96,0.92,0.88)
      acond = function(x,y)(abs(1-x/y)>0.025)
      if(acond(rel_lev_val*reflev,val))  arrows(xs[1],rel_lev_val*reflev,xs[1],val,0.15,col="blue",lwd=2)
      if(val>1E5){if(acond(val, update*reflev*rel_lev_val)){ arrows(xs[2],val,xs[2], update*reflev*rel_lev_val,0.15,col="green",lwd=2)}}
      if(acond(update*reflev*rel_lev_val, update2*reflev*rel_lev_val)) arrows(xs[3],update*reflev*rel_lev_val,xs[3],update2*reflev*rel_lev_val,0.15,col="red",lwd=2)
      mtext(paste0("HCR #",hh,",  ",Time_steps[length(Data@Years)]),line=0.3,font=2,cex=0.9)
  }

  delta
}

gethistind = function(Data, ns=4, lag=2)  match(Data@YearLH, Data@Years) - (lag-1)*ns - ns:1 + 1

mapseason = function(Data, season, ns = 4, lag=1){ # seasonal time steps of the last 4 seasons of a Dataset with the lag specified year
  nt = length(Data@Years)
  (nt-season-ns*(lag-1)+1)-(ns-1):0
}


mapseason_lastyear = function(Data, season, ns = 4, lag=1){ # final season time step of the lag specified year
  nt = length(Data@Years)
  (nt-season-ns*(lag-1)+1)
}

mapseason_index = function(Index, season, ns=4, lag=1){ # seasonal time steps for the last 4 seasons of a time series with the lag specified year (lags often build into indices already)
  nt = length(Index)
  (nt-season-ns*(lag-1)+1)-(ns-1):0
}


Index_smooth<-function(xx, enp.mult, plot=F, plotname=""){
  tofill<-!is.na(xx)
  xx[xx==0]<-1E3
  predout<-rep(NA,length(xx))
  dat<-data.frame(x=1:length(xx),y=log(xx))
  enp.target<-sum(tofill)*enp.mult
  out<-loess(y~x,dat=dat,enp.target=enp.target)
  predout[tofill]<-exp(predict(out))
  if(plot){
    ts_ind = (1:length(xx))[!is.na(xx)]
    plot(xx,type="p",xlab="Year",ylab="Index",ylim=c(0,max(xx,na.rm=T)*1.05),main="",col="black");grid()
    legend('topright',legend=c("Index",paste0("Smoothed (",enp.mult,")")),text.col=c("black","red"),bty='n')
    lines((1:length(xx))[!is.na(xx)],predout[ts_ind],col="#ff000090",lwd=2)
    mtext(plotname,line=0.3,cex=0.9,font=2)
  }
  predout
}


