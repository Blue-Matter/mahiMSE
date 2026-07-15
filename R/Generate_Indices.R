
# convert a standard normal distribution (mean 0, sd 1) to a log-normal with mean 1 and cv specified
conv_snd = function(nd, cv) exp((nd * MSEtool::sdconv(mu,cv)) + MSEtool::mconv(mu,cv))

# OMlab = Data@Misc$DataOM@OM@Name
getOMno = function(OMlab){
  firstext = strsplit(OMlab,"#")[[1]][2]
  sectext = strsplit(firstext,":")[[1]][1]
  as.integer(sectext)
}

# Ipred = hvar; Iobs = Ind; betamin = -1; betamax = 1; nbeta = 100

getBeta = function(Iobs, Ipred, betamin = -1, betamax = 1, nbeta = 100){
  nts = length(Ipred)
  betavec = exp(seq(betamin,betamax,length.out=nbeta))
  getobj=function(beta, Iobs, Ipred){sum((MSEtool:::CalcIndexResiduals(Iobs, matrix(Ipred^beta,nrow=1))$LogResiduals)^2)}
  #getobj=function(beta, Iobs, Ipred){cor(Iobs,Ipred^beta)^2}

  ressdvec = sapply(betavec,getobj,Iobs=Iobs, Ipred=Ipred)
  betavec[which.min(ressdvec)]
  # plot(betavec,ressdvec)
}


# Data = Example_Data_3; Index_number = 1; Index_cv = rep(NA,4); Index_ac = rep(NA,4); Index_lag = c(1,2,1,2); Index_seed = rep(1,4); Index_beta = c(F,F,F,F); return_stats=F

gen_ind_2 = function(Data, Index_names, Index_cv, Index_ac, Index_lag, Index_seed, Index_beta, plot=T, plot_ts = 100, refyrs = 4, return_stats=F){

  if(!is.na(Index_names[1])){ # only if indices are needed
    # Get season and pertinent MSE dimensions
    season = get_season(Data)
    lhy = match(Data@YearLH, Data@Years) # last historical time step
    FDeadArea = Data@Misc$DataOM@FDeadArea$Dolphinfish[1,,lhy,,]
    na = dim(FDeadArea)[1]; nf = dim(FDeadArea)[2]; nr = dim(FDeadArea)[3]
    N = Data@Misc$DataOM@Number$Dolphinfish[1,,,]
    nt = dim(N)[2]  # Total time steps (historical and projected)
    lhobs = length(Data@Years) # last observation

    # Get index-specific info
    subIdat = Indices[Indices$Name %in% Index_names,]
    INames = Index_names
    rowind = match(INames, subIdat$Name)
    INo = subIdat$Number[rowind]
    ISeason = subIdat$Season[rowind]
    IArea = subIdat$Area[rowind]
    IType = subIdat$Type[rowind]
    nI = length(INames)

    # For Random number generation (must be consistent across projected years)
    OMno = getOMno(Data@Misc$DataOM@OM@Name)
    lhrec = N[1, match(Data@YearLH,Data@Years),1] # rec in last historical year for area 1 used for random seeding

    sIndices = stats = obsyrs = obs = allobs = allobsyrs = pred = allpred = resid = fullyrs = full = list() # Recording for plots

    for(ii in 1:nI){

      if(IType[ii] == "B")   var = Data@Misc$DataOM@Biomass[1,1,] /1E6
      if(IType[ii] == "SSB") var = Data@Misc$DataOM@SBiomass[1,1,]/1E6
      if(grepl("VB",IType[ii])){
        ft = as.integer(strsplit(IType[ii],"VB")[[1]][2])
        sel = Data@Misc$DataOM@OM@Fleet[[1]][[ft]]@Selectivity@MeanAtAge[1,,,] # age, ts, area
        wta = Data@Misc$DataOM@OM@Stock$Dolphinfish@Weight@MeanAtAge[1,,]      # age, ts
        I_area_vec = rep(0,nr); I_area_vec[IArea[[ii]]] = 1
        VBA = N * array(sel, c(na, nt, nr)) * array(rep(I_area_vec, each = na * nt), c(na, nt, nr)) * array(wta, c(na, nt, nr))
        var = apply(VBA,2, sum) / 1E6
      }

      # Statistics
      Ind_yrs = Indices$Year[Indices$Name == INames[ii]]

      if(ISeason[ii] == 0){ # all seasons
        var_yrs = floor(Time_steps[1:lhy])
        keep = (1:lhy)[var_yrs %in% Ind_yrs]
        hvar = aggregate(var[keep],by=list(yr = var_yrs[keep]),mean)$x
      }else{               # A specific season specified
        Ind_yrs_s = Ind_yrs + (ISeason[ii]/4)-1/4
        keep = match(Ind_yrs_s,Time_steps)
        keep = keep[keep<=lhy] # only historical years for stat calcs
        hvar = var[keep]
      }

      # Calc errors based on historical model correspondence (ie observed indices and OM quantity to 2022)
      Ind0 = Indices$Index[Indices$Name == INames[ii]]
      Ind = Ind0[Ind_yrs <= Data@YearLH]

      if(Index_beta[ii]) beta = getBeta(Ind,hvar)
      if(!Index_beta[ii]) beta = 1
      hvar = hvar^beta
      res = MSEtool:::CalcIndexResiduals(Ind/mean(Ind),matrix(hvar/mean(hvar),nrow=1))$LogResiduals[1,]
      muhInd = mean(Ind)
      muhvar = mean(hvar)

      cop = muhInd/muhvar
      df = MSEtool:::Calc_Stats(res)
      df$beta = beta

      if(!is.na(Index_cv[ii]))df$SD = Index_cv[ii] # Manual override
      if(!is.na(Index_ac[ii]))df$AC = Index_ac[ii] # Manual override

      # Projection years to simulate data
      nhy = lhy / 4
      ny = nt / 4
      sI = rep(NA,ny) # full time series
      OMyrs = unique(floor(Time_steps))
      sI[match(Ind_yrs,OMyrs)] = Ind0                # all observed historical data
      pind = (match(max(Ind_yrs),OMyrs)+1):ny       # projection years
      np = length(pind)                             # number of projection years

      # get OM variable for all years
      if(ISeason[ii] == 0){ # all seasons
        var_yrs = floor(Time_steps)
        pvar = aggregate(var,by=list(yr = var_yrs),mean)$x ^ beta
      }else{               # A specific season specified
        var_yrs_s = OMyrs + (ISeason[ii]/4)-1/4
        keep = match(var_yrs_s,Time_steps)
        pvar = var[keep] ^ beta
      }

      # get lst.err - to correctly project errors with lag-1 AC you need the last log deviation which may be in the future if observed indices (e.g. to 2024) are longer than the historical time period (2022)
      dfproj = df
      lastobsy = pvar[pind[1]-1]
      lastobs = Ind0[length(Ind0)]
      dfproj$lst.err = log(lastobs/muhInd)-log(lastobsy/muhvar)

      set.seed(lhrec+Index_seed[ii]+OMno+ii) # future simulations differ according to simulation and can be shared among Indices (for testing rules with same data)
      perr = MSEtool:::Gen_Residuals(dfproj, 1, np)[1,]

      sI[pind] = pvar[pind] * perr * cop

      # apply lag to final available index
      currentyear = floor(Time_steps[length(Data@Years)])
      lagkeep = OMyrs <= currentyear-Index_lag[ii]
      sIndices[[ii]] = sI[lagkeep]
      names(sIndices[[ii]]) = OMyrs[OMyrs <= currentyear-Index_lag[ii]]

      if(plot|return_stats){ full[[ii]] =  sIndices[[ii]]; stats[[ii]] = dfproj;
                resid[[ii]] = res; obs[[ii]] = Ind; allobs[[ii]] = Ind0;
                allobsyrs[[ii]] = Ind_yrs; pred[[ii]] = hvar*cop;
                obsyrs[[ii]] = Ind_yrs[Ind_yrs <= Data@YearLH];
                fullyrs[[ii]] = OMyrs[OMyrs <= currentyear-Index_lag[ii]];
                allpred[[ii]] = pvar[1:max(pind)]*cop} # recording for plots

    }
    names(sIndices) = INames

    #  jpeg("C:/Users/tcar_/Dropbox/temp/Dolphin/comms/Index_first_look.jpg",res=400,width=12,height=10,units="in")
    if(plot){
      par(mfrow=c(nI,4),mai=c(0.5,0.5,0.25,0.1))
      xline=2.2;yline=2.2;labcex = 0.7
      for(ii in 1:nI){
        plot(obsyrs[[ii]],obs[[ii]],ylim=c(0,max(obs[[ii]])),xlab="",ylab="",pch=19);grid()
        lines(obsyrs[[ii]],pred[[ii]],col="red")
        mtext("Year",1,line=xline,cex=labcex); mtext("Index",2,line=yline,cex=labcex)
        legend('topright',legend=c("Obs.","Pred."),text.col=c("black","red"),bty="n")
        mtext(INames[[ii]],line=0.3, font=2,cex=0.9)

        plot(obs[[ii]],pred[[ii]],xlab="",ylab="",pch=19,xlim=c(0,max(obs[[ii]])),ylim=c(0,max(pred[[ii]])));grid();
        lines(c(-1E10,1E10),c(-1E10,1E10),col="blue",lty=2)
        mtext("Observed",1,line=xline,cex=labcex); mtext("Predicted",2,line=yline,cex=labcex)
        legend('bottomright',legend=c(paste0("CV = ",round(stats[[ii]]$SD,3)),paste0("AC = ",round(stats[[ii]]$AC,3)),paste0("Beta = ",round(stats[[ii]]$beta,3))),bty="n")


        plot(fullyrs[[ii]],full[[ii]],xlab="",ylab="",col="red",pch=19);grid()
        points(allobsyrs[[ii]],allobs[[ii]],col="black",pch=19,cex=1.05)
        ap = allpred[[ii]]
        yrp = OMyrs[1:length(ap)]
        lines(yrp,ap,col="red")
        mtext("Year",1,line=xline,cex=labcex); mtext("Index",2,line=yline,cex=labcex)
        legend('topright',legend=c("Obs.","Sim."),text.col=c("black","red"),bty="n")
        abline(v=max(obsyrs[[ii]]+0.5),col="blue",lty=2)

        plot(fullyrs[[ii]],log(full[[ii]]),xlab="",ylab="",col="red",pch=19);grid()
        points(allobsyrs[[ii]],log(allobs[[ii]]),col="black",pch=19,cex=1.05)
        lines(yrp,log(ap),col="red")
        mtext("Year",1,line=xline,cex=labcex); mtext("Log Index",2,line=yline,cex=labcex)
        legend('topright',legend=c("Obs.","Sim."),text.col=c("black","red"),bty="n")
        abline(v=max(obsyrs[[ii]]+0.5),col="blue",lty=2)

      }
    }
    #dev.off()
    if(return_stats){
      return(list(sIndices = sIndices, obsyrs = obsyrs, obs = obs, allobs = allobs,
                  allobsyrs = allobsyrs, pred = pred, allpred = allpred, resid = resid,
                  fullyrs = fullyrs, full = full)) # return stat

    }else{
      return(sIndices)
    }
  } # end of 'only if Index_number != NA
}


# Data = Example_Data_3; Index_type = c("VB1","SSB", "VAST"); Index_cv = c(0.01, 0.1, 0.4); Index_ac = c(0, 0.95, 0.95); Index_lag = c(1,2,3); Index_nt = c(4,8,4); Index_seed = rep(1,3); refyrs = 3
gen_ind = function(Data, Index_names, Index_cv, Index_ac, Index_lag, Index_nt, Index_seed, plot=T, plot_ts = 100, refyrs = 4){

  season = get_season(Data)
  nI = length(Index_names)
  N = Data@Misc$DataOM@Number$Dolphinfish[1,,,]
  nt = dim(N)[2]
  OMno = getOMno(Data@Misc$DataOM@OM@Name)
  lhrec = N[1, match(Data@YearLH,Data@Years),1] # rec in last historical year for area 1 used for random seeding
  lhobs = length(Data@Years)
  lhy = match(Data@YearLH, Data@Years)

  FDeadArea = Data@Misc$DataOM@FDeadArea$Dolphinfish[1,,lhy,,]
  for(i in 1:(refyrs-1)) FDeadArea = FDeadArea + Data@Misc$DataOM@FDeadArea$Dolphinfish[[lhy-i]]
  na = dim(FDeadArea)[1]; nf = dim(FDeadArea)[2]; nr = dim(FDeadArea)[3]

  Indices = list()

  for(ii in 1:nI){

    # Get the index to the current year - 1
    if(Index_type[ii] == "VAST"){

      ind = Data@Survey@Value[,1]

      if(!is.na(Index_cv[ii]) | !is.na(Index_ac[ii])){ # if user specified stats properties

        wta = Data@Misc$DataOM@OM@Stock$Dolphinfish@Weight@MeanAtAge
        VBA = N * array(wta, c(na, nt, nr))
        var = apply(VBA, 2, sum) / 1E6
        hvar = var[1:lhy]
        muvar = mean(hvar)
        hind = ind[1:lhy] # historical VAST observations
        muind = mean(hind)
        res = log(hind/muind)-log(hvar/mean(hvar))
        df = MSEtool:::Calc_Stats(res)
        if(!is.na(Index_cv[ii]))df$SD = Index_cv[ii]
        if(!is.na(Index_ac[ii]))df$AC = Index_ac[ii]
        set.seed(lhrec+Index_seed[ii]+OMno) # future simulations differ according to simulation and can be shared among Indices (for testing rules with same data)
        perr = MSEtool:::Gen_Residuals(df, 1, lhobs-lhy)[1,]
        pvar = var[(lhy+1):lhobs] / muvar * muind
        ind = c(hind,pvar*perr)# Joined
      }

    }else if(grepl("sVAST",Index_type[ii])){ # Area-based vast index

      rr = as.integer(strsplit(Index_type[ii],"sVAST")[[1]][2]) # area
      wta = Data@Misc$DataOM@OM@Stock$Dolphinfish@Weight@MeanAtAge
      VBA = N * array(wta, c(na, nt, nr))
      var = apply(VBA[,,rr], 2, sum) / 1E6
      hvar = var[1:lhy]
      muvar = mean(hvar)
      hind = sVAST[,rr]
      muind = mean(hind)
      res = log(hind/muind)-log(hvar/mean(hvar))
      df = MSEtool:::Calc_Stats(res)
      if(!is.na(Index_cv[ii]))df$SD = Index_cv[ii]
      if(!is.na(Index_ac[ii]))df$AC = Index_ac[ii]
      set.seed(lhrec+Index_seed[ii]+OMno) # future simulations differ according to simulation and can be shared among Indices (for testing rules with same data)
      perr = MSEtool:::Gen_Residuals(df, 1, lhobs-lhy)[1,]
      pvar = var[(lhy+1):lhobs] / muvar * muind
      ind = c(hind,pvar*perr)# Joined

    }else{

      if(Index_type[ii] == "B")   var = Data@Misc$DataOM@Biomass[1,1,] /1E6
      if(Index_type[ii] == "SSB") var = Data@Misc$DataOM@SBiomass[1,1,]/1E6
      if(grepl("VB",Index_type[ii])){
        ft = as.integer(strsplit(Index_type[ii],"VB")[[1]][2])
        sel = Data@Misc$DataOM@OM@Fleet[[1]][[ft]]@Selectivity@MeanAtAge[1,,,] # age, ts, area
        wta = Data@Misc$DataOM@OM@Stock$Dolphinfish@Weight@MeanAtAge[1,,]      # age, ts
        F_area_last = as.integer(FDeadArea[na,ft,]>1E-2) # what areas are fished by this fleet?
        VBA = N * array(sel, c(na, nt, nr)) * array(rep(F_area_last, each = na * nt), c(na, nt, nr)) * array(wta, c(na, nt, nr))
        var = apply(VBA,2, sum) / 1E6
        hvar = var[1:lhy]
        muvar = mean(hvar)
        hind = ind[1:lhy] # historical VAST observations
        muind = mean(hind)
        res = log(hind/muind)-log(hvar/mean(hvar))
        df = MSEtool:::Calc_Stats(res)
      }

      # Calc errors
      # Historical errors
      set.seed(OMno) # historical errors are the same for all simulations of OM (ie real observations)
      df = data.frame(SD=Index_cv[ii], AC=Index_ac[ii], lst.err=0)
      herr = MSEtool:::Gen_Residuals(df, 1, lhy)[1,]
      # Projected errors
      df$lst.err = log(herr[lhy])
      set.seed(lhrec+Index_seed[ii]+OMno) # future simulations differ according to simulation and can be shared among Indices (for testing rules with same data)
      perr = MSEtool:::Gen_Residuals(df, 1, nt-lhy)[1,]
      # Joined
      ind = (var * c(herr, perr))[1:lhobs]

    }

    # apply lag to final available index
    Indices[[ii]] = ind[1:mapseason_lastyear(Data, season, ns=4, lag=Index_lag[ii])]

  }


  if(plot){
    # par(mfrow=c(ceiling(nI/2),2),mai=c(0.8,0.8,0.5,0.1))

    for(ii in 1:nI){
      ts_ind = (nYear(Data@Misc$DataOM@OM) - (plot_ts-1)):length(Indices[[ii]])
      plot(Time_steps[ts_ind],Indices[[ii]][ts_ind],pch=19,xlab="Time step",ylab=Index_type[ii],ylim=c(0,max(Indices[[ii]])*1.05)); grid()
      abline(v=lhy+0.5,col="red")
      mtext(paste0("Index Sim. (",ii,") CV = ", Index_cv[ii], " AC = ", Index_ac[ii]),cex=0.9,line=0.3, font=2)
    }
  }

  Indices

}
