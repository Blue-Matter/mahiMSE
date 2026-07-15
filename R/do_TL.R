do_TL = function(Data, TL, plot=F, TL_empirical, ad, ref_ts = 20, nseasons = 4){ # ref_ts is the last 5 years of data used to characterize CR distributions in 6_Dolphin_MP_prereq.R

  if(any(!is.na(TL))){ # only if something has to be calculated

    nage = nAge(Data@Misc$DataOM@OM)
    #Ret0 = Retention(MeanAtAge = rep(1,nage))
    #Retention_list = rep(list(Ret0),nFleet(Data)) #
    Retention_list = ad@Retention

    season = get_season(Data)
    find = (1:length(TL))[!is.na(TL)]
    refind = (match(Data@YearLH, Data@Years)-((ref_ts-1):0))[rep(1:nseasons,ref_ts/nseasons)==season] # the correct season time step over last ref_ts timesteps
    dat_ts = length(Data@Years)
    narea = dim(C_yfr)[3]

    ret = list()
    j = 0

    for(ff in find){

      VB = N = Data@Misc$DataOM@Number$Dolphinfish[1,,,] # nage, nts, narea
      Wtage = Data@Misc$DataOM@OM@Stock$Dolphinfish@Weight@MeanAtAge[1,,] # age, nts
      Selage = Data@Misc$DataOM@OM@Fleet[[1]][[ff]]@Selectivity@MeanAtAge[1,,,] # age, nts
      Nind = MSEtool:::TEG(dim(N))
      VB[Nind] = N[Nind] * Wtage[Nind[,1:2]] * Wtage[Nind[,1:2]]
      VBr = apply(VB,2:3,sum)
      areas = (1:narea)[apply(C_yfr[refind,ff,],2,sum)>0]

      for(rr in areas){
        j=j+1
        VBratio = VBr[dat_ts+1,rr] / mean(VBr[refind,rr])
        if(TL_empirical) retval = calc_ret_emp(VBratio, season, ff, rr, TL[ff], plot = plot)
        if(!TL_empirical) retval = calc_ret_fit(VBratio, season, ff, rr, TL[ff], plot=plot)
        ret[[j]] = c(season=season,ff=ff,rr=rr,retval=retval)
      }

    }

    # now overwrite where applicable
    for(j in 1:length(ret)){
      ff = ret[[j]]['ff']
      rr = ret[[j]]['rr']
      retval = ret[[j]]['retval']
      if(!is.na(retval)){
        Retention_list[[ff]]@MeanAtLength[,rr] = Retention_list[[ff]]@MeanAtLength[,rr] * retval
      }
    }

    ad@Retention = Retention_list

  } # end of if any(!is.na(TAC))

  ad

}


calc_ret_emp=function(VBratio, season, ff, rr, TLff, plot=F){
  #  VBratio = 0.5742388; ff = 3; rr = 1; TLff = 10
  retffrr = NA
  lu = TL_info$lookup
  nTLdat = nrow(lu)
  ind = lu$Q == season & lu$Ft == ff & lu$A == rr

  if(any(ind)){ # only carry on if bag limit data are available
    TLind = (1:nTLdat)[ind]
    edens = TL_info$emp_dens[[TLind]]
    xs = xbl = edens$x * VBratio # shrink / expand x axis to get expected catch rate (keeps CV the same, adjusts mean only)
    totdens = sum(xs * edens$y) # assuming all is kept
    xbl[xbl>TLff] = TLff
    retdens = sum(xbl* edens$y)
    retffrr = retdens / totdens
    if(plot){
      plot(edens,type="l",col="darkgreen", main=paste0("Tp. Lim. Ret. Q = ", season, ", Ft = ",ff,", Area = ",rr))
      lines(xs,edens$y,col="red"); abline(v=TLff,col="blue")
      legend('topright',legend=c("2018-2022",
                                 "In Projection",
                                 paste0("Trip lim. = ", TLff),
                                 paste0("Retn. = ",round(retffrr,3)),
                                 paste0("VBratio =", round(VBratio,3))),
                                  text.col =c("darkgreen","red","blue","black","black"))
    }
  }

  retffrr

}



calc_ret_fit=function(VBratio, season, ff, rr, TLff, maxCR=100, plot=F){
  #  VBratio = 0.5742388; ff = 3; rr = 1; TLff = 10
  retffrr = NA
  lu = TL_info$lookup
  nTLdat = nrow(lu)
  ind = lu$Q == season & lu$Ft == ff & lu$A == rr

  if(any(ind)){ # only carry on if bag limit data are available
    TLind = (1:nTLdat)[ind]
    size = TL_info$fitsize[[TLind]]
    mu = TL_info$fitsize[[TLind]] * VBratio
    CRs = CRsbl = 1:maxCR
    dens = dnbinom(CRs,size=size, mu=mu)
    totdens = sum(CRs * dens) # assuming all is kept
    CRsbl[CRsbl > TLff] = TLff
    retdens = sum(CRsbl* dens)
    retffrr = retdens / totdens
    if(plot){
      plot(CRs,dens,type="l",main=paste0("Q = ", season, ", Ft = ",ff,", Area = ",rr))
      lines(CRsbl,dens,col="red"); abline(v=TLff,col="blue")
      legend('topright',legend=c(paste0("Trip lim. = ", TLff),
                                 paste0("Retn. = ",round(retffrr,3)),
                                 paste0("VBratio =", round(VBratio,3))))
    }
  }

  retffrr

}
