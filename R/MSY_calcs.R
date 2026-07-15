# main function
calcmahi_MSY = function(MSEobj,stock=1){

  dms = getMSEdims(list(MSEobj))
  MSYin = makeMSYinputs(MSEobj, dms)
  logFMSY = sapply(1:dms$nsim, doF_MSYcalc, MSYin = MSYin)
  MSYin$logF_s = logFMSY
  MSYout = sapply(1:dms$nsim,doMSYcalc, MSYin)
  eqSB0 = MSEobj@Unfished@Equilibrium@SBiomass[,stock,dms$ny]
  eqB0 = MSEobj@Unfished@Equilibrium@Biomass[,stock,dms$ny]
  eqSBMSY = eqSB0 * MSYout[4,]
  UB = MSYout[1,] / MSYout[6,]
  UVB = MSYout[1,] / MSYout[7,]
  MSYout = rbind(MSYout, eqSB0, eqB0, eqSBMSY, UB, UVB)
  dimnames(MSYout)[[2]] = paste0("sim", 1:dms$nsim)
  as.data.frame(t(MSYout))

}

calcSSB0 = function(OM){

  OMp = OM
  OMp@Stock$Dolphinfish@SRR@RecDevProj[] = exp(log(OMp@Stock$Dolphinfish@SRR@RecDevProj[])*1E-4)
  hist=Simulate(OMp)

  dUB = hist@Unfished@Dynamic@Biomass[,1,]
  dUSSB = hist@Unfished@Dynamic@SBiomass[,1,]
  eUB = hist@Unfished@Equilibrium@Biomass[,1,]
  eUSSB = hist@Unfished@Equilibrium@SBiomass[,1,]

  nofish = function(Data){
   Advice(TAC=1E-5)
  }
  class(nofish) = 'mp'

  proj = Project(hist,"nofish",nSim=4)

  pB = proj@Biomass[1:6,1,77:80,1]
  # eUB[1:6,145:148]
  #apply(dUB,1,mean)
  #dUB[1:6,145:148]

  pSSB = proj@SBiomass[1:4,1,77:80,1]
  #eUSSB[1:6,145:148]
  #dUSSB[1:6,145:148]

}

doMSYcalc=function(x, MSYin){
  out = MSEtool:::MSYCalcs(
    logF = MSYin$logF_s[x],
    M_at_Age = MSYin$M_at_Age_s[x,],
    Wt_at_Age = MSYin$Wt_at_Age_s[x,],
    Mat_at_Age = MSYin$Mat_at_Age_s[x,],
    Fec_at_Age = MSYin$Fec_at_Age_s[x,],
    V_at_Age = MSYin$V_at_Age_s[x,],
    Wt_at_Age_C = MSYin$Wt_at_Age_s[x,],
    maxage = MSYin$maxage,
    relRfun = function(){},
    SRRpars=list(),
    R0x = MSYin$R0x_s[x],
    SRrelx = 1,
    hx = MSYin$hx_s[x],
    SSBpR = MSYin$SSBpR_s[x],
    opt=2L)
  out

}


doF_MSYcalc=function(x, MSYin, boundsF = c(0.1,1.5)){

  logF = MSYin$logF_s[x]
  M_at_Age = MSYin$M_at_Age_s[x,]
  Wt_at_Age = MSYin$Wt_at_Age_s[x,]
  Mat_at_Age = MSYin$Mat_at_Age_s[x,]
  Fec_at_Age = MSYin$Fec_at_Age_s[x,]
  V_at_Age = MSYin$V_at_Age_s[x,]
  Wt_at_Age_C = MSYin$Wt_at_Age_s[x,]
  maxage = MSYin$maxage
  relRfun = function(){}
  SRRpars=list()
  R0x = MSYin$R0x_s[x]
  SRrelx = 1
  hx = MSYin$hx_s[x]
  SSBpR = MSYin$SSBpR_s[x]

  doopt <- optimise(MSEtool:::MSYCalcs, log(boundsF),
                    M_at_Age,
                    Wt_at_Age,
                    Mat_at_Age,
                    Fec_at_Age,
                    V_at_Age,
                    Wt_at_Age_C,
                    maxage,
                    relRfun,
                    SRRpars,
                    R0x,
                    SRrelx,
                    hx,
                    SSBpR,
                    opt=1)
  doopt$minimum

}


makeMSYinputs = function(MSEobj, dms, ns = 4, stock = 1){

  muF_s = apply(MSEobj@Hist@Effort[,dms$nt-0:(ns-1),],c(1,3),mean)
  logF_s = rep(-0.5,dms$nsim) # log(apply(muF,1,sum))
  M_at_Age_s = MSEobj@OM@Stock[[stock]]@NaturalMortality@MeanAtAge[,,dms$nt]
  Wt_at_Age_s = MSEobj@OM@Stock[[stock]]@Weight@MeanAtAge[,,dms$nt]
  Mat_at_Age_s = MSEobj@OM@Stock[[stock]]@Maturity@MeanAtAge[,,dms$nt]
  Fec_at_Age_s = MSEobj@OM@Stock[[stock]]@Fecundity@MeanAtAge[,,dms$nt]
  allV = aperm(array(unlist(lapply(MSEobj@OM@Fleet[[stock]],function(x,dms)x@Selectivity@MeanAtAge[,,dms$nt,1],dms=dms)),c(dms$nsim,dms$na,dms$nf)),c(1,3,2)) # sim, fleet, age #area 1
  CbF = MSEobj@Hist@Landings[,stock,dms$nt,]
  totF = apply(CbF,1,sum)
  wV = apply(allV * array(CbF,dim(allV)),c(1,3),sum)
  V_at_Age_s = wV / apply(wV,1,max)
  surv_s = t(apply(M_at_Age_s,1,function(x)exp(cumsum(-x))))

  rec_s = apply(MSEobj@Hist@Number$Dolphinfish[,1,,],1:2,sum)
  SSB_s = MSEobj@Hist@SBiomass[,stock,]
  hx_s = MSEobj@OM@Stock[[stock]]@SRR@Pars$h[,dms$nt]
  if(length(hx_s)<dms$nsim)hx_s = rep(hx_s,dms$nsim)
  SSBpR_s = apply(surv_s * Mat_at_Age_s * Wt_at_Age_s, 1, sum)
  R0x_s = rep(1,dms$nsim) #MSEobj@OM@Stock[[stock]]@SRR@R0[,dms$nt]
  # R0_guess_s =apply(rec_s,1,mean)*1.25
  # R0x_season = sapply(1:dms$nsim, R0solve, rec_s=rec_s, SSB_s=SSB_s, hx_s=hx_s, SSBpR_s=SSBpR_s, R0_guess_s=R0_guess_s)
  # R0x_s = R0x_season * ns
  maxage = dms$na-1

  MSYin = list(logF_s = logF_s, M_at_Age_s = M_at_Age_s, Wt_at_Age_s= Wt_at_Age_s,
               Mat_at_Age_s = Mat_at_Age_s, Fec_at_Age_s=Fec_at_Age_s, V_at_Age_s=V_at_Age_s,
               R0x_s=R0x_s, hx_s=hx_s, surv_s=surv_s, SSBpR_s = SSBpR_s, maxage=maxage)
  MSYin
}


R0solve = function(x, rec_s, SSB_s, hx_s, SSBpR_s, R0_guess_s){

  rec = rec_s[x,]
  SSB = SSB_s[x,]
  hx = hx_s[x]
  SSBpR = SSBpR_s[x]
  R0_guess = R0_guess_s[x]
  opt<-optim(log(R0_guess),getBH,method="L-BFGS-B",lower=log(R0_guess/2),upper=R0_guess*2, hessian=F,
             SSB=SSB,rec=rec,SSBpR=SSBpR,mode=1, plot=T,h=hx)
  exp(opt$par)

}

getBH<-function(pars,SSB,rec,SSBpR,mode=1,plot=F,h){

  R0 = exp(pars)

  recpred<-((0.8*R0*h*SSB)/(0.2*SSBpR*R0*(1-h)+(h-0.2)*SSB))

  if(plot){
    recentrec<-rec[length(rec)]
    ord<-order(SSB)
    cols=rep("#0000ff95",length(rec))
    cols[length(rec)]<-"#ff000095"
    plot(SSB[ord],rec[ord],ylim=c(0,max(rec,R0)),xlim=c(0,max(SSB,R0*SSBpR)),xlab="",ylab="",col=cols[ord],pch=19)
    SSB2<-seq(0,max(R0*SSBpR,max(SSB)),length.out=500)
    recpred2<-((0.8*R0*h*SSB2)/(0.2*SSBpR*R0*(1-h)+(h-0.2)*SSB2))
    lines(SSB2,recpred2,col='blue')
    abline(v=c(0.2*R0*SSBpR,R0*SSBpR),lty=2,col='red')
    abline(h=c(R0,R0*h),lty=2,col='red')
    legend('topright',legend=c(paste0("h = ",round(h,3)),paste0("lnR0 = ",round(log(R0),3))),bty='n')
  }

  if(mode==1){
    #return(sum(((recpred-rec)/10000)^2))
    return(-sum(dnorm(log(recpred)-log(rec),0,0.5,log=T))-dnorm(pars[1],0,2,log=T)) # add a vague prior on h = 0.8
    #return(-sum(dnorm(recpred,rec,rec*0.5,log=T)))
  }else{
    return(rec-recpred)
  }

}

