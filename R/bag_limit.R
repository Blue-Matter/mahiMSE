# bag limit analysis and code


RRposthoc = function(dat,TL=c(54, 50, 45, 40, 35, 30)){

  RRdat = NULL
  totcat = dat$Catch
  for(i in 1:length(TL)){
    retcat = dat$Catch; retcat[retcat>TL[i]] = TL[i]
    tot = aggregate(totcat,by=list(dat$Year, dat$Area),sum)
    ret = aggregate(retcat,by=list(dat$Year, dat$Area),sum)
    rel = (1-ret$x/tot$x) * 100
    tempdat = data.frame(Year = tot$Group.1, State = tot$Group.2, TripLim = rep(TL[i],nrow(tot)), RR = rel)
    RRdat = rbind(RRdat, tempdat)
  }

  RRdat$TripLim = as.factor(RRdat$TripLim)

  ggplot(RRdat, aes(x=Year, y = RR, group=TripLim,colour=TripLim))+
    geom_line()+ facet_wrap(~ State)+
    ylab("Predicted Release Rate (%)")

}

process_CR_RR = function(dir, Targ=F){ # process catch rate vs release rate

  if(grepl(".csv",dir)){
    d4d = read.csv(dir)
  }else{
    d4d = as.data.frame(read_xlsx(dir))
  }

  dat = data.frame(Year = d4d$YEAR, Month=d4d$month,
                   Mode = d4d$MODE_F,
                   Gear = d4d$GEAR, Area = d4d$ST,
                   Targ = d4d$prim1_common, Targ_alt = d4d$prim2_common,
                   Ret = d4d$CLAIM_UNADJ,
                   Rel = d4d$RELEASE_UNADJ,
                   Catch = d4d$tot_cat,
                   RR = d4d$RELEASE_UNADJ / (d4d$CLAIM_UNADJ+d4d$RELEASE_UNADJ)) #d4d$tot_cat)

  if(Targ)dat = subset(dat, dat$Targ == "DOLPHIN")


  plot(dat$Catch,dat$RR,pch=19,col="#99999999",xlim=c(0,60),xlab = "Catch per trip (MRIP 1981 - 2024)", ylab ="Releases per trip (MRIP 1981 - 2024)");grid(col="red",lty=2)
  for(i in 1:60)lines(1:100,i/(1:100),col="blue",lty=2)


  #ggplot(dat, aes(x=Catch,y=RR)) +
   # geom_point(color="black", fill="white",binwidth=1)+
    #geom_vline(xintercept = 10,linetype= "dashed", col = "black", size = 0.7)+
    #facet_wrap(~ Area,scales = "free")+xlim(0, 30)


  dat

}

#catint = ceiling(catchy)
#k = x = 0:50
#size = r = 5
#prob = p = 0.2
#mu = r*(1-p)/p
#var = r*(1-p)/(p*p)
#sd = var^0.5
#dd = dnbinom(k,r,p)
#plot(k,dd)
#dd2 = dnbinom(k,r, mu = mu)
#lines(k,dd2,col="red")


nbinfit = function(lnsize, lnmu, catchy){
  size=exp(lnsize); mu=exp(lnmu)
  catint = ceiling(catchy)
  sum(-dnbinom(catint,size=size,mu=mu,log=T))
}


nbinfit2 = function(pars, catchy){
  size=exp(pars[1]); mu=exp(pars[2])
  catint = ceiling(catchy)
  sum(-dnbinom(catint,size=size,mu=mu,log=T))
}

nbplot = function(size,mu,catchy,maxplot,col="red"){
  k = 1:maxplot
  dd = dnbinom(k,size=size,mu=mu)
  lines(k,dd,col=col)
}

calc_emp = function(dens, demoBL = 54){ # calculates probability of being above catch rate
  ix = floor(dens$x)       # integer catch rate (bag limit, trip limit)

  # probability trip over
  sumy = aggregate(dens$y,by=list(ix),sum) # sum of empirical density for each node
  CRvals = sumy$Group.1
  csumy = cumsum(rev(sumy$x))              # cumulative sum in reverse (amount above CR)
  prob = rev(csumy/max(csumy))             # phrased as probability of being above CR
  prob = prob / max(prob)

  # fraction of fish
  sumf = sum(dens$y*dens$x) # sum of empirical density for each node
  retx = dens$x; retx[retx > demoBL] = demoBL
  sumret = sum(dens$y*retx)
  ret = sumret/sumf

  data.frame(BL = CRvals, PgBL = prob, ret = ret) # return tentative bag limits and probs of being over
}

calc_fit = function(dens,size,mu, demoBL = 54){

  # probability trip over
  ix = min(floor(dens$x)):max(floor(dens$x))       # integer trip limit
  prob = pnbinom(ix,size=size,mu=mu,lower.tail = F)
  prob = prob / max(prob)

  # probability fish over
  densfit = dnbinom(ix,size=size,mu=mu)
  sumf=sum(densfit*ix)
  retix = ix; retix[retix>demoBL] = demoBL
  sumret = sum(densfit*retix)
  ret = sumret/sumf

  data.frame(BL = ix, PgBL = prob, ret=ret) # return tentative bag limits and probs of being over
}

fit_mod_20_24_S = function(dat, cutoff = 20, maxplot= 70,demoBL = 54){
  dat1 = dat[dat$Area %in% c("FLORIDA","GEORGIA","SOUTH CAROLINA","NORTH CAROLINA") & dat$Year%in%(2020:2024), ]
  Area = c(1,3)[as.numeric(dat1$Area == "FLORIDA")+1] # area assignment
  Fleet = c(5,3)[as.numeric(dat1$Mode =="PRIVATE/RENTAL")+1] # fleet assignment
  Fleetnams = rep(NA,8); Fleetnams[3] = "RecS"; Fleetnams[5] = "HireS"
  Areanams = rep(NA,5); Areanams[1] = "CAR+FLK"; Areanams[3] = "NCFL"
  Quarter = floor((as.numeric(dat1$Month)/1.01)/3)+1
  CRmu = aggregate(dat1$Catch,by=list(Q = Quarter, Ft = Fleet, A = Area),FUN=mean)
  CRsd = aggregate(dat1$Catch,by=list(Q = Quarter, Ft = Fleet, A = Area),FUN=sd)
  CRn =  aggregate(rep(1,nrow(dat1)),by=list(Q = Quarter, Ft = Fleet, A = Area),FUN=sum)
  keep = CRn$x >= cutoff
  CRmu = CRmu[keep,]
  CRsd = CRsd[keep,]
  nfit = sum(keep)
  fitmu = fitcv = fitsize = fitnll = rep(NA,nfit)
  breaks = seq(0,maxplot,2)
  par(mfrow=c(ceiling(nfit/3),3),mai = c(0.3,0.3,0.45,0.01), omi=c(0.4,0.4,0.01,0.01))
  emp_prob = fit_prob = emp_dens = list()
  for(i in 1:nfit){
    qq = CRmu$Q[i]; ff = CRmu$Ft[i]; aa = CRmu$A[i]
    catchy = dat1$Catch[Quarter == qq & Fleet == ff & Area == aa]
    cplot=catchy[catchy<=maxplot]
    dens = density(cplot,from=1,to=maxplot,adjust=2)
    emp_dens[[i]] = dens
    emp_prob[[i]] = calc_emp(dens,demoBL=demoBL)
    plot(dens, ylim=c(0,max(dens$y)),main=paste0("Q = ", qq, ", ",Fleetnams[ff], ", ",Areanams[aa])); grid()
    mu =fitmu[i] = CRmu$x[i]
    opt = optimize(nbinfit, c(-3,5), lnmu = log(CRmu$x[i]), catchy=catchy)
    size = exp(opt$minimum)
    nbplot(size,mu=mu,catchy,maxplot)
    fitnll[i] = opt$objective
    r = fitsize[i] = size
    p = r/(r+mu)
    var = r*(1-p)/(p*p)
    sd = var^0.5
    cv = fitcv[i] = sd/mu
    fit_prob[[i]] = calc_fit(dens, size, mu, demoBL = demoBL)
    legend('topright',legend=c(paste0("Mean = ", round(CRmu$x[i],2), ", sd = ",round(CRsd$x[i],2), ", CV = ",round(CRsd$x[i]/CRmu$x[i],2)),
                               paste0("Fit: Size(r) = ", round(size,2),", sd = ", round(sd,2), ", CV = ",round(cv,2))),
                            bty="n",text.col=c("black","red"))
    legend('right',legend=paste0("n = ",length(catchy)),text.col="blue",bty='n')
    mtext(paste0("(",letters[i],")"),adj=0.02,line=0.1,cex=0.75)
    abline(v=demoBL ,col="darkgreen", lty=2)

  }
  mtext("Relative Frequency",2,line = 0.8,outer=T)
  mtext("Catch per trip (MRIP 2020-2024)",1, outer=T)

  list(dat1 = dat1, CRmu = CRmu, CRsd = CRsd, CRn = CRn, fitmu = fitmu, fitcv = fitcv, fitsize = fitsize, emp_dens=emp_dens, emp_prob = emp_prob, fit_prob = fit_prob)

}

probplot = function(fits, demoBL = 54){
  nfit = nrow(fits$CRmu)
  respind = 2

  par(mfrow=c(ceiling(nfit/3),3),mai = c(0.3,0.3,0.45,0.01), omi=c(0.4,0.4,0.01,0.01))
  for(i in 1:nfit){
    ep = fits$emp_prob[[i]][,c(1,respind)]
    fp = fits$fit_prob[[i]][,c(1,respind)]
    eret =  fits$emp_prob[[i]][1,3]
    fret = fits$fit_prob[[i]][1,3]
    qq = fits$CRmu$Q[i]; ff = fits$CRmu$Ft[i]; aa = fits$CRmu$A[i]
    plot(ep, main=paste0("Q = ", qq, ", ",Fleets[ff], ", ",Areas[aa]),type="l"); grid()
    lines(fp,col="red")
    mtext(paste0("(",letters[i],")"),adj=0.02,line=0.1,cex=0.75)
    if(i==1)legend('bottomright',legend=c("Empirical","Neg. Bin model"),text.col=c("black","red"),bty="n")
    abline(v=demoBL ,col="darkgreen", lty=2)
    emp_pBL = ep[ep[,1]==demoBL,2]
    fit_pBL = fp[fp[,1]==demoBL,2]
    legend('topright',legend=c(paste0("P(trip CR)>BL = ",round(emp_pBL*100,2),"%"),
                               paste0("P(trip CR)>BL = ",round(fit_pBL*100,2),"%"),
                               paste0("Frac. Fish Rel.>BL = ",round((1-eret)*100,2),"%"),
                               paste0("Frac. Fish Rel.>BL = ",round((1-fret)*100,2),"%")),text.col=c("black","red"),bty="n")


  }
  mtext("Prob. a trip exceeds catch rate",2,line=0.8,outer=T)
  #if(plotfish)mtext("Expected fraction of fish released if bag limit set at catch per trip",2,line=0.8,outer=T)
  mtext("Catch per trip (MRIP 2020-2024)",1, outer=T)
}


