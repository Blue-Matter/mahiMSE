# bag limit analysis and code


process_CR_RR = function(dir, Targ=F){ # process catch rate vs release rate

  d4d = read_xlsx(dir)

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

nbplot = function(size,mu,catchy,col="red"){
  k = 1:50

  dd = dnbinom(k,size=size,mu=mu)
  lines(k,dd,col=col)
}


fit_mod_20_24_S = function(dat, cutoff = 20){
  dat1 = dat[dat$Area %in% c("FLORIDA","GEORGIA","SOUTH CAROLINA","NORTH CAROLINA") & dat$Year%in%(2020:2024), ]
  Area = (4:5)[as.numeric(dat1$Area == "FLORIDA")+1] # area assignment
  Fleet = c(5,3)[as.numeric(dat1$Mode =="PRIVATE/RENTAL")+1] # fleet assignment
  Fleetnams = rep(NA,8); Fleetnams[3] = "RecS"; Fleetnams[5] = "HireS"
  Areanams = rep(NA,7); Areanams[4] = "FLK"; Areanams[5] = "NCFL"
  Quarter = floor((as.numeric(dat1$Month)/1.01)/3)+1
  CRmu = aggregate(dat1$Catch,by=list(Q = Quarter, Ft = Fleet, A = Area),FUN=mean)
  CRsd = aggregate(dat1$Catch,by=list(Q = Quarter, Ft = Fleet, A = Area),FUN=sd)
  CRn =  aggregate(rep(1,nrow(dat1)),by=list(Q = Quarter, Ft = Fleet, A = Area),FUN=sum)
  keep = CRn$x >= cutoff
  CRmu = CRmu[keep,]
  CRsd = CRsd[keep,]
  nfit = sum(keep)
  fitmu = fitcv = fitsize = fitnll = rep(NA,nfit)
  maxplot = 50
  breaks = seq(0,maxplot,2)
  par(mfrow=c(ceiling(nfit/4),4),mai = c(0.3,0.3,0.3,0.01), omi=c(0.4,0.4,0.01,0.01))
  for(i in 1:nfit){
    qq = CRmu$Q[i]; ff = CRmu$Ft[i]; aa = CRmu$A[i]
    catchy = dat1$Catch[Quarter == qq & Fleet == ff & Area == aa]
    cplot=catchy[catchy<=maxplot]
    plot(density(cplot,from=1,to=50,adjust=2), main=paste0("Quarter = ", qq, ", ",Fleetnams[ff], ", ",Areanams[aa]))
         #main = paste0("Q = ",qq, ", F = ", ff, ", A = ",aa))#,breaks=breaks)

    mu =fitmu[i] = CRmu$x[i]
    opt = optimize(nbinfit, c(-3,5), lnmu = log(CRmu$x[i]), catchy=catchy)
    size = exp(opt$minimum)
    nbplot(size,mu=mu,catchy)
    #opt2 = optim(c(opt$minimum,log(CRmu$x[i])),nbinfit2,catchy=catchy)
    fitnll[i] = opt$objective
    #nbplot(exp(opt2$par[1]),mu=exp(opt2$par[2]),catchy,col="blue")
    r = fitsize[i] = size #exp(opt2$par[1])
    p = r/(r+mu)
    var = r*(1-p)/(p*p)
    sd = var^0.5
    cv = fitcv[i] = sd/mu

    legend('topright',legend=c(paste0("Mean = ", round(CRmu$x[i],2), ", sd = ",round(CRsd$x[i],2), ", CV = ",round(CRsd$x[i]/CRmu$x[i],2)),
                               paste0("Fit: Size(r) = ", round(size,2),", sd = ", round(sd,2), ", CV = ",round(cv,2))),
                            bty="n",text.col=c("black","red"))


  }
  mtext("Relative Frequency",2,outer=T)
  mtext("Catch per trip (MRIP 2020-2024)",1, outer=T)

  list(dat1 = dat1, CRmu = CRmu, CRsd = CRsd, CRn = CRn, fitmu = fitmu, fitcv = fitcv, fitsize = fitsize)

}
