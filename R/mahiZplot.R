
mahiZplot = function(MM,
                     Measures =   c("Brel","Frel","USCom Landings","RecN Landings","RecS Landings","HireN Landings","HireS Landings"),
                     TimePeriod = c("all", "all", "all",            "all",          "all",             "all",             "all",              "all"),
                     Lab =        c("SSB/SSBMSY", "F/FMSY", "US Com. Landings (t)", "Rec N. Landings (t)","Rec S. Landings (t)","Nire N. Landings (t)","Hire S. Landings (t)"),
                     cols = c("cornflowerblue","orange","red2","forestgreen","black","darkblue","tan4","indianred1","darkgrey","magenta"),
                     leg=T
){


  nm = length(Measures)
  ncol = ceiling(nm^0.5)
  nrow = ceiling(nm / ncol)
  par(mfrow=c(nrow,ncol),mai=c(0.8,0.6,0.1,0.1),omi=c(0,0,0,0))
  Tcols = MSEtool::makeTransparent(cols,90)

  for(mm in 1:nm){
    MPnames = rownames(MM[[Measures[mm]]][[TimePeriod[mm]]])
    denom=1
    if(grepl("Landings",Measures[mm]))denom=1000
    y5 = MM[[Measures[mm]]][[TimePeriod[mm]]][,"5%"]/denom
    y25 = MM[[Measures[mm]]][[TimePeriod[mm]]][,"25%"]/denom
    y50 = MM[[Measures[mm]]][[TimePeriod[mm]]][,"50%"]/denom
    y75 = MM[[Measures[mm]]][[TimePeriod[mm]]][,"75%"]/denom
    y95 = MM[[Measures[mm]]][[TimePeriod[mm]]][,"95%"]/denom

    nmp = length(y5)
    plot(c(0.5,nmp+0.5),range(y95,y5,na.rm=T),col="white",xlab="",ylab="",axes=F); grid()
    axis(2)
    axis(1,c(-1000,1000),rep("",2))
    axis(1,at=1:nmp, MPnames, las = 2)
    points(1:nmp,y50,pch=19,cex=1.1,col=Tcols)
    for(mp in 1:nmp){
      lines(rep(mp,2),c(y5[mp],y95[mp]),col=Tcols[mp],lwd=2)
      lines(rep(mp,2),c(y25[mp],y75[mp]),col=Tcols[mp],lwd=4)
    }
    mtext(Lab[mm],2,line=2.2,cex=0.8)

  }


}
