

mahiTplot = function(MM,
                     Measures =   c("Brel","Frel","USCom Landings","RecS Landings","HireS Catch Rate","USCom Catch Dif","USCom Mean length","HireS Mean length"),
                     TimePeriod = c("all", "all", "all",            "all",          "all",             "all",             "all",              "all"),
                     Lab =        c("SSB/SSBMSY", "F/FMSY", "US Com. Landings (t)", "Rec S. Landings (t)",
                                    "Hire S. Catch Rate (rel. 2022)", "US Com. Av, Ann. Var, Yield", "US Com. Mean Length (rel. 2022)",
                                    "Hire S. Mean Length (rel. 2022)"),
                     cols = c("cornflowerblue","orange","red2","forestgreen","black","darkblue","tan4","indianred1","darkgrey","magenta"),
                     leg=T
                     ){


    nm = length(Measures)
    np = nm/2
    ncol = ceiling(np^0.5)
    nrow = ceiling(np / ncol)
    par(mfrow=c(nrow,ncol),mai=c(0.8,0.8,0.1,0.1))
    Tcols = MSEtool::makeTransparent(cols,90)
    mind = 1:(nm/2)*2

    for(mm in mind){
      MPnames = rownames(MM[[Measures[mm-1]]][[TimePeriod[mm-1]]])
      denom=1
      if(grepl("Landings",Measures[mm-1]))denom=1000
      xmu = MM[[Measures[mm-1]]][[TimePeriod[mm-1]]][,"50%"]/denom
      xLB = MM[[Measures[mm-1]]][[TimePeriod[mm-1]]][,"5%"]/denom
      xUB = MM[[Measures[mm-1]]][[TimePeriod[mm-1]]][,"95%"]/denom
      denom=1
      if(grepl("Landings",Measures[mm]))denom=1000
      ymu = MM[[Measures[mm]]][[TimePeriod[mm-1]]][,"50%"]/denom
      yLB = MM[[Measures[mm]]][[TimePeriod[mm-1]]][,"5%"]/denom
      yUB = MM[[Measures[mm]]][[TimePeriod[mm-1]]][,"95%"]/denom

      plot(range(xLB,xUB),range(yLB,yUB),col="white",xlab="",ylab=""); grid()
      points(xmu,ymu,pch=19,cex=1.1,col=Tcols)
      nmp = length(xmu)
      for(mp in 1:nmp){
        lines(c(xLB[mp],xUB[mp]),rep(ymu[mp],2),col=Tcols[mp])
        lines(rep(xmu[mp],2),c(yLB[mp],yUB[mp]),col=Tcols[mp])
      }
      if(mm==2 & leg)legend('topright',legend=MPnames,text.col=cols,bty='n')
      mtext(Lab[mm-1],1,line=2.2)
      mtext(Lab[mm],2,line=2.2)
    }

}
