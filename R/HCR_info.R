HCR_info = function(anMP, HCRno = 1, Iref = 1.5, leg=T, xmax=NA, ymax=NA, domain=T){

  if(is.na(formals(anMP)$HCR_ICP)) stop("No HCR specified in this MP")

  Ilab = formals(anMP)$Index_names[formals(anMP)$HCR_input]
  ICP = formals(anMP)$HCR_ICP[[HCRno]]
  LCP = formals(anMP)$HCR_LCP[[HCRno]]
  Lev = formals(anMP)$HCR_lever[HCRno]

  xlim = c(0,max(2,   ICP*1.25, xmax, na.rm=T))
  ylim = c(0,max(1.5, LCP*1.25, ymax, na.rm=T))

  plot(xlim,ylim,col='white',xlab="",ylab=""); grid()
  abline(v=ICP,col="blue",lty=2)
  abline(h=LCP,col="green",lty=2)
  lines(c(0,ICP[1]),rep(LCP[1],2),lwd=4)
  lines(ICP,LCP,lwd=4)
  lines(c(ICP[2],1000),rep(LCP[2],2),lwd=4)

  abline(v=Iref,col="red")
  abline(h=1,col="red")

  mtext(paste0(Ilab," Index"),1,line=2.2)
  mtext(paste0("Fraction of 2022 ",Lev),2,line=2.2)
  if(leg)legend('topleft',c("Index control points",paste0(Lev, " control points"), "2022 Levels", "HCR"),text.col=c("blue","green","red","black"))
  if(domain)mtext(deparse(substitute(anMP)),line=0.5,font=2,cex=1.2)

}
