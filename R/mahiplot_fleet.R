
# ftno =1; varnam = "Landings_byfleet"; dorange=T; varlab = "Landings (t)"; MPs = 4:5; denom = 1E3; stock= 1; fcols = c("#10109950","#99000050"); ns = 4; annual = T
MP_comp_fleet = function(MR, ftno, varnam = "Landings_byfleet", varlab = "Landings (t)", MPs = 1:2, denom = 1E3, stock= 1, fcols = c("#0000ff50","#ff000050","#00ff0050"),ns, annual = T,dorange=T, addn=T, addmed=T, addl = T){

  dms = MR$dms
  np = dms$np; nt = dms$nt; lhy = MR$CurrentYear; py = np/ns; nsim = dms$nsim; nmp = dms$nMP; nom=dms$nMSE
  var0 = MR[[varnam]] /denom
  tsind = nt/ns + 1:(np/ns)
  npy = length(tsind)
  var1 = var0[,,tsind,ftno,]
  var = array(var1,c(nom*nsim,npy,nmp))
  lab = lhy + 1:py

  ylim = c(0,quantile(var,0.9975))
  posMPs = 1:nmp
  MPs = MPs[MPs%in%posMPs]
  for(mp in MPs){
    mat = var[,,mp]
    proj_plot(mat,lab, ylab = "",add=mp!=MPs[1],fcol=fcols[match(mp,MPs)],
              ylim = ylim, dorange = dorange, addn = addn, addmed = addmed, addl = addl)
    #abline(v=MSEobj@OM@CurrentYear+0.5,lty=2)
  }
  for(mp in MPs){
    mat = var[,,mp]
    matplot(lab,apply(mat,2,median),type="l",lwd=2,lwy=1,add=T,col=fcols[match(mp,MPs)])
    matplot(lab,apply(mat,2,median),type="l",lwd=2,lwy=1,add=T,col=fcols[match(mp,MPs)])
  }

  mtext(Fleets[ftno],line=0.4,cex=1,font=2)
}

#' Plot overall yield for various fleets across MPs
#'
#' @param MSEobj A class of object 'MSE'
#' @param MPs An integer vector of MPs to plot
#' @param fts An interger veector of Fleets to plot (see object Fleets)
#' @param annual Boolean, should annual or seasonal results be plotted?
#' @param ns Integer, number of seasons.
#' @param MPcols Character vector of colors for plotting MP projections
#' @param npy Positive integer, number of plot years - No. historical years to plot before projection
#' @param MPnam Character vector, optional, a character vector of MP names
#' @param leg Boolean, add a legend for the MPs?
#' @param dorange Boolean, add the quantile ranges?
#' @param addn Boolean, add narrow interquartile range shading?
#' @param addmed Boolean, add a line for the median?
#' @param addl Boolean, add individual simulation lines?
#' @examples
#' mahiplot_fleet(myMSE)
#' @author T. Carruthers
#' @export
mahiplot_fleet = function(MR, MPs = 1:3, fts = 1:5, annual = T, ns = 4, stock = 1,  MPcols = c("#0000ff50","#ff000050","#00ff0050"),
                    npy=20, MPnams = NA,leg=T, dorange = T, addn=F, addmed=T, addl = F){
  # MR =mahiResults(myMSE); MPs = 1:2; fts = 1:5; dorange=T; annual = T; ns = 4; stock = 1;  MPcols = c("red","green","blue","orange","grey","purple"); npy=20; MPnams = NA; leg=F; addn=F; addmed=T; addl = F
  nf = length(fts); fnams = Fleets[fts]
  ncol = ceiling(nf^0.5); nrow = ceiling(nf/ncol)
  par(mfrow=c(nrow,ncol),mai=c(0.25,0.25,0.35,0.05),omi=c(0.5,0.5,0.01,0.01))
  for(ff in 1:nf)  MP_comp_fleet(MR, fts[ff],ns=ns,stock=stock, MPs=MPs, dorange=dorange, addn = addn, addmed = addmed, addl = addl, annual = annual,fcols=MPcols)

  if(leg){
    if(is.na(MPnams[1]))MPnams = names(MR$MPs)[MPs]
    if(ncol*nrow>nf){
      plot(x=0,y=0,axes=F,col="white",xlab="",ylab="")
      legend('center',legend = MPnams,text.col=MPcols,bty='n',text.font=2,cex=1.2)
      legend('center',legend = MPnams,text.col=MPcols,bty='n',text.font=2,cex=1.2)
    }else{
      legend('bottomright',legend = MPnams,text.col=MPcols,bty='n',text.font=2)
      legend('bottomright',legend = MPnams,text.col=MPcols,bty='n',text.font=2)
    }
  }
  mtext("Projection Year",1,outer=T,line=1.2)
  mtext("Landings (t)",2,outer=T,line=1.2)

}
