#' Plots projected exploitation rate and spawning biomass relative to MSY across MPs
#'
#' @param MR A mahiresults object produced by mahiResults(MSElist)
#' @param MPs An integer vector of MPs to plot
#' @param MPcols Character vector of colors for plotting MP projections
#' @param MPnam Character vector, optional, a character vector of MP names
#' @param leg Where to place the legend. NA means no legend will be added.
#' @param dorange Boolean, add quantile ranges?
#' @param addn Boolean, add narrow interquartile range shading?
#' @param addmed Boolean, add a line for the median?
#' @param addl Boolean, add individual simulation lines?
#' @param firstyr Calendar year, first year to plot results
#' @param lastyr Calendar year, last year to plot results
#' @examples
#' mahiplot_MSY(mahiResults(myMSE))
#' @author T. Carruthers
#' @export
mahiplot_MSY = function(MR, MPs = 1:2, MPcols = c("#0000ff50","#ff000050","#00ff0050"),
                     MPnams = NA,leg='topleft', dorange = T, addn=F, addmed=T, addl = F, firstyr = 2005, lastyr = 2042){
  # MR = mahiResults(list(myMSE, myMSE)); MPs = 1:2; MPcols = c("#0000ff50","#ff000050","#00ff0050"); MPnams = NA; leg=T; addn=F; addmed=T; addl = F; firstyr = 2005; lastyr = 2042

  ns = 4
  par(mfrow=c(1,2),mai=c(0.25,0.85,0.15,0.05),omi=c(0.5,0.01,0.01,0.01))
  ysteps = unique(floor(Time_steps))
  MP_comp2(MR, varnam = "Brel", varlab ="SSB / SSBMSY", xlab = ysteps, ns=ns, MPs = MPs, dorange=dorange, addn = addn, addmed = addmed, addl = addl, annual=annual, fcols = MPcols, firstyr = firstyr, lastyr=lastyr)
  abline(h=1,lty=2)
  abline(v=MR$CurrentYear+0.5,lty=2)
  MP_comp2(MR, varnam = "Frel", varlab = "F / FMSY", xlab = ysteps,ns=ns, MPs=MPs, dorange=dorange, addn = addn, addmed = addmed, addl = addl, annual = annual,fcols=MPcols, firstyr = firstyr, lastyr=lastyr)
  abline(h=1,lty=2)
  abline(v=MR$CurrentYear+0.5,lty=2)
  if(!is.na(leg)){
    legend(leg,legend = names(MR$MPs)[MPs],text.col=MPcols,bty='n',text.font=2)
    legend(leg,legend = names(MR$MPs)[MPs],text.col=MPcols,bty='n',text.font=2)
  }
  mtext("Projection Year",1,outer=T,line=1.2)
}


# varnam = "Brel"; varlab = "SSB / SSBMSY"; xlab = Time_steps; MPs = 1:2; denom = 1; fcols = c("#10109950","#99000050"); ns = 4; annual = T; addn=T; addmed=T; addl = T; firstyr = 2000; lastyr = 2042
MP_comp2 = function(MR, varnam = "Brel", varlab = "SSB / SSBMSY", xlab, MPs = 1:2, denom = 1, fcols = c("#0000ff50","#ff000050","#00ff0050"),ns, annual = T,dorange = T, addn=T, addmed=T, addl = T, firstyr = 2000, lastyr = 2042){

  dms = MR$dms
  nmp =  dms$nMP; nom = dms$nMSE; nsim = dms$nsim
  var0 = MR[[varnam]]
  ny = dim(var0)[3]

  var = array(var0,c(nsim*nom,ny,nmp))

  ylim = c(0,quantile(var,0.9975))
  posMPs = 1:nmp
  MPs = MPs[MPs%in%posMPs]
  keep = !(xlab < firstyr | xlab > lastyr)

  i = 0
  for(mp in MPs){
    i=i+1
    mat = var[,,mp]
    proj_plot(mat[,keep], lab = xlab[keep], ylab = varlab, add=mp!=MPs[1], fcol=fcols[i],
              ylim = ylim, dorange = dorange, addn = addn, addmed = addmed, addl = addl)
  }

  i=0
  for(mp in MPs){
    i=i+1
    mat = var[,,mp]
    matplot(xlab[keep],apply(mat[,keep],2,median),type="l",lwd=2,lwy=1,add=T,col=fcols[i])
    matplot(xlab[keep],apply(mat[,keep],2,median),type="l",lwd=2,lwy=1,add=T,col=fcols[i])
  }


}
