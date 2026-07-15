# mahiMSE plotting functions
# MSE = miniMSE

# MSE = Project(Hist, c("mahiMP","mahiMP2"))

# sims = 1:2; qsw = c(0.05,0.95); qsn = c(0.25,0.75); y0 = TRUE; ylim = NA; addn=T; addmed=T; addl = T; xlab = "Projection Year"; ylab = "Value"; fcol= "#10109950"
proj_plot = function(mat, lab, sims = 1:2, qsw = c(0.05,0.95), qsn = c(0.25,0.75), y0 = TRUE, ylim = NA, dorange = T, addn=T, addmed=T, addl = T,
                     xlab = "Projection Year", ylab = "Value",fcol= "#10109950", add = F){
  yw = apply(mat,2,quantile,p=qsw)
  yn = apply(mat,2,quantile,p=qsn)
  med = apply(mat,2,median)
  if(is.na(ylim[1])){
    ylim = range(yw)
    if(y0) ylim[1] = 0
  }
  if(!add)matplot(range(lab),ylim,col="white",xlab = xlab, ylab = ylab)
  grid()
  if(dorange)polygon(c(lab,rev(lab)),c(yw[1,],rev(yw[2,])),col=fcol, border=fcol)
  if(addn) polygon(c(lab,rev(lab)),c(yn[1,],rev(yn[2,])),col=fcol, border=fcol)
  if(addmed) lines(lab,med,col="white",lwd=2)
  if(addl) matplot(lab,t(mat[sims,]),col="black",lty=1:2,lwd=1,type="l",add=T)

}

# varnam = "SBiomass"; varlab = "Spawning Stock Biomass (t)"; MPs = 1:2; denom = 1E3; stock= 1; fcols = c("#10109950","#99000050"); ns = 4; annual = T
MP_comp = function(MR, varnam = "SSB", varlab = "Spawning Stock Biomass (t)", MPs = 1:2, denom = 1E3,  fcols = c("#0000ff50","#ff000050","#00ff0050"),ns, annual = T, dorange = T, addn=T, addmed=T, addl = T, firstyr = 2005, lastyr = 2042){

  dms = MR$dms
  np = dms$np; nt = dms$nt; lhy = MR$CurrentYear; py = np/ns; nsim = dms$nsim; nmp = dms$nMP; nOM = dms$nMSE; ny=dms$ny
  var0 = MR[[varnam]]/denom

  var = array(var0,c(nOM*nsim,ny,nmp))

  if(annual){
    lab =  unique(floor(Time_steps))
    newvar = array(var, c(nsim, ns, ny/ns, nmp))
    var = apply(newvar,c(1,3,4),mean)
   }else{
    lab = Time_steps #rep(lhy + 0:(py-1),each=ns) + seq(1/ns, 1, length.out = ns)
   }

  ylim = c(0,quantile(var,0.9975))
  #par(mai = c(0.2,0.7,0.15,0.15),omi=)
  posMPs = 1:nmp
  MPs = MPs[MPs%in%posMPs]
  keep = !(lab < firstyr | lab > lastyr)

  for(mp in MPs){
    mat = var[,,mp]
    proj_plot(mat[,keep],lab[keep], ylab = varlab,add=mp!=MPs[1],fcol=fcols[match(mp,MPs)],
              ylim = ylim, dorange=dorange, addn = addn, addmed = addmed, addl = addl)
  }
  for(mp in MPs){
    mat = var[,,mp]
    matplot(lab[keep],apply(mat[,keep],2,median),type="l",lwd=2,lwy=1,add=T,col=fcols[match(mp,MPs)])
    matplot(lab[keep],apply(mat[,keep],2,median),type="l",lwd=2,lwy=1,add=T,col=fcols[match(mp,MPs)])
  }


}

#' Plot overall yield and biomass projections across MPs
#'
#' @param MR A class of object 'mahiResults' made by the function mahiResults(MSElist)
#' @param MPs An integer vector of MPs to plot
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
#' mahiplot(mahiResults(myMSE))
#' @author T. Carruthers
#' @export
mahiplot = function(MR, MPs = 1:2, annual = T, ns = 4, stock = 1,  MPcols = c("#0000ff50","#ff000050","#00ff0050"),
                    npy=20, MPnams = NA,leg=T, dorange = T, addn=F, addmed=T, addl = F, firstyr = 2005, lastyr = 2042){
  # MR =mahiResults(myMSE); MPs = 1:2; annual = T; ns = 4; stock = 1;  MPcols = c("red","green","blue","orange","grey","purple"); npy=20; MPnams = NA; leg=F; addn=F; addmed=T; addl = F
  par(mfrow=c(1,2),mai=c(0.25,0.85,0.15,0.05),omi=c(0.5,0.01,0.01,0.01))
  MP_comp(MR, ns=ns, MPs = MPs, dorange=dorange, addn = addn, addmed = addmed, addl = addl, annual=annual, fcols = MPcols, firstyr = firstyr, lastyr = lastyr)
  abline(v=MR$CurrentYear+0.5,lty=2)
  MP_comp(MR, varnam = "Landings", varlab = "Landings (t)",ns=ns, MPs=MPs, dorange=dorange, addn = addn, addmed = addmed, addl = addl, annual = annual,fcols=MPcols, firstyr = firstyr, lastyr = lastyr)
  abline(v=MR$CurrentYear+0.5,lty=2)
  if(leg){
    legend('bottomright',legend = names(MR$MPs)[MPs],text.col=MPcols,bty='n',text.font=2)
    legend('bottomright',legend = names(MR$MPs)[MPs],text.col=MPcols,bty='n',text.font=2)
  }
  mtext("Projection Year",1,outer=T,line=1.2)
}


#' @export
mahiexplanatory = function(MSEobj, ny=2, ny_proj=12, stock = 1, MP = 1, sims = 1:2, ns = 4, simcol = MPcols,doylab=T,doxlab=T, ws=1.15){
  # ny = 2; ny_proj=15; stock = 1; MP = 1; sims = 1:2; ns = 4; simcol = c("#ff0000","#0000ff"); ws=1.05

  pt = MSEobj@OM@pYear; nt = MSEobj@OM@nYear; lhy = MSEobj@OM@Data$Dolphinfish@YearLH; py = pt/ns; nsim = nSim(MSEobj); nmp = length(MSEobj@MPs)

  labs = lhy +1:py

  B1 = MSEobj@Biomass[sims, stock, ,MP, drop=T]
  Bagg = apply(array(B1,c(length(sims),4,py)),c(1,3),mean)/1E6
  B = Bagg[,1:ny]
  Bylim = c(0,max(Bagg[,1:ny_proj])*ws)

  L1 = apply(MSEobj@Landings[sims, stock, ,,MP, drop=T],c(1,2),sum)
  Lagg = apply(array(L1,c(length(sims),4,py)),c(1,3),sum)/1E3
  L = Lagg[,1:(ny+1)]
  Lylim = c(0,max(Lagg[,1:ny_proj])*ws)

  ylab = "Biomass Index"; if(!doylab)ylab=""
  xlab = "Projection Year"; if(!doxlab)xlab=""
  matplot(labs[1:ny_proj], t(Bagg[,1:ny_proj]), type="l", col=simcol, lty=2, ylab=ylab, xlab = xlab, ylim=Bylim); grid()
  matplot(labs[1:ny], t(B), type="l", col=simcol, lty = 1, lwd=2 , add=T)
  points(rep(labs[ny],2),B[,ny],pch=19,cex=1.4,col=simcol)
  legend('top', legend=paste0("Biomass Index in ",labs[ny]),bty="n")

  ylab = "TAC (t)"; if(!doylab)ylab=""
  xlab = "Projection Year"; if(!doxlab)xlab=""
  matplot(labs[1:ny_proj], t(Lagg[, 1:ny_proj]), type="l", col=simcol, lty=2, ylab=ylab, xlab = xlab,ylim=Lylim); grid()
  matplot(labs[1:(ny+1)], t(L), type="l", col=simcol, lty = 1, lwd=2 , add=T)
  points(rep(labs[ny+1],2),L[,ny+1],pch=19,cex=1.4,col=simcol)
  legend('top', legend=paste0("Response in TAC in ",labs[ny+1]),bty="n")

}

