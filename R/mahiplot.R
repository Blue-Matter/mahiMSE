# mahiMSE plotting functions
# MSE = miniMSE


getChist = function(MSE){
  mHist = MSE@multiHist[[1]]
  Cr = lapply(mHist,function(x)apply(x@TSdata$Landings,1:2,sum))
  nr = length(Cr);  nsim = nrow(Cr[[1]]); nt=ncol(Cr[[1]])
  array(unlist(Cr),c(nsim,nt,nr))
}


mahiqplot = function(mu, sims, xval, vline, ny, MPcols, ylab, npy, iqr = 0.70){
  ay1 = dim(mu)[2];   keepy = (ny-npy+1):ay1;
  mu1=mu[,keepy]; sims1 = sims[,,keepy]; xval1 = xval[keepy]; ny1=npy
  nmp = dim(mu1)[1]; ay = dim(mu1)[2]
  qs = apply(sims1,2:3,quantile,p=c(((1-iqr)/2),(1-(1-iqr)/2)))
  matplot(xval1,t(mu1),col="white",type="l",ylim=c(0,max(qs)),yaxs="i", xlab="", ylab=""); grid()
  polygon(xval1[c(1:ny1,ny1:1)],c(qs[1,1,1:ny1],qs[2,1,ny1:1]),col="darkgrey",border=NA)
  pys = ny1:ay
  for(mp in 1:nmp)polygon(xval1[c(pys,rev(pys))],c(qs[1,mp,pys],qs[2,mp,rev(pys)]),col=makeTransparent(MPcols[mp],70),border=NA)
  matplot(xval1,t(mu1),col=MPcols,type="l",add=T, lty=1,lwd=2)
  lines(xval1[1:ny1], mu1[1,1:ny1],col="black",lwd=2)
  abline(v=vline,lty=2)
  mtext(ylab,2,line=2.2,cex=0.95)
}


#' Plot overall yield and biomass projections across MPs
#'
#' @param MSE A class of object 'MMSE' multiMSE
#' @param iqr Positive fraction, the interquartile range - default is 0.8 - intervals shown are the 10 percent and 90 percent quantiles of simulations
#' @param MPcols Character vector of colors for plotting MP projections
#' @param npy Positive integer, number of plot years - No. historical years to plot before projection
#' @examples
#' mahiplot(myMSE)
#' @author T. Carruthers
#' @export
mahiplot = function(MSE, iqr = 0.8, MPcols = c("red","green","blue","orange","grey","purple"),npy=20){

  nsubyr = 4
  MPnams = unlist(MSE@MPs)
  CurrentYr = MSE@OM[[1]][[1]]$CurrentYr[1]
  nsim = MSE@nsim; nMP = MSE@nMPs
  nt = MSE@nyears;  pt = MSE@proyears;   at = nt+pt;   ny = nt/nsubyr; py = pt/nsubyr; ay = at/nsubyr
  yind = rep(1:ny, each=4)
  ylab = CurrentYr + (-ny +1):py

  # Calcs
  histy = MSE@multiHist[[1]][[1]]

  # SSB calcs
  SSB = MSE@SSB[,1,,] /1E6 # kt only one stock SMT Simulation, MP, Timestep
  hSSB = apply(histy@TSdata$SBiomass,1:2,sum) / 1E6# ST , Simulation Time step
  ASSB = array(NA, c(nsim, nMP, at)) # SMT
  for(MP in 1:nMP)ASSB[,MP,1:nt] = hSSB # historical SSB by MP fill
  ASSB[,,nt+1:pt] = SSB              # projected SSB
  ASSBy = apply(array(ASSB, c(nsim, nMP, nsubyr, ay)),c(1,2,4),mean)
  muASSBy = apply(ASSBy,2:3,mean)

  # Yield calcs
  C = apply(MSE@Catch, c(1,4,5), sum) /1E6 # SMT
  hCr = getChist(MSE) # STR
  hC = apply(hCr,1:2,sum) / 1E6  # ST
  AC = array(NA, c(nsim, nMP, at))  # SMT
  for(MP in 1:nMP)AC[,MP,1:nt] = hC # historical C by MP fill
  AC[,,nt+1:pt] = C                 # projected C
  ACy = apply(array(AC, c(nsim, nMP, nsubyr, ay)),c(1,2,4),sum)
  muACy = apply(ACy,2:3,mean)

  par(mfrow=c(1,2),mai=c(0.25,0.65,0.03,0.03),omi=c(0.5,0.01,0.01,0.01))
  mahiqplot(muACy, ACy, ylab, CurrentYr+0.5, ny, MPcols,"Yield (kt)",npy,iqr=iqr)
  legend('topright',legend=MPnams,text.col=MPcols,text.font=2,cex=0.85,bty='n')
  mahiqplot(mu=muASSBy,sims=ASSBy,xval = ylab, vline = CurrentYr+0.5, ny=ny, MPcols = MPcols, "SSB (kt)",npy, iqr=iqr)
  mtext("Year",1,line=1,outer=T)
}
