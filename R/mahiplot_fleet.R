

get_cat_area = function(MSE){
  nsubyr = 4
  MPnams = unlist(MSE@MPs)
  CurrentYr = MSE@OM[[1]][[1]]$CurrentYr[1]
  nsim = MSE@nsim; nMP = MSE@nMPs
  nt = MSE@nyears;  pt = MSE@proyears;   at = nt+pt;   ny = nt/nsubyr; py = pt/nsubyr; ay = at/nsubyr
  yind = rep(1:ny, each=4)
  ylab = CurrentYr + (-ny +1):py

  ED = EffDist$C_qfr
  C = MSE@Catch[,1,,,]
  hCr = getChist(MSE)
  nr = dim(hCr)[3]
  AC = array(NA, c(nsim, nMP, nr, at))  # SMT
  for(MP in 1:nMP)AC[,MP,,1:nt] = aperm(hCr,c(1,3,2)) # historical C by MP fill
  AC[,,,nt+1:pt] = aperm(C,c(1,3,2,4))                 # projected C
  AC
}


mahiplot_area = function(MSE,iqr = 0.8, MPcols = c("red","green","blue","orange","grey","purple"),npy=20, MPnams=NA){
  nsubyr = 4
  if(is.na(MPnams[1]))MPnams = unlist(MSE@MPs)
  CurrentYr = MSE@OM[[1]][[1]]$CurrentYr[1]
  nsim = MSE@nsim; nMP = MSE@nMPs
  nt = MSE@nyears;  pt = MSE@proyears;   at = nt+pt;   ny = nt/nsubyr; py = pt/nsubyr; ay = at/nsubyr
  yind = rep(1:ny, each=4)
  ylab = CurrentYr + (-ny +1):py

  Cp = get_cat_area(MSE)/1E6 # s mp, r, time
  nr = dim(Cp)[3]
  Cy = apply(array(Cp, c(nsim, nMP, nr, nsubyr, ay)),c(1,2,3,5),sum)
  mus = apply(Cy,2:4,mean) # MP, nr, ay

  par(mfrow=c(ceiling(nr/2),2),mai=c(0.25,0.35,0.03,0.03),omi=c(0.3,0.3,0.01,0.01))
  for(rr in 1:nr){

    mahiqplot(mu = mus[,rr,],sims= Cy[,,rr,], xval = ylab, vline = CurrentYr+0.5, ny, MPcols,ylab = "",npy,iqr=iqr)
    legend('topright',legend=Areas[rr],bty="n",cex=0.85)
    if(rr ==1)legend('topleft',legend=MPnams,text.col=MPcols,text.font=2,cex=0.85,bty='n')

  }
  mtext("Year",1,outer=T,line=0.4)
  mtext("Yield (kt)",2,outer=T,line=0.4)
}

mahiplot_fleet = function(MSE){




}
