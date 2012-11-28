
rgbvolmap <-
function(fbase=NULL, rg=c(1,1), bview="coronal",
 texture=FALSE, transparent=FALSE,  saveplot=FALSE,
 pngfile=paste(tempdir(),"/rgbmap.png",sep="") )
{
  bviews <- c("sagittal", "coronal", "axial")
  kv <- grep(bview, bviews)
  ##----------------------------
  img.nifti  <- readniidata(fbase=fbase, filename="data_gfa_gqi.nii.gz")
  gfavol <- img.nifti@.Data
  v1.nifti  <- readniidata(fbase=fbase, filename="data_V1_gqi.nii.gz")
  V1 <- v1.nifti@.Data
  ##----------------------------
  d <- dim(gfavol)
  switch(kv,
    { nr <- d[2]; nc <- d[3]}, # sagittal,
    { nr <- d[1]; nc <- d[3]}, # coronal
    { nr <- d[1]; nc <- d[2]})  # axial
  if(is.null(rg)) {
    switch(kv,
      { nslices <- d[1]}, # sagittal,
      { nslices <- d[2]}, # coronal
      { nslices <- d[3]})  # axial
    first <- 1; last <- nslices
    rg <- c(first,last)
  }
  else { first <- rg[1]; last <- rg[2] }
  ## Initialize Display
  if(transparent)
    x11()
  else
    x11(bg="black")
  displist <- vector("list",diff(rg)+1)
  j <- 0
  for (sl in (first:last)) {
    cat(sl,"")
    j <- j+1
    #---------------
    # RGB channels
    switch(kv,
      gfaslice <- gfavol[sl,,], # sagittal
      gfaslice <- gfavol[,sl,], # coronal
      gfaslice <- gfavol[,,sl]) # axial
    ix <- which(gfaslice != 0, arr.ind=TRUE)
    gfas <- gfaslice[ix]
    tmp <- matrix(0, nrow=nr, ncol=nc)  
    switch(kv,
      vtmp <- V1[sl,,,], # sagittal
      vtmp <- V1[,sl,,], # coronal
      vtmp <- V1[,,sl,]) # axial
    v <- matrix(0,nrow=dim(ix)[1],  ncol=3)
    for(k in 1:3) {
      va <- vtmp[,,k] 
      v[,k] <- va[ix]
    }
    cvm <- gfas*abs(v) ## choose ONE color as GFA*|Vpeak| 
    zp <- array(0, dim=c(nr,nc,3))
    tmp <- matrix(0, nrow=nr, ncol=nc)  
    for(k in 1:3) {
      tmp[ix] <- cvm[,k]
      zp[,,k] <- tmp 
      tmp[ix] <- 0
    }
    if(max(zp) > 0)
      zp <- zp/max(zp)
    zpr <- as.raster(zp)
    ##----------
    ## display
    mask <- gfaslice
    if(transparent) 
      zpr[ mask == 0  ] <- "transparent"
    zpr2 <- array(0, dim=c(nc,nr))
    if(texture){ ## for use as texture
      nm <- max(nc,nr)
      zpr2 <- array(0, dim=c(nm,nm)) # using square texture
      zpr2[nc:1, ] <- zpr[ ,nc:1] 
      saveplot <- TRUE
    }
    else {
      zpr2 <- array(0, dim=c(nc,nr))
      zpr2[1:nc, ] <- zpr[ ,nc:1]
    }
    r1 <- rasterGrob(zpr2, interpolate=TRUE)
    displist[j] <- list(r1)
  }
  cat("\n")
  ## display gmaps 
  do.call(grid.arrange,  displist) 
  if(saveplot) {
    savePlot(pngfile)
    cat("saved",pngfile,"\n")
    if(texture){ Sys.sleep(0.25); dev.off() }
  }
}

