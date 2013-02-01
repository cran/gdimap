gqi.odfvmf <-
function(gdi="gqi", run=TRUE, fbase=NULL, rg=NULL, swap=FALSE, lambda=NULL, depth=3, btoption=2,
 threshold=0.4, showglyph=FALSE, bview="coronal", savedir=tempdir())
{
  gdimethods <- c("gqi", "gqi2")
  gdimethod <- match(gdi, gdimethods)
  bviews <- c("sagittal", "coronal", "axial")
  kv <- match(bview, bviews)
  stopifnot(is.na(kv) != TRUE)
  ## movMF options
  startctl=list(E="softmax", minalpha=8, start="s", maxiter=200) 
  ## generate S2 grid
  s2 <- s2tessel.zorder(depth=depth, viewgrid=FALSE)
  odfvertices <- s2$pc
  tcsurf <- s2$tcsurf
  ##----------------
  ## Read data
  testfilexist(fbase=fbase, btoption=btoption)
  if(btoption == 1){ ## Option 1: S2-shell (DSI 203-point 3mm)
    btable <- as.matrix(readtable(fbase=fbase, filename="btable.txt"))
  } else {
    if(btoption == 2) { ## Option 2: 3D-dsi grid 
      bval <- scantable(fbase=fbase, filename="data.bval")
      # bvec <- readtable(fbase=fbase, filename="data.bvec")
      bvec <- scantable(fbase=fbase, filename="data.bvec")
      bvec <- matrix(bvec, ncol=3)
      btable <- cbind(bval,bvec)
      rm(bval, bvec)
    }
    else stop()
  }
  gc()
  cat("Reading data ...\n") 
  ptm <- proc.time()
  img.nifti  <- readniidata(fbase=fbase, filename="data.nii.gz")
  volimg <- img.nifti@.Data  
  mask.nifti <- readniidata(fbase=fbase, filename="data_brain_mask.nii.gz")
  volmask <- mask.nifti@.Data  
  print(proc.time() - ptm)
  rm(img.nifti, mask.nifti)
  gc()
  ##----------------
  d <- dim(volmask)
  volgfa <- array(0, dim=d)   ## gfas map
  V1 <- array(0, dim=c(d, 3)) ## V1 direction
  V2 <- array(0, dim=c(d, 3)) ## V2 direction
  if(is.null(rg)) {
    switch(kv,
      { nslices <- d[1]}, # sagittal,
      { nslices <- d[2]}, # coronal
      { nslices <- d[3]})  # axial
    first <- 1; last <- nslices
  }
  else { first <- rg[1]; last <- rg[2] }
  cat("\n")
  ##-----------------------------
  ## "gdimethod" process
  cat("Estimating slice odfs ...\n")
  switch(gdimethod,
    q2odf <- gqifn(odfvert=odfvertices, btable=btable,
                   lambda=lambda),
    q2odf <- gqifn2(odfvert=odfvertices, btable=btable,
                   lambda=lambda) )
  ##-----------------------------
  ## store 1st vector directions for each non-thresholded voxel 
  ## v1list: vector of lists
  nv1 <- length(first:last)
  v1list <- vector(mode="list", nv1)
  v1count <- 0
  # rglinit()
  npar1 <- 7
  npar2 <- 15
  for (sl in (first:last)) {
    cat("slice",sl,"\n")
  	ptm <- proc.time()
    if(run) {
      slicedata <- read.slice(img=volimg, mask=volmask, slice=sl,
       swap=swap, bview=bview)
      ymaskdata <- premask(slicedata)
      if(ymaskdata$empty) next # empty mask
      ## odfs
      odfs <- q2odf %*% (ymaskdata$yn)
      odfs <- apply(odfs, 2, norm01) ## normalize 
      ## gfas
      gfas <- apply(odfs, 2, genfa)
      gfas <- norm01(gfas) ##?
      z2d <- ymaskdata$kin
      ## mask out thresholded values
      zx <- which(gfas <= threshold)
      if(length(zx)) {
        z2d <- z2d[-zx,]
        gfas <- gfas[-zx]
        odfs <- odfs[,-zx]
      }
      if(is.null(dim(z2d))) next
      if(length(gfas) < 8) next 
      lix <- dim(z2d)[1]
      v1perslice <- matrix(0, nrow=lix,ncol=3) # v1 directions 
      v2perslice <- matrix(0, nrow=lix,ncol=3) # v2 directions 
      nullvectors <- NULL
      tperc <- c(20, 40, 60, 80)
      tline <- floor(c(0.2,0.4,0.6,0.8)*lix) 
      cat("processing", lix,"voxels \n")
      for(m in 1:lix) {
        tt <- which(tline == m)
        if(length(tt) != 0) {
          cat(paste(tperc[tt],"% ", sep="")); cflush() }
        odf <- odfs[,m]
        ## Find peaks based on clusters
        ith <- which(odf < threshold) 
        vx <- odfvertices[-ith,]
        n <- dim(vx)[1]
        ## Fit a vMF mixture  with k=2
        y1 <- movMF::movMF(vx, k=2, control=startctl)
        ## Inspect the fitted parameters:
        par1 <- logb(n)*npar1
        bic1 <- 2*logLik(y1) - par1
        ## Fit a vMF mixture  with k=4
        y2 <- movMF::movMF(vx, k=4, control=startctl)
        par2 <- logb(n)*npar2
        bic2 <- 2*logLik(y2) - par2
        if(bic1 >= bic2) {
          yy <- y1
          np <- 2
        }
        else {
          yy <- y2
          np <- dim(yy$theta)[1]
        }
        # no scaling
        # mx <- apply(abs(yy$theta), 1, max)
        # pcoords <- yy$theta/mx
        # pk <- list(np=np , pcoords=t(pcoords))
        pk <- list(np=np, pcoords=t(yy$theta))
        v1perslice[m,] <- pk$pcoords[,1]
        if(np == 4) {
          if(all.equal(abs(pk$pcoords[,1]), abs(pk$pcoords[,2]),
               toler=0.01) == TRUE) {
            v2perslice[m,] <- pk$pcoords[,3]
          }
          else {
            v2perslice[m,] <- pk$pcoords[,2]
          }
        }
        if(showglyph) {
          ## normalize for vizualization
          mx <- apply(abs(yy$theta), 1, max)
          pcoordsn <- yy$theta/mx
          pkn <- list(np=np , pcoords=t(pcoordsn))
          if(pk$np > 2) {
            plotglyph(odfs[,m], odfvertices, pkn, kdir=4)
            pp <- readline(
                "\nmore crossing-fiber glyphs  ? ('n' to exit) ") 
            if(pp == "n" ) { showglyph <- FALSE; }
            else { rgl.clear( type = "shapes" ) }
          }
        }
      }
      cat("100% completed\n")
      ## remove null pk vectors
      nvl <- lix
      nnv <- length(nullvectors)
      if(nnv > 0) {
        nvl <- nvl-nnv
        v1perslice <- v1perslice[-nullvectors,]
        v2perslice <- v2perslice[-nullvectors,]
        z2d <- z2d[-nullvectors,]
        gfas <- gfas[-nullvectors]
      }
      if(is.null(dim(z2d))) next
      ## V volumes
      for(k in 1:3) {
        switch(kv,
        { mx <- matrix(0, d[2],d[3])
        mx[z2d] <- v1perslice[,k]
        V1[sl,,,k] <- mx
        mx <- matrix(0, d[2],d[3])
        mx[z2d] <- v2perslice[,k]
        V2[sl,,,k] <- mx }, # sagittal
        { mx <- matrix(0, d[1],d[3])
        mx[z2d] <- v1perslice[,k]
        V1[,sl,,k] <- mx
        mx <- matrix(0, d[1],d[3])
        mx[z2d] <- v2perslice[,k]
        V2[,sl,,k] <- mx }, # coronal
        { mx <- matrix(0, d[1],d[2])
        mx[z2d] <- v1perslice[,k]
        V1[,,sl,k] <- mx
        mx <- matrix(0, d[1],d[2])
        mx[z2d] <- v2perslice[,k]
        V2[,,sl,k] <- mx } ) # axial
      }
      ## gfas volume
      fsave <- paste(savedir,"/vol",sl,".RData",sep="")
      switch(kv,
      { mx <- matrix(0, d[2],d[3])
      mx[z2d] <- gfas
      volgfa[sl,,] <- mx 
      res <- list(kv=kv, gfa=volgfa[sl,,],  v1=V1[sl,,,], v2=V2[sl,,,],
        file=fsave) }, # sagittal
      { mx <- matrix(0, d[1],d[3])
      mx[z2d] <- gfas
      volgfa[,sl,] <- mx
      res <- list(kv=kv, gfa=volgfa[,sl,],  v1=V1[,sl,,], v2=V2[,sl,,],
        file=fsave) }, # coronal
      { mx <- matrix(0, d[1],d[2])
      mx[z2d] <- gfas
      volgfa[,,sl] <- mx
      res <- list(kv=kv, gfa=volgfa[,,sl],  v1=V1[,,sl,], v2=V2[,,sl,],
        file=fsave) } ) # axial
      ##
      save(res, file=fsave)
      cat("wrote", fsave,"\n")
    } else {
      fsave <- paste(savedir,"/vol",sl,".RData",sep="")
      load(fsave)
      cat("loaded", fsave, "\n")
      switch(res$kv,
			{ V1[sl,,,] <- res$v1
        V2[sl,,,] <- res$v2
        volgfa[sl,,] <- res$gfa },
			{ V1[,sl,,] <- res$v1
        V2[,sl,,] <- res$v2
        volgfa[,sl,] <- res$gfa },
			{ volgfa[,,sl] <- res$gfa
        V1[,,sl,] <- res$v1
        V2[,,sl,] <- res$v2 } )
    }
  	print(proc.time() - ptm)
  }
  f <- paste(savedir,"/data_gfa",sep="")
  writeNIfTI(volgfa, filename=f, verbose=TRUE)
  cat("wrote",f,"\n")
  f <- paste(savedir,"/data_V1",sep="")
  writeNIfTI(V1, filename=f, verbose=TRUE)
  cat("wrote",f,"\n")
  f <- paste(savedir,"/data_V2",sep="")
  writeNIfTI(V2, filename=f, verbose=TRUE)
  cat("wrote",f,"\n")
}

