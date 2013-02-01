sph.odfvmf <-
function(run=TRUE, fbase=NULL, rg=NULL, swap=FALSE, btoption=1, threshold=0.4, showglyph=FALSE, bview="coronal", savedir=tempdir(), order=4)
{
  bviews <- c("sagittal", "coronal", "axial")
  kv <- match(bview, bviews)
  stopifnot(is.na(kv) != TRUE)
  ## movMF options
  startctl=list(E="softmax", minalpha=8, start="s", maxiter=200) 
  ##-----------
  ## Read data
  testfilexist(fbase=fbase, btoption=btoption)
  if(btoption == 1) { ## Option 1: S2-shell (DSI 203-point 3mm)
    btable <- as.matrix(
      readtable(fbase=fbase, filename="btable.txt"))
  }
  else {
    if(btoption == 2) { 
      if(is.null(fbase)) {
        cat("Data files 'data.bval' and 'data.bvec' unspecified !\n") 
        stop()
      } else {
        bval <- scantable(fbase=fbase, filename="data.bval")
        # bvec <- readtable(fbase=fbase, filename="data.bvec")
        bvec <- scantable(fbase=fbase, filename="data.bvec")
        bvec <- matrix(bvec, ncol=3)
        btable <- cbind(bval,bvec)
        rm(bval, bvec)
      }
    }
    else stop()
  }
  b0 <- which(btable[,1] == 0)
  odfvertices <- matrix(btable[-b0,2:4], ncol=3)
  tc <-  geometry::delaunayn(odfvertices)
  tcsurf <- t( surf.tri(odfvertices,tc))  
  ##----------------------------
  gc()
  cat("Reading data ...")
  img.nifti  <- readniidata(fbase=fbase, filename="data.nii.gz")
  volimg <- img.nifti@.Data  
  mask.nifti <- readniidata(fbase=fbase, filename="data_brain_mask.nii.gz")
  volmask <- mask.nifti@.Data  
  rm(img.nifti, mask.nifti)
  gc()
  ##----------------------------
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
  ## SPH process preparation
  gradient <- t(odfvertices)
  z <- design.spheven(order,gradient,lambda=0.006)
  plz <- plzero(order)/2/pi
  ngrad <- dim(gradient)[2]
  ngrad0 <- ngrad
  lord <- rep(seq(0,order,2),2*seq(0,order,2)+1)
  while(length(lord)>=ngrad0){
    order <- order-2
    lord <- rep(seq(0,order,2),2*seq(0,order,2)+1)
     cat("Reduced order of spherical harmonics to",order,"\n")
   }
  cat("Using",length(lord),"spherical harmonics\n")
  L <- -diag(lord*(lord+1)) 
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
    if(run) {
      slicedata <- read.slice(img=volimg, mask=volmask, slice=sl,
       swap=swap, bview=bview)
      ymaskdata <- premask(slicedata)
      if(ymaskdata$empty) next # empty mask
      maxslicedata <- max(slicedata$niislicets) ##????
      S <- ymaskdata$yn[-b0,]
       S <- S / maxslicedata
      s0 <- 1
      si <- apply(S, 2, datatrans, s0)
      sicoef <- z$matrix%*% si
      sphcoef <- plz%*%L%*%sicoef
      coef0 <- sphcoef[1,]
      sphcoef[1,] <- 1/2/sqrt(pi)
      sphcoef[-1,] <- sphcoef[-1,]/8/pi
      ## odfs
      odfs <- t(z$design) %*% sphcoef
      odfs <- apply(odfs, 2, norm01) 
      ## gfas
      gfas <- apply(odfs, 2, genfa)
      gfas <- norm01(gfas) ## ??
      z2d <- ymaskdata$kin
      ## mask out thresholded values
      zx <- which(gfas <= threshold)
      if(length(zx)) {
        z2d <- z2d[-zx,]
        gfas <- gfas[-zx]
        odfs <- odfs[,-zx]
      }
      if(is.null(dim(z2d))) next
      if(length(gfas) < 8) next # 8 elements as minimum number
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
          cat(paste(tperc[tt],"% ", sep="")); cflush()
        }
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
        mx <- apply(abs(yy$theta), 1, max)
        pcoords <- yy$theta/mx
        pk <- list(np=np , pcoords=t(pcoords))
        v1perslice[m,] <- pk$pcoords[,1]
        if(np == 4) {
          if(all.equal(abs(pk$pcoords[,1]), abs(pk$pcoords[,2]),
                 toler=0.01) == TRUE) 
            v2perslice[m,] <- as.numeric(pk$pcoords[,3])
           else 
            v2perslice[m,] <- as.numeric(pk$pcoords[,2])
        }
        if(showglyph) {
          if(pk$np > 2) {
          plotglyph(odfs[,m], odfvertices, pk, kdir=4)
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
  }
  cat("\n")
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

