simul.fandtasia <-
function(gridsz=32, b=4000, depth=3, sigma=0.01, threshold=0.5, showglyph=FALSE, savedir=tempdir(), snapshot=FALSE, pngfig="fandtasia")
{
  ## S2 shell grid
  s2 <- s2tessel.zorder(depth=depth)
  g0 <- s2$pc
  simul.fandtasiaSignal(g=s2$pc, gridsz=gridsz, b=b, sigma=sigma,
    savedir=savedir)
  sfield <- readniidata(fbase=savedir, filename="simfield.nii.gz")
  ## Estimate ODFs
  fielddata <- sfield@.Data
  odfs <- field.gqi(g0, fielddata, b=b, mddratio=1.2)
  ## with vomMF clustering
  plotodfvmf(s2, odfs, threshold=threshold, showglyph=showglyph)
  if(snapshot) {
    sp <- paste(savedir,"/",pngfig,".png", sep="")
    Sys.sleep(2)
    rgl.snapshot(sp)
    cat("snapshot file", sp,"\n")
  }
}

##----------------------

field.gqi <-
function(grad, fielddata, b=4000, mddratio=1.2)
{
  cat("Estimating field odfs ...\n")
  sz <- dim(fielddata)[4]-1
  dv <- dim(fielddata)
  odfs <- array(0, dim=c(dv[1:3], dv[4]-1))
  bn <- rep(b, dim(grad)[1])
  btable <- cbind(bn, grad)
  q2odf <- gqifn(odfvert=grad, btable=btable, mddratio=mddratio)
  for(i in 1:dv[1]) { 
    for(j in 1:dv[2]) {
       S <- fielddata[i,j,1,]
      S <- S[2:dv[4]]
      odf <- q2odf%*%S
      odf <- odf - min(odf)
      odfs[i,j,1,] <- odf
    }
  }
  invisible(odfs)
}

##----------------------
# Overlay segments and GFA (anisotropy threshold)
plotodfvmf <-
function(s2, odfsdata, threshold=0.5, showglyph=FALSE) 
{
  startctl=list(E="softmax", minalpha=6, start="s", maxiter=200) # movMF inits
  odfvertices <- s2$pc
  dm <- dim(odfsdata)
  odfs.mat <- array(odfsdata, dim=c(dm[1]*dm[2], dm[4]))
  odfs.reg <- t(apply(odfs.mat, 1 , norm01))
  gfas <- apply(odfs.reg, 1, genfa)
  odfs.regarray <- array(odfs.reg, dim=dm)
  nn <- 8*dm[1]*dm[2]
  v <- matrix(0, nrow=nn, ncol=3)
  ck <- numeric(nn)
  m <- 1
  q <- 1
  npar1 <- 7
  npar2 <- 15
  tperc <- c(20, 40, 60, 80)
  tline <- floor(c(0.2,0.4,0.6,0.8)*dm[1])
  cat("vMF estimation for ", dm[1]*dm[2], "voxels, ...\n")
  for(i in 1:dm[1]) {
    tt <- which(tline == i)
    if(length(tt) != 0) {
      cat(paste(tperc[tt],"% ", sep="")); cflush() }
     for(j in 1:dm[2]) {
      odf <- odfs.regarray[i,j,1,]
      gk <- gfas[m]
      ith <- which(odf < threshold) 
      odf <- odf - min(odf)
      vx <- odfvertices[-ith,]
      n <- dim(vx)[1]
      y1 <- movMF(vx, k=2, control=startctl) 
      par1 <- logb(n)*npar1
      bic1 <- 2*logLik(y1) - par1
      y2 <- movMF(vx, k=4, control=startctl) 
      par2 <- logb(n)*npar2
      bic2 <- 2*logLik(y2) - par2
      if(bic1 >= bic2) {
        fitted <- predict(y1) }
      else {
        fitted <- predict(y2) }
      np <- length(unique(fitted))
      pcoords <- matrix(0,np,3)
      for(ii in 1:np) {
        ca <- which(fitted == ii)
        va <- vx[ca,]
        pcoords[ii,] <- colMeans(va)
      }
      pk <- list(np=np , pcoords=t(pcoords))
      ## optional visualization of crossing-fiber glyphs
      if(showglyph) {
        if(np == 4) {
          if(rgl.cur() == 0) open3d()
          plotglyph(odf, odfvertices, pk, kdir=4)
          pp <- readline(
            "\ncontinue showing crossing-fiber glyphs ? ('n' to exit) ") 
          if(pp == "n" ) { showglyph <- FALSE; }
          else { rgl.clear( type = "shapes" ) }
        }
      }
      ##----------------------------
      ## 1st direction of max odf values for  odfs
      pc <- odfvertices * as.vector(odf) 
      pc <- pc / (2*max(pc))
      pos <- c(i,j,0)
      coords <- pk$pcoords
      for(k in 1:min(pk$np, 4)) {
        zch <- coords[,k] * gk
        zch <- t(norm01(abs(zch)))
        ck[q] <- rgb(zch)
        ck[q+1] <- ck[q]
        pp <- pk$pcoords[,k]/2
        # v[q,] <- -pp + pos
        v[q,] <- pos
        v[q+1,]  <-  pp + pos
        q <- q+2
      }
      m <- m+1
    }
  }
  cat("100% completed\n")
  cat("\nplotting ... \n")
  open3d()
  segments3d(v[1:(q-1),], col=ck[1:(q-1)], lwd=2, alpha=1)
  rgl.viewpoint(0,0)
  par3d("windowRect"= c(100, 100, 700, 700), zoom=0.68)
  rgl.bringtotop()
}


