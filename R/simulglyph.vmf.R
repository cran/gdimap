simulglyph.vmf <-
function(s2grid=NULL, angles=c(20,100), depth=3, b=3000, mddratio=1.24, sigma=NULL, threshold=0.4, snapshot=FALSE, savedir=tempdir(), pngfig="glyph1", showglyph=TRUE)
{
  if(is.null(s2grid)) {
    ## S2 grid
    s2 <- s2tessel.zorder(depth=depth)
    g0 <- s2$pc
  }
	else 
    g0 <- s2grid
  ## synthetize diffusion signal
  ## open3d()
  S <- synthfiberss2z(g0=g0, angles=angles, logplot=TRUE, b=b, S0=1,
    sigma=sigma, showglyph=showglyph)
  if(snapshot) {
    ## sp <- tempfile(pattern="gsim", tmpdir=savedir, fileext=".png")
    sp <- paste(savedir,"/",pngfig,"a.png", sep="")
    rgl.bringtotop();  rgl.snapshot(sp)
    cat("snapshot file", sp,"\n")
  }
  ##
  ## ODF reconstruction 
  br <- rep(b, length(S))
  btable <- cbind(br, g0)
  # cat("Estimating slice odfs ...\n")
  q2odf <- gqifn(odfvert=g0, btable=btable, mddratio=mddratio)
  odf <- q2odf%*%S
  odf <- odf - min(odf)
  odf <- (odf - min(odf))/diff(range(odf)) 
  ## vmf maxima estimation
  a <- gdi.vmf(odf=odf, odfvertices=g0, threshold=threshold,
    showglyph=showglyph)
  if(snapshot) {
    ## sp <- tempfile(pattern="gvmf", tmpdir=savedir, fileext=".png")
    sp <- paste(savedir,"/",pngfig,"b.png", sep="")
    rgl.bringtotop();  rgl.snapshot(sp)
    cat("snapshot file", sp,"\n")
  }
  invisible(a)
}


gdi.vmf <- 
function(odf, odfvertices, threshold=0.5, showglyph=TRUE)
{
  startctl=list(E="softmax", minalpha=8, start="s", maxiter=200) # movMF inits
  ith <- which(odf < threshold) 
  vx <- odfvertices[-ith,]
  n <- dim(vx)[1]
  nc <- dim(vx)[2]
  kc <- 1:8
  npar <- nc*kc+kc-1
  bic <- -1.0e+10; nf <- 0;; yy <- NULL
  for(k in seq(2,6,by=2)) {
    y2 <- movMF(vx, k=k, control=startctl) 
    par <- logb(n)*npar[k]
    bic2 <- 2*logLik(y2) - par
    if(bic2 > bic) {
      bic <- bic2
      nf <- k
      yy <- y2
    }  
  }
  # print(nf)
  np <- dim(yy$theta)[1]
  pcoords <- yy$theta/max(yy$theta)
  pk <- list(np=np , pcoords=t(pcoords))
  if(showglyph) {
    open3d()
    plotglyph(odf, odfvertices, pk, kdir=6)
  }
  ##-----------
  ## min. angle in rad
  a <- 2*pi
  for(i in 2:pk$np) {
    a1 <- anglev(pcoords[1,],pcoords[i,])
    if(a1 < a) a <- a1
  }
  invisible(c(np, a)) # a=min angle rad
  # invisible(c(np, a*180/pi)) # a=min angle degrees
}

anglev <-
function(a , b)
{
  normvf <- function(x) { norm(matrix(x,length(x),1),"f") }
  cn2 <- normvf(a)*normvf(b)
  cpn2 <- as.numeric(crossprod(a,b)/cn2)
  if(cpn2 >= 0)
    cpn2 <- min(cpn2,1)
  else
    cpn2 <- max(cpn2,-1)
  theta <- Re(acos(cpn2))
  invisible(theta) # in rad
}


