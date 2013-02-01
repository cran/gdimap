
simulglyph.vmf <-
function(gdi="gqi", s2grid=NULL, angles=c(20,100), depth=3, b=3000, lambda=NULL, order=4, sigma=NULL, threshold=0.4, snapshot=FALSE, savedir=tempdir(), pngfig="glyph1", showglyph=TRUE, aniso=NULL, logplot=TRUE, wi=NULL)
{
  gdimethods <- c("gqi", "gqi2", "sph")
  gdimethod <- match(gdi, gdimethods)
  stopifnot(is.na(gdimethod) != TRUE)
  if(is.null(s2grid)) {
    ## S2 grid
    s2 <- s2tessel.zorder(depth=depth)
    g0 <- s2$pc
  }
  else 
    g0 <- s2grid
  ## synthetize diffusion signal
  ## open3d()
	if(is.null(wi)) {
	  l <- length(angles)
	  wi <- rep(1/l, times=l)
	}
	else 
		stopifnot(length(wi) == length(angles))
  ##
  S <- synthfiberss2z(g0=g0, angles=angles, logplot=logplot, b=b, S0=1,
    sigma=sigma, showglyph=showglyph, wi=wi)
  if(snapshot) {
    ## sp <- tempfile(pattern="gsim", tmpdir=savedir, fileext=".png")
    sp <- paste(savedir,"/",pngfig,"a.png", sep="")
    rgl.bringtotop();  rgl.snapshot(sp)
    cat("snapshot file", sp,"\n")
  }
  ##-------------------
  ## ODF reconstruction 
  # cat("Estimating slice odfs ...\n")
  switch(gdimethod,
    { br <- rep(b, length(S))
      btable <- cbind(br, g0)
      q2odf <- gqifn(odfvert=g0, btable=btable, lambda=lambda)
      odf <- q2odf%*%S }, # GQI
    { br <- rep(b, length(S)) 
      btable <- cbind(br, g0)
      q2odf <- gqifn2(odfvert=g0, btable=btable, lambda=lambda)
      odf <- q2odf%*%S }, # GQI@
    {  odf <- sphodf(g0=g0, S=S, order=order, lambda=lambda) } # SPH.w
  )
  ##--------------------
  if(!is.null(aniso)) {
    stopifnot(aniso >= 0 & aniso < 1)
    mx <- max(odf)  
    mn <- min(odf)  
    odf <- pmax(odf-mn-aniso*(mx-mn),0)
  }
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

sphodf <-
function(g0, S, order=4, lambda=0.006) 
{
  ## SPH process preparation
  cat("Estimating slice odfs ...\n")
  # gradient <- t(odfvertices)
  gradient <- t(g0)
  ## SPH: get SH design matrix ...

  z <- design.spheven(order,gradient,lambda=lambda) # OK
  # z <- design.wspheven(order,gradient,lambda=lambda) ## blobbier!!

  plz <- plzero(order)/2/pi
  ngrad <- dim(gradient)[2]
  ngrad0 <- ngrad
  lord <- rep(seq(0,order,2),2*seq(0,order,2)+1)
  while(length(lord)>=ngrad0){
     order <- order-2
     lord <- rep(seq(0,order,2),2*seq(0,order,2)+1)
     cat("Reduced order of spherical harmonics to",order,"\n")
  }
  cat("Using",length(lord),"sperical harmonics\n")
  L <- -diag(lord*(lord+1))
  ## SPH process 
  sicoef <- z$matrix%*% S
  sphcoef <- plz%*%L%*%sicoef
  coef0 <- sphcoef[1,]
  sphcoef[1,] <- 1/2/sqrt(pi)
  sphcoef[-1,] <- sphcoef[-1,]/8/pi
  sphcoef <- plzero(order)%*%sicoef
  odf <- t(z$design) %*% sphcoef#####
  invisible(odf)
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
    y2 <- movMF::movMF(vx, k=k, control=startctl) 
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


