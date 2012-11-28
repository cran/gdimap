simul.simplefield <-
function(ang=60, b=3000, sigma=NULL, threshold=0.5)
{
  ## S2 shell grid
  s2 <- s2tessel.zorder(depth=3, viewgrid=FALSE)
  g0 <- s2$pc
  ## field simulation 
  field <- myglyph.synthsimul(g0, ang=ang, b=b, sigma=sigma)
  ## Estimate ODFs
  odfs <- fieldtestodf.gqi(g0, field,  b=b, mddratio=1.2)
  ## Visualize grid of glyphs with color
  plotodfvxgrid(g0, field=field, odfsdata=odfs)
  ## Using movMF (von Mises for peak detection)
  plotodflines.vmf.colors(g0, field=field, odfsdata=odfs,
    showglyph=FALSE, threshold=threshold)
}

##--------------------------

myglyph.synthsimul <-
function(g0, ang=60, b=3000, sigma=NULL)
{
  ng0 <- dim(g0)[1]
  stopifnot(ang <= 90)
  if(ang < 70 )
    a <- glyph.2x2()
  else
    a <- glyph.cross4()
  ## swap for visualization !!
  as <- a
  as[4:1,1:4] <- a 
  print(as)
  nn <- length(which(a != 0))
  ## synthesis
  S <- matrix(0, nn, ng0)
  k <- 0
  open3d()
  for(i in 1:4) {
    for(j in 1:4) {
      pos <- 2 * c(j,i,0)
      if(as[i,j] == 1) {
        sv <- synthfiberss2z(g0=g0, angles=0, b=b, sigma=sigma, pos=pos,
          showglyph=TRUE, new=FALSE)
        k <- k+1
        S[k,] <- sv
      }
      else {
        if(as[i,j] == 3) {
          sv <- synthfiberss2z(g0=g0, angles=c(0,ang), b=b, sigma=sigma,
          pos=pos, showglyph=TRUE, new=FALSE)
          k <- k+1
          S[k,] <- sv
        }
        else {
          if(as[i,j] == 2) {
            sv <- synthfiberss2z(g0=g0, angles=ang, b=b, sigma=sigma,
              pos=pos, showglyph=TRUE, new=FALSE)
            k <- k+1
            S[k,] <- sv
          }
        }
      }
    }
  }
  list(S=S, mask=a)
}

##--------------------------

glyph.2x2<-
function()
# function(ang=90, linefibs=TRUE)
{
  a <- matrix(2,4,4)
  a[2,1:4] <- 3
  a[3,1:4] <- 3
  a[1,1] <- 0
  a[1,2] <- 0
  a[4,3] <- 0
  a[4,4] <- 0
  a[2,1] <- 1
  a[2,4] <- 1
  a[3,1] <- 1
  a[3,4] <- 1
  ##
  # a[4:1,1:4] <- a # required  
  invisible(a)
}


glyph.cross4<-
function(ang=90)
{
  a <- matrix(2,4,4)
  a[2,1:4] <- 3
  a[3,1:4] <- 3
  a[1,1] <- 0
  a[1,4] <- 0
  a[4,1] <- 0
  a[4,4] <- 0
  a[2,1] <- 1
  a[2,4] <- 1
  a[3,1] <- 1
  a[3,4] <- 1
  ##
  invisible(a)
}

#--------------------------

fieldtestodf.gqi <-
function(grad, field, b=3000, mddratio=1.2)
{
  cat("estimating field odfs ...\n")
  sfield <- field$S
  mask <- field$mask
  dv <- dim(sfield)
  odfs <- matrix(0, dv[1], dv[2])
  bn <- rep(b, dim(grad)[1])
  btable <- cbind(bn, grad)
  q2odf <- gqifn(odfvert=grad, btable=btable, mddratio=mddratio)
  k <- 0
  for(i in 1:4) { 
    for(j in 1:4) {
      if(!mask[i,j]) next
      k <- k+1
      odf <- as.vector(q2odf%*%sfield[k,])
      odf <- odf - min(odf)
      odfs[k,] <- odf
    }
  }
  invisible(odfs)
}

#--------------------------

plotodfvxgrid <-
function(pc0, field, odfsdata)
{
  mask <- field$mask
  ## GFA
  odfs.reg <- t(apply(odfsdata, 1 , norm01))
  gfas <- apply(odfs.reg, 1, genfa)
  tc <-  delaunayn(pc0)
  tc.surf <- t( surf.tri(pc0,tc) )
  dt2 <- dim(pc0[tc.surf,])[1]
  d1 <- dim(odfsdata)[1]
  sgrid <- matrix(0, nrow=d1*dt2, ncol=3)
  vcolors <- matrix(0, nrow=d1, ncol=dt2)
  mask[1:4,4:1] <- mask # turn sweep upside-down 
  k <- 0
  cat("Running ...\n")
  for(i in 1:4) {
     for(j in 1:4) {
      if(!mask[i,j]) next
      k <- k+1
      odf <- odfs.reg[k,]
      gk <- gfas[k]
      ## RGB channels
      zch <- pc0*gk
      zch <- t(apply(abs(zch),1,norm01))
      ck <- rgb(zch)
      pc <- pc0 * as.vector(odf) 
      pc <- pc / (2*max(pc))
      # pc <- pc / max(pc)
      pos <- c(j,i,0)
      pcsurf <- cbind(
        pc[tc.surf,1]+pos[1], pc[tc.surf,2]+pos[2], pc[tc.surf,3]+pos[3])
      b <- (k-1)*dt2; e <- b+dt2
      sgrid[(b+1):e, ] <- pcsurf
      vcolors[k,] <- ck[tc.surf]
    }
  }
  cat("Plotting ...\n")
  # rgl.open()
  open3d()
  rgl.viewpoint(theta=0, phi=0)
  rgl.triangles(sgrid[,1], sgrid[,2], sgrid[,3], col=t(vcolors))
  rgl.viewpoint(0,0)
}


#--------------------------
#
# Overlay segments and GFA (anisotropy threshold)
#
plotodflines.vmf.colors <-
function(pc0, field, odfsdata, showglyph=FALSE, threshold=0.5)
{
  startctl=list(E="softmax", minalpha=6, start="s", maxiter=200) # movMF inits
  mask <- field$mask
  ## GFA
  odfs.reg <- t(apply(odfsdata, 1 , norm01))
  gfas <- apply(odfs.reg, 1, genfa)
  tc <-  delaunayn(pc0)
  tc.surf <- t( surf.tri(pc0,tc) )
  ## ------------
  d1 <- dim(odfsdata)[1]
  nn <- 8*d1
  v <- matrix(0, nrow=nn, ncol=3)
  ck <- numeric(nn)
  mask[1:4,4:1] <- mask # turn sweep 
  q <- 1
  m <- 0
  cat("running ... \n")
  for(i in 1:4) {
    for(j in 1:4) {
      if(!mask[i,j]) next
      m <- m+1
      odf <- odfs.reg[m,]
      gk <- gfas[m]
      ith <- which(odf < threshold) 
      vx <- pc0[-ith,]
      n <- dim(vx)[1]
      nc <- dim(vx)[2]
      kc <- 1:8
      npar <- nc*kc+kc-1
      bic <- -1.0e+10; nf <- 0; yy <- NULL
      for(k in seq(2,4,by=2)) {
        y2 <- movMF(vx, k=k, control=startctl) 
        par <- logb(n)*npar[k]
        bic2 <- 2*logLik(y2) - par
        if(bic2 > bic) {
          bic <- bic2
          nf <- k
          yy <- y2
        }  
      }
      np <- dim(yy$theta)[1]
      pcoords <- yy$theta/max(yy$theta)
      pk <- list(np=np , pcoords=t(pcoords))
      ##----------------------------
      ## 1st direction of max odf values for  odfs
      pc <- pc0 * as.vector(odf) 
      pc <- pc / (2*max(pc))
      pos <- c(j,i,0)
      coords <- pk$pcoords
      for(k in 1:min(pk$np, 4)) {
        zch <- coords[,k] * gk
        zch <- t(norm01(abs(zch)))
        ck[q] <- rgb(zch)
        ck[q+1] <- ck[q]
        pp <- pk$pcoords[,k]/2
        v[q,] <- pos
        v[q+1,]  <-  pp/2 + pos
        q <- q+2
      }
    }
  }
  cat("\n")
  open3d()
  cat("plotting ... \n")
  segments3d(v[1:(q-1),], col=ck[1:(q-1)], lwd=2, alpha=1)
  rgl.viewpoint(0,0)
}

