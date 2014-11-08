## $Id: readutils.R,v 1.1 2014/07/16 15:06:52 arfs Exp arfs $
##
## niislicets: slice time-series
## mask: slice mask

read.slice <-
function (img, mask, slice, swap=FALSE, reorient=FALSE, bview="axial") 
{
  bviews <- c("sagittal", "coronal", "axial")
  kv <- match(bview, bviews)
  stopifnot(is.na(kv) != TRUE)
  X <- nrow(img)
  Xm <- nrow(mask)
  if (swap) { ## swap=T to be represented as in fslview
    switch(kv, ##?
    { niislicets <- img[slice, Xm:1, , ]
      mask <- mask[slice, Xm:1, ]}, 
    { niislicets <- img[Xm:1, slice, , ]  
      mask <- mask[Xm:1, slice, ]},
    { niislicets <- img[X:1, , slice, ]
      mask <- mask[Xm:1, , slice]})
  }
  else {
    switch(kv,
    { niislicets <- img[slice, , , ]
      mask <- mask[slice, , ]}, 
    { niislicets <- img[, slice , , ]  
      mask <- mask[, slice , ]},
     {  niislicets <- img[, , slice , ]
      mask <- mask[, , slice]} )
  }
  invisible(list(slice=slice, niislicets=niislicets, mask=mask, swap=swap))
}

##----------------------------------------------------------
readniidata <-
function(fbase=NULL, filename)
{
  if(is.null(fbase)) {
    file <- system.file(file.path("extdata", filename), package = "gdimap")
  }
  else {
	  file <- file.path(fbase,filename)
  }
	invisible(file)
}


##--------------------------------------
#
# Mask out slice times series and keep indices
#
premask <-
function (slicedata) 
{
  slice <- slicedata$slice
  niislicets <- slicedata$niislicets
  mask <- slicedata$mask
  # kin <- which(mask == 1, arr.ind = T) # indices of pixels in mask 
  kin <- which(mask >= 1, arr.ind = T) # indices of pixels in mask 
  kin <- matrix(kin, ncol=2)
  d <- dim(kin)
  if (d[1] < 2 ) { # minimum 2 dots in mask
    return(list(empty=TRUE)) }
  yn <- matrix(0, nrow=dim(niislicets)[3], ncol=d[1])
  rx <- numeric(d[1])
  ## do not include null time series even if mask is 1 
  for (i in 1:d[1]) {
    yx <- niislicets[kin[i, 1], kin[i, 2], ]
    if (sd(yx)) 
      yn[,i] <- yx
    else 
      rx[i] <- 1
  }
  ri <- which(rx == 1)
  if(length(ri) != 0) {
    mask[kin[ri, 1], kin[ri, 2]] <- 0
    kin <- kin[-ri,]
    yn <- yn[,-ri]
  }
  # stdf <- function(y) { return((y - mean(y))/sd(y)); }
  # yn <- apply(yn, 2, stdf)
  nobs <- slicedata$nobs
  stopifnot(nobs == nrow(yn))
  nreg <- ncol(yn)
  invisible(list(yn = yn, kin = kin, nreg = nreg, empty=FALSE))
}

##----------------------------------------------------------
scantable <-
function(fbase=NULL, filename)
{
  if(is.null(fbase)) 
    file <- system.file(file.path("extdata", filename), package = "gdimap")
  else 
    file <- file.path(fbase,filename)
  invisible(scan(file))
}

##----------------------------------------------------------
readtable <-
function(fbase=NULL, filename)
{
  if(is.null(fbase)) 
    file <- system.file(file.path("extdata",filename), package = "gdimap")
  else 
    file <- file.path(fbase,filename)
  invisible(as.matrix(read.table(file)))
}


#-------------------------
testfilexist <-
function(fbase=getwd(), btoption=2)
{
  getfilename <- function(filename, fbase=getwd()) {
    if(is.null(fbase)) 
      file <- system.file(file.path("extdata",filename), package = "gdimap")
    else 
      file <- file.path(fbase,filename)
    invisible(file)
  }
  data.nii.gz <- getfilename(filename="data.nii.gz", fbase=fbase)
  data.nii <- getfilename(filename="data.nii", fbase=fbase)
  stopifnot(file.exists(data.nii.gz) | file.exists(data.nii))
  data_brain_mask.nii.gz <- getfilename(filename="data_brain_mask.nii.gz", fbase=fbase)
  data_brain_mask.nii <- getfilename(filename="data_brain_mask.nii", fbase=fbase)
  stopifnot(file.exists(data_brain_mask.nii.gz) | file.exists(data_brain_mask.nii))
  ##
  if(btoption == 1) {
  btable.txt <- getfilename(filename="btable.txt", fbase=fbase)
  stopifnot(file.exists(btable.txt))
  }
  if(btoption == 2) {
    data.bvec <- getfilename(filename="data.bvec", fbase=fbase)
    stopifnot(file.exists(data.bvec))
    data.bval <- getfilename(filename="data.bval", fbase=fbase)
    stopifnot(file.exists(data.bval))
  }
}

#-------------------------
rglstart <- 
function(bg="white")
{
  rgl.open()
  rgl.clear("all")
  # colorlut <- terrain.colors(12) # height color lookup table
  # rgl.bg(color=c(colorlut[1],"white"))
  # rgl.light()
  colorlut <- colorspace::terrain_hcl(12, h = c(0, -100), c. = c(40, 80), l = c(75, 40), power = 1)
  bg3d(col=bg)
  light3d()  
}

#-------------------------
rglinit <- 
function()
{
  rgl.open()
  rgl.clear("all")
  rgl.bg(color="white")
  rgl.light()
}

#-------------------------
cflush <-
function() 
{
  if (Sys.info()[1] == "Windows") flush.console()
  return()
}

