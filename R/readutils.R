## $Id: read.niftivol.R,v 1.2 2012/09/15 17:48:50 arfs Exp $
##
## Default data: see DSI_Studio
## data may be prefiltered by FSL 
##   niislicets: slice time-series
##    mask: slice mask


read.slice <-
function (img, mask, slice, swap=FALSE, reorient=FALSE, bview="axial") 
{
		bviews <- c("sagittal", "coronal", "axial")
		kv <- grep(bview, bviews)
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
       	{	niislicets <- img[, , slice , ]
        	mask <- mask[, , slice]} )
    }
    invisible(list(slice=slice, niislicets=niislicets, mask=mask, swap=swap))
}

##----------------------------------------------------------
readniidata <-
function(fbase=NULL, filename)
{
  if(is.null(fbase)) {
		file <- system.file(paste("extdata/", fbase, filename, sep = ""), package = "gdimap")
  }
	else {
		file <- paste(fbase,"/",filename, sep="")
	}
	options(niftiAuditTrail = FALSE)
  vol.nifti <- readNIfTI(file, reorient=FALSE)
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
#     kin <- which(mask == 1, arr.ind = T) # indices of pixels in mask 
    kin <- which(mask >= 1, arr.ind = T) # indices of pixels in mask 
    if (dim(kin)[1] < 2 ) { # minimum 2 dots in mask
    # if (!length(kin)) 
        # cat("\n slice", slice, ":\tempty slice mask - nothing to do\n")
#         return()
         return(list(empty=TRUE))
    }
    d <- dim(kin)
		yn <- matrix(0, nrow=dim(niislicets)[3], ncol=d[1])
    for (i in 1:d[1]) {
        yx <- niislicets[kin[i, 1], kin[i, 2], ]
        if (sd(yx)) {  # do not include null time series even if mask is 1 
						yn[,i] <- yx
        }
        else { # remove form mask 
            mask[kin[i, 1], kin[i, 2]] <- 0
        }
    }
#    kin <- which(mask == 1, arr.ind = T)  # update indices of pixels in mask
    kin <- which(mask >= 1, arr.ind = T)  # update indices of pixels in mask
		###
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
  if(is.null(fbase)) {
		file <- system.file(paste("extdata/", fbase, filename, sep = ""), package = "gdimap")
  }
	else {
		file <- paste(fbase,"/",filename, sep="")
	}
  invisible(scan(file))
}

##----------------------------------------------------------
readtable <-
function(fbase=NULL, filename)
{
  if(is.null(fbase)) {
		file <- system.file(paste("extdata/", fbase, filename, sep = ""), package = "gdimap")
  }
	else {
		file <- paste(fbase,"/",filename, sep="")
	}
  invisible(as.matrix(read.table(file)))
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
	colorlut <- terrain_hcl(12, h = c(0, -100), c. = c(40, 80), l = c(75, 40), power = 1)
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

#-------------------------
