##
## Normalize PD values in NIfTI files.
##
niinorm <- 
function(srcdir=tempdir(), filename="data_V1",
         savedir=tempdir())
{
  normvf <- function(x) { norm(matrix(x,length(x),1),"f") }
  cat("Reading data ...\n")
  niifile  <- readniidata(fbase=srcdir, filename=filename)
  FEvol <- nifti.image.read(niifile)
  d <- dim(FEvol)
  FE <- array(0, dim=d)
  FE[] <- FEvol[]
  cat("Data normalization ...\n")
  nnx <- apply(FE, c(1,2,3), normvf)
  nnx <- drop(nnx)
  for(i in 1:d[4])
    FE[,,,i] <- FE[,,,i]/nnx
  zi <- which(is.nan(FE))
  FE[zi] <- 0
  f <- strsplit(filename,".nii")[[1]][1]
  f <- paste(f,"n",sep="")
  niifile <- niisetup(savedir=savedir, filename=f, dim=dim(FE))
  niifile[] <- FE[]
  nifti.image.write(niifile)
  cat("wrote",file.path(savedir,f),"\n")
}

