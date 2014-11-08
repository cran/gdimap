## Setup for writing NIfTI file using Rniftilib
##
niisetup <- 
function(savedir=tempdir(),filename, dim)
{
  niifile <- nifti.image.new()
  #f <- paste(savedir,paste(filename,"nii.gz", sep=""),sep="")
  f <- file.path(savedir, paste(filename,".nii.gz", sep=""))
  if(file.exists(f)) file.remove(f)
  nifti.set.filenames(niifile, f)
  niifile$dim <- dim
  # niifile$pixdim <- volimg$pixdim[1:3]
  invisible(niifile)
}

