\name{data}
\alias{data}
\docType{data}
\title{A Real Dataset for Diffusion MRI Analysis}
\description{
The data set \file{data.nii.gz} is a gzipped NIfTI data set converted from the original DICOM data set provided by the \emph{Advanced Biomedical MRI Lab, National Taiwan University Hospital}, under the designation \file{DSI 203-point 3mm}.

For reasons of limited CRAN storage space, a trimmed down data set with just 4 coronal slices was selected. Coronal slices (29:32) were extracted from the original data set. Hence, the argument \code{bview="coronal"} is used in the examples, and the range of coronal slices is (1:4). User specified complete data sets may be used without restrictions. 
}

\format{
	The file \file{data.nii.gz} is in gzipped NIfTI format.
	The R-package \pkg{oro.nifti} is required to read NIfTI files.
}

\details{

The data set is included in the \pkg{DSI Studio} package, publicly available from the NITRC repository.
This data set is from a normal 24-year-old male volunteer, and has been provided as a demonstration data set in connection with the \emph{DSI Studio} software for diffusion MR image analysis. 
The data set was obtained with an echo planar imaging diffusion sequence with twice-refocused echo, dimension 64x64x40, and slice thickness 2.9mm.
Further details on the data set specification are available from the NITRC repository.
}

\references{
Whitcher, B., Schmid, V. J., and Thornton, A. Working with the DICOM and NIfTI data standards in R. \emph{Journal of Statistical Software} 44, 6 (2011), 1-28.

Yeh, F.-C., Wedeen, V. J., and Tseng, W.-Y. I. Generalized q-Sampling Imaging. \emph{IEEE Transactions on Medical Imaging} 29, 9 (2010),
1626--1635.

NITRC repository,  \url{http://www.nitrc.org}.
}

\seealso{
  \code{\link{data_brain}},
  \code{\link{data_brain_mask}},
  \code{\link{data.bval}},
  \code{\link{data.bvec}},
  \code{\link{data_V1_gqi}},
  \code{\link{data_V2_gqi}},
  \code{\link{data_gfa_gqi}},
  \code{\link{dsi203_bmax4000}}
}


\keyword{ data }

