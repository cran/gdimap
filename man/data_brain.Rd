\name{data_brain}
\alias{data_brain}
\docType{data}
\title{Extracted Brain Data Set}
\description{
The \file{data_brain.nii.gz} dataset was obtained by applying the \emph{FST/BET brain extraction tool} to the diffusion data set \file{data.nii.gz}. \emph{BET} deletes non-brain tissue from an image of the whole head.
}

\format{
	The file \file{data_brain.nii.gz} is in gzipped NIfTI format.
	The R-package \pkg{oro.nifti} is required to read NIfTI files.
}

\details{
For reasons of limited CRAN storage space, a trimmed down data set with just 4 coronal slices was selected. Coronal slices (29:32) were extracted from the original data set. Hence, the argument \code{bview="coronal"} is used in the examples, and the range of coronal slices is (1:4). User specified complete data sets may be used without restrictions. 
}

\references{
Yeh, F.-C., Wedeen, V. J., and Tseng, W.-Y. I. Generalized q-Sampling Imaging. \emph{IEEE Transactions on Medical Imaging} 29, 9 (2010),
1626--1635.

Jenkinson, M., Beckmann, C. F., Behrens, T. E., Woolrich, M. W., and Smith, S. M. Fsl. \emph{NeuroImage} 62, 2 (2012), 782-790.

NITRC repository,  \url{http://www.nitrc.org}.

Whitcher, B., Schmid, V. J., and Thornton, A. Working with the DICOM and NIfTI data standards in R. \emph{Journal of Statistical Software} 44, 6 (2011), 1-28.
}

\seealso{
  \code{\link{data}},
  \code{\link{data_brain_mask}},
  \code{\link{data.bval}},
  \code{\link{data.bvec}},
  \code{\link{data_V1_gqi}},
  \code{\link{data_V2_gqi}},
  \code{\link{data_gfa_gqi}},
  \code{\link{dsi203_bmax4000}}
}

\keyword{ data }

