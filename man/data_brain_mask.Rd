\name{data_brain_mask}
\alias{data_brain_mask}
\docType{data}
\title{Example of Mask File Used in Diffusion MRI Processing}
\description{
The \file{data_brain_mask.nii.gz} dataset defines the mask used in the diffusion MRI analysis for \file{data.nii.gz}.
}

\format{
	The file \file{data_brain_mask.nii.gz} is in gzipped NIfTI format.
	The R-package \pkg{oro.nifti} is required to read NIfTI files.
	This file may be obtained from \file{data.nii.gz} by applying FSL/BET pre-processing tools.
}

\details{
For reasons of limited CRAN storage space, a trimmed down data set with just 4 coronal slices was selected. Coronal slices {29-32} were extracted from the original data set.

The mask used here is an all-brain mask; it just removes non-brain regions, as the result of applying a brain extraction tool to the \file{data.nii.gz} dataset.
However, other masks may be defined to select regions of interest (ROIs).
}

\references{
Yeh, F.-C., Wedeen, V. J., and Tseng, W.-Y. I. Generalized q-Sampling Imaging. \emph{IEEE Transactions on Medical Imaging} 29, 9 (2010),
1626--1635.

Jenkinson, M., Beckmann, C. F., Behrens, T. E., Woolrich, M. W., and Smith, S. M. Fsl. \emph{NeuroImage} 62, 2 (2012), 782-790.

NITRC repository,  \url{http://www.nitrc.org}.
}

\seealso{
  \code{\link{data}},
  \code{\link{data_brain}},
  \code{\link{data.bval}},
  \code{\link{data.bvec}},
  \code{\link{data_V1_gqi}},
  \code{\link{data_V2_gqi}},
  \code{\link{data_gfa_gqi}},
  \code{\link{dsi203_bmax4000.Rd}}
}

\keyword{ data }

