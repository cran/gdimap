\name{data.bvec}
\alias{data.bvec}
\docType{data}
\title{Shell Grid File}
\description{
Shell grid file used in diffusion data acquisition for file \file{data.nii.gz}.
}

\format{
	The file \file{data.nii.gz} is in gzipped NIfTI format.
	The R-package \pkg{oro.nifti} is required to read NIfTI files.
}

\details{
The S2 shell has 203 points uniformly distributed on a 3D shell limited to the volume of the unit sphere.
The data set \file{data.nii.gz} is a DICOM data set provided by the \emph{Advanced Biomedical MRI Lab, National Taiwan University Hospital}, under the designation \file{DSI 203-point 3mm}.
}

\references{
Yeh, F.-C., Wedeen, V. J., and Tseng, W.-Y. I. Generalized q-Sampling Imaging. \emph{IEEE Transactions on Medical Imaging} 29, 9 (2010),
1626--1635.
}

\keyword{ data }

