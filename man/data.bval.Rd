\name{data.bval}
\alias{data.bval}
\docType{data}
\title{b-Table File}
\description{
b-table file for the 3D-DSI sampling scheme used in the DICOM data acquisition process for file \file{data.nii.gz}.
}

\format{
	The file \file{data.nii.gz} is in gzipped NIfTI format.
	The R-package \pkg{oro.nifti} is required to read NIfTI files.
}

\details{
The data set \file{data.nii.gz} is a DICOM data set provided by the \emph{Advanced Biomedical MRI Lab, National Taiwan University Hospital}, under the designation \file{DSI 203-point 3mm}.
This b-table has 203 points uniformly distributed on a 3D grid limited to the volume of the unit sphere.
The b-values range from 0 to 4000.
}


\references{
Yeh, F.-C., Wedeen, V. J., and Tseng, W.-Y. I. Generalized q-Sampling Imaging. \emph{IEEE Transactions on Medical Imaging} 29, 9 (2010),
1626--1635.
}

\keyword{ data }

