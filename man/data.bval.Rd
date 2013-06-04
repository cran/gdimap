\name{data.bval}
\alias{data.bval}
\docType{data}
\title{b-Table File}
\description{
b-values file for the 3D-DSI sampling scheme used in the DICOM data acquisition process for file \file{data.nii.gz}.
}

\format{
	\file{data.bval} is in ASCII format
}

\details{
The b-values in \file{data.bval} range from 0 to 4000.
The data set \file{data.nii.gz} is NIfTI-converted from the DICOM data set provided by the \emph{Advanced Biomedical MRI Lab, National Taiwan University Hospital}, under the designation \file{DSI 203-point 3mm}.
The corresponding file \file{data.bvec} has 203 points uniformly distributed on a 3D (DSI) grid limited to the volume of the unit sphere.
}


\references{
Yeh, F.-C., Wedeen, V. J., and Tseng, W.-Y. I. Generalized q-Sampling Imaging. \emph{IEEE Transactions on Medical Imaging} 29, 9 (2010),
1626--1635.
}

\keyword{ data }

