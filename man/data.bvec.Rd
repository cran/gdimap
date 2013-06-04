\name{data.bvec}
\alias{data.bvec}
\docType{data}
\title{3D b-table Vectors for Data Acquisition}
\description{
3D grid file containing the sampling points used in diffusion data acquisition for file \file{data.nii.gz}.
}

\format{
	\file{data.bvec} is in ASCII format, specified as a sequence of values,
  as produced by the 'dcm2nii' tool, and compatible with the formats used
  in the 'DSI_Studio' and 'Dipy' toolboxes.
}

\details{
The file \file{data.bvec} has 203 points uniformly distributed on a 3D (DSI) grid limited to the volume of the unit sphere.
The data set \file{data.nii.gz} is  NIfTI-converted from the DICOM data set provided by the \emph{Advanced Biomedical MRI Lab, National Taiwan University Hospital}, under the designation \file{DSI 203-point 3mm}.
}

\references{
Yeh, F.-C., Wedeen, V. J., and Tseng, W.-Y. I. Generalized q-Sampling Imaging. \emph{IEEE Transactions on Medical Imaging} 29, 9 (2010),
1626--1635.
}

\keyword{ data }

