\name{data_V2_gqi}
\alias{data_V2_gqi}
\docType{data}

\title{ODF Second Principal Directions File}
\description{
The \file{data_V2_gqi} contains the second principal directions of the Orientation Distribution Function (ODF) of the fibres at each \file{data} voxel.
}

\format{
	The file \file{data.nii.gz} is in gzipped NIfTI format.
	The R-package \pkg{oro.nifti} is required to read NIfTI files.
}

\details{
The file \file{data_V2_gqi} was produced by the function \code{gqi.odfvmf}  (using the default arguments) included in \pkg{gdimap}.

For reasons of limited CRAN storage space, the source data set contains just 4 coronal slices as indicated in \code{\link{data}}.
}

\references{
Ferreira da Silva, A. Facing the challenge of estimating human brain white matter pathways. \emph{In Proc. of NCTA 2012 - Int. Conf. on Neural Computation Theory and Applications} (Oct. 2012), K. Madani, J. Kacprzyk, and J. Filipe, Eds., SciTePress.

Jenkinson, M., Beckmann, C. F., Behrens, T. E., Woolrich, M. W., and Smith, S. M. Fsl. \emph{NeuroImage} 62, 2 (2012), 782-790.
}

\seealso{
  \code{\link{data_V1_gqi}},
  \code{\link{data_gfa_gqi}},
}
\keyword{ data }

