\name{data_gfa}
\alias{data_gfa}
\docType{data}
\title{Generalized Fractional Anisotropy (GFA) File}
\description{
The \file{data_gfa} summarizes anisotropic properties of the Orientation Distribution Function (ODF) of the fibres using the generalized fractional anisotropy (GFA) metric, based on the GQI method.
}

\format{
	The file \file{data.nii.gz} is in gzipped NIfTI format.
	The R-package \pkg{oro.nifti} is required to read NIfTI files.
}

\details{
The file \file{data_gfa} was produced by the function \code{gqi.odfvmf}  (using the default arguments) included in \pkg{gdimap}.
Functions for visualization of RGB maps and line-maps are included in \pkg{gdimap}.
The output files \file{data_V1.nii.gz} and \file{data_gfa.nii.gz} are compatible with the "FSL/fslview" tool for RGB map and line map visualization.

For reasons of limited CRAN storage space, the source data set contains just 4 coronal slices as indicated in \code{\link{data}}.
}

\references{
Ferreira da Silva, A. R. Computational Representation of White Matter Fiber Orientations, \emph{International Journal of Biomedical Imaging}, Vol. 2013, Article ID 232143, Hindawi Publishing Corporation \url{http://dx.doi.org/10.1155/2013/232143}.

Ferreira da Silva, A. Facing the challenge of estimating human brain white matter pathways. \emph{In Proc. of NCTA 2012 - Int. Conf. on Neural Computation Theory and Applications} (Oct. 2012), K. Madani, J. Kacprzyk, and J. Filipe, Eds., SciTePress.

Whitcher, B., Schmid, V. J., and Thornton, A. Working with the DICOM and NIfTI data standards in R. \emph{Journal of Statistical Software} 44, 6 (2011), 1-28.

Jenkinson, M., Beckmann, C. F., Behrens, T. E., Woolrich, M. W., and Smith, S. M. Fsl. \emph{NeuroImage} 62, 2 (2012), 782-790.
}

\seealso{
  \code{\link{data_V1}},
  \code{\link{data_V2}}
}

\keyword{ data }

