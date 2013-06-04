\name{gdimap-package}
\Rdversion{0.1-1}
\alias{gdimap-package}
\docType{package}
\title{
Generalized Diffusion Magnetic Resonance Imaging
}

\description{
The package implements algorithms to estimate and visualize the orientation of neuronal pathways in model-free methods (q-space imaging methods).
For fibre orientation estimation based on mixtures of von Mises-Fischer (vMF) distributions see \code{\link{gqi.odfvmflines}}, \code{\link{gqi.odfvmf}}, \code{\link{sph.odfvmflines}}, and \code{\link{sph.odfvmf}}.
For fibre orientation estimation based on local peak detection procedures see \code{\link{gqi.odfpeaklines}}, \code{\link{gqi.odfpeaks}}, \code{\link{sph.odfpeaklines}}, and \code{\link{sph.odfpeaks}}.
In these names, the prefix \code{gqi} references functions using the standard GQI method (Yeh et.al. 2010), or the GQI2 variant (Garyfallidis 2012) for ODF reconstruction.
The prefix \code{sph} references functions using spherical harmonics approximations based on Aganj's ODF reconstruction (CSA-QBI).
Visualization of RGB-maps and glyph-maps is implemented through \code{\link{rgbvolmap}} and \code{\link{gqi.odfvxgrid}}, respectively.
Various types of simulations of diffusion and ODF profiles, illustrating the application of vMF mixtures for fibre orientation estimation are implemented in \code{\link{simulglyph.vmf}}, \code{\link{simul.simplefield}}, \code{\link{simul.fandtasia}}.
Data sets used in the examples are documented.
}

\details{
For a complete list of functions type \code{library(help=gdimap)}.
}


\note{
For reasons of limited CRAN storage space, a trimmed down data set with just 4 coronal slices was selected. Coronal slices (29:32) were extracted from the original data set. Hence, the argument \code{bview="coronal"} is used in the examples, and the range of coronal slices is (1:4). User specified complete data sets may be used without restrictions. 

Commonly used acronyms:\cr
GFA - Generalized Fractional Anisotropy\cr
ODF - Orientation Distribution Function\cr
GQI - Generalized q-Sampling Imaging\cr
GQI2 - Generalized q-Sampling Imaging variant\cr
QBI - Q-ball Imaging\cr
CSA-QBI - Constant Solid Angle QBI
MRI - Magnetic Resonance Imaging\cr
dMRI - Diffusion Magnetic Resonance Imaging\cr
RGB maps - Red:Green:Blue color maps\cr
HARDI - High Angular Resolution Diffusion Imaging\cr
vMF - von Mises-Fisher\cr
}

\references{
Ferreira da Silva, A. R. Facing the Challenge of Estimating Human Brain White Matter Pathways. In \emph{Proc. of the 4th International Joint Conference on Computational Intelligence} (Oct. 2012), K. Madani, J. Kacprzyk, and J. Filipe, Eds., SciTePress, pp. 709-714.

Hornik, K., and Gruen, B. \pkg{movMF}: Mixtures of von Mises-Fisher Distributions, 2012. R package version 0.1-0.

Hankin, R. K. S. Special functions in R: introducing the \pkg{gsl} package. R News 6 (October 2006).

Adler, D., and Murdoch, D. \pkg{rgl}: 3D visualization device system (OpenGL), 2012. R package version 0.92.880.

Auguie, B. \pkg{gridExtra}: functions in Grid graphics, 2012. R package version 0.9.1.

Barber, C. B., Habel, K., Grasman, R., Gramacy, R. B., Stahel, A., and Sterratt, D. C. \pkg{geometry}: Mesh generation and surface tesselation, 2012. R package version 0.3-2.

R Core Team. \pkg{R}: A Language and Environment for Statistical Computing.  R Foundation for Statistical Computing, Vienna, Austria, 2012. ISBN 3-900051-07-0.

Whitcher, B., Schmid, V. J., and Thornton, A. Working with the DICOM and NIfTI data standards in R. \emph{Journal of Statistical Software} 44, 6 (2011), 1-28.

Jenkinson, M., Beckmann, C. F., Behrens, T. E., Woolrich, M. W., and Smith, S. M. FSL. \emph{NeuroImage} 62, 2 (2012), 782-790.

Tuch D. S., Q-Ball Imaging, \emph{Magnetic Resonance in Medicine} 52 (2004), 1358-1372. 

Tabelow K., Polzehl J.: \pkg{dti}: DTI/DWI Analysis, 2012. R package version 1.1-0.

Polzehl J., Tabelow K., dti: Beyond the Gaussian Model in Diffusion-Weighted Imaging. \emph{Journal of Statistical Software} 44, 12 (2011).

Yeh, F.-C., Wedeen, V. J., and Tseng, W.-Y. I. Generalized q-Sampling Imaging. \emph{IEEE Transactions on Medical Imaging} 29, 9 (2010), 1626-1635.

Garyfallidis E., Towards an Accurate Brain Tractography, 2012, PhD Thesis, University of Cambridge.

Aganj I., Lenglet C., Sapiro G., Yacoub E., Ugurbill K. and Harel N. Reconstruction of the orientation distribution function in single- and multiple-shell q-ball imaging within constant solid angle, \emph{Magnetic Resonance in Medicine}, 64(2) (2010), 554-566.

}


\author{ Adelino Ferreira da Silva, Universidade Nova de Lisboa,
  Faculdade de Ciencias e Tecnologia, Portugal, \email{afs at fct.unl.pt}
}

\keyword{ package }

