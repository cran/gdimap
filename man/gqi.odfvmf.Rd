\name{gqi.odfvmf}
\alias{gqi.odfvmf}
\concept{Diffusion Magnetic Resonance}
\concept{GQI Reconstruction}
\concept{von Mises distributions}
\concept{Orientation Distribution Function}
\concept{RGB maps}


\title{GQI Reconstruction and Fibre Orientation Estimation Based on von Mises Distributions}
\description{
Fibre orientations in multiple fibre voxels are estimated using a mixture of von Mises-Fisher (vMF) distributions.
This statistical estimation procedure is used to resolve crossing fibre mappings.
}

\usage{
gqi.odfvmf(run=TRUE, fbase=NULL, rg=NULL, swap=FALSE, mddratio=1.2, depth=3, btoption=2, 
threshold=0.4, showglyph=FALSE, bview="coronal", savedir=tempdir())
}

\arguments{
\item{run}{ logical variable enabling loading previously processed data (default: \code{TRUE}). } 
\item{fbase}{ Directory where the required input data files are located. } 
\item{rg}{ range of slices to process; default option \code{rg=NULL} processes all slices. }
\item{swap}{ toggle radiological/neurological orientation  (default: \code{FALSE}). }
\item{mddratio}{ GQI parameter (default: 1.2). }
\item{depth}{ sampling density on the hemisphere used in simulation (default N=321; depth=3). }
\item{btoption}{ b-table selection between \file{dsi203_bmax4000.txt} (\code{btoption=1}), and the 3D-DSI grid b-table extracted from the diffusion data set (\file{data.bvec} and \file{data.bval}). By default, the 3D-DSI grid b-table is used (\code{btoption=2}). }
\item{threshold}{ thresholding generalized fractional anisotropy (GFA) value at each voxel (default: 0.4). }
\item{bview}{ MRI slice view selection in \{\code{axial, coronal, sagittal}\} (default: \code{"axial"}). }
\item{showglyph}{ logical variable controlling visualization of voxel glyphs (default: \code{FALSE}). }
\item{savedir}{ directory for saving/loading processed results (default: \code{tempdir()}. }
}

\details{
Generalized q-Sampling Imaging (GQI) is used for orientation distribution function (ODF) reconstruction.
GQI specifies an operational sampling scheme in q-space from which the ODF can be estimated.
For directional clustering estimation \code{gqi.odfvmf} uses a mixture of 2 and 4 von Mises-Fisher (vMF) distributions that serves as a model for directional ODF profile data, corresponding to multiple fibre orientations.
Statistical orientation estimation in \code{gqi.odfvmf} is based on von Mises clustering procedures provided by the R-package \pkg{movMF}, by Kurt Hornik and Bettina Gruen.

Starting with the raw diffusion signal acquired on a grid of q-space, the ODF profile is estimated at each voxel, considering a sampling density of unit vectors on a unit S2 grid.
When a threshold is applied to the estimated ODF at each voxel, the non-thresholded unit vectors provide directional statistics information about the estimated  ODF profile.
The main ODF orientations at each voxel relevant for fibre tracking may be estimated by clustering the non-thresholded unit vectors.

The main diffusion data set used in the examples is a DICOM data set provided by the "Advanced Biomedical MRI Lab, National Taiwan University Hospital",  which is included in the "DSI Studio" package,  publicly available from the NITRC repository (\url{http://www.nitrc.org}).
Two b-tables defining the acquisition setup are specified.
One is a b-table for a S2-like grid denoted by \file{dsi203_bmax4000.txt}. 
The other is the b-table for the 3D-DSI sampling scheme used in the DICOM data acquisition.
This b-table has 203 points uniformly distributed on a 3D grid limited to the volume of the unit sphere.
In both tables, the b-values range from 0 to 4000.
Sampling densities of N=81 (depth=2) and N=321 (depth=3) on the hemisphere are often used in ODF profile reconstruction from diffusion acquisitions. 

The output files \file{data_V1_gqi.nii.gz}, \file{data_V2_gqi.nii.gz} and \file{data_gfa_gqi.nii.gz} are used to visualize RGB maps through \code{rgbvolmap()} or via the "FSL/fslview" tool.
These files may be used to perform white matter fibre tractography. 
}

\value{
\code{gqi.odfvmf} outputs three data files in NIfTI format named \file{data_V1_gqi.nii.gz}, \file{data_V2_gqi.nii.gz}, and \file{data_gfa_gqi.nii.gz}. The first and second main fibre directions per voxel are contained in \file{data_V1_gqi.nii.gz}, \file{data_V2_gqi.nii.gz}, respectively. The file \file{data_gfa_gqi.nii.gz} contains the GFA metric per voxel.
}

\note{
Commonly used acronyms:\cr
GFA - Generalized Fractional Anisotropy\cr
ODF - Orientation Distribution Function\cr
GQI - Generalized q-Sampling Imaging\cr
MRI - Magnetic Resonance Imaging\cr
dMRI - Diffusion Magnetic Resonance Imaging\cr
RGB maps - Red:Green:Blue color maps\cr
HARDI - High Angular Resolution Diffusion Imaging\cr
vMF - von Mises-Fisher\cr
}

\references{
Ferreira da Silva, A. R. Facing the Challenge of Estimating Human Brain White Matter Pathways. In \emph{Proc. of the 4th International Joint Conference on Computational Intelligence} (Oct. 2012), K. Madani, J. Kacprzyk, and J. Filipe, Eds., SciTePress, pp. 709-714.

Hornik, K., and Gruen, B. \pkg{movMF}: Mixtures of von Mises-Fisher Distributions, 2012. R package version 0.1-0.

Yeh, F.-C., Wedeen, V. J., and Tseng, W.-Y. I. Generalized q-Sampling Imaging. \emph{IEEE Transactions on Medical Imaging} 29, 9 (2010), 1626-1635.

Hankin, R. K. S. Special functions in R: introducing the \pkg{gsl} package. R News 6 (October 2006).

Adler, D., and Murdoch, D. \pkg{rgl}: 3D visualization device system (OpenGL), 2012. R package version 0.92.880.

Auguie, B. \pkg{gridExtra}: functions in Grid graphics, 2012. R package version 0.9.1.

Barber, C. B., Habel, K., Grasman, R., Gramacy, R. B., Stahel, A., and Sterratt, D. C. \pkg{geometry}: Mesh generation and surface tesselation, 2012. R package version 0.3-2.

R Core Team. \pkg{R}: A Language and Environment for Statistical Computing.  R Foundation for Statistical Computing, Vienna, Austria, 2012. ISBN 3-900051-07-0.

Whitcher, B., Schmid, V. J., and Thornton, A. Working with the DICOM and NIfTI data standards in R. \emph{Journal of Statistical Software} 44, 6 (2011), 1-28.

Jenkinson, M., Beckmann, C. F., Behrens, T. E., Woolrich, M. W., and Smith, S. M. FSL. \emph{NeuroImage} 62, 2 (2012), 782-790.
}


\author{ Adelino Ferreira da Silva, Universidade Nova de Lisboa,
  Faculdade de Ciencias e Tecnologia, Portugal, \email{afs at fct.unl.pt}
}

\seealso{
  \code{\link{gqi.odfvmflines}},
  \code{\link{gqi.odfpeaklines}},
  \code{\link{gqi.odfvxgrid}},
  \code{\link{rgbvolmap}},
  \code{\link{gqi.odfpeaks}},
  \code{\link{s2tessel.zorder}},
  \code{\link{plotglyph}},
  \code{\link{simulglyph.vmf}},
  \code{\link{simul.fandtasia}},
  \code{\link{simul.simplefield}},
  \code{\link{data}},
  \code{\link{data_brain}},
  \code{\link{data.bval}},
  \code{\link{data.bvec}},
  \code{\link{dsi203_bmax4000.Rd}}
}

\examples{
\dontrun{
    ## Generate ODF volumes (GQI volume processing)
    ## for a range of slices using von Mises clustering
    gqi.odfvmf(depth=2, showglyph=FALSE, threshold=0.5, savedir=tempdir())
    ## RGB maps for range of slices processed by gqi.odfvmf()
    rgbvolmap(fbase=tempdir(), rg=c(1,4), bview="coronal")
    ##-------------
    ## Show reconstructed glyphs in ODF processing 
    ## for first and second main fibre direction determination
    gqi.odfvmf(rg=c(1,1), bview="coronal", depth=3, showglyph=TRUE,
      threshold=0.5)
}
}

\keyword{ methods }
\keyword{ clustering }


