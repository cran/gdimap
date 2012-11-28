\name{gqi.odfvmflines}
\alias{gqi.odfvmflines}
\concept{Diffusion Magnetic Resonance}
\concept{GQI Reconstruction}
\concept{von Mises Distributions}
\concept{Fibre Orientation Mapping}

\title{Fibre Orientation Mapping Based on von Mises-Fisher Distributions}
\description{
In order to enable mapping complex white matter fibres in the brain, \code{gqi.odfvmflines} implements a new methodology based on directional statistics to estimate fibre profiles from high angular resolution diffusion imaging data.
Statistical orientation estimation in \code{gqi.odfvmf} is based on von Mises-Fisher clustering procedures provided by the R-package \pkg{movMF}, by Kurt Hornik and Bettina Gruen.
}

\usage{
gqi.odfvmflines(run=TRUE, fbase=NULL, roi=NULL,  rg=c(1,1), swap=FALSE, mddratio=1.2,
 depth=3, btoption=2, threshold=0.4, kdir=4, zfactor=5, showglyph=FALSE, snapshot=FALSE,
 showimage="linesgfa", bview="coronal", savedir=tempdir(), pngfig="odfvmf", bg="white",
 texturefile=NA)
}

\arguments{
\item{run}{ logical variable enabling loading previously processed data (default: \code{TRUE}). } 
\item{fbase}{ Directory where the required input data files are located. } 
\item{roi}{ Region of Interest (ROI) to use as mask; default mask (\code{roi=NULL}) uses an all brain mask for the supplied data set. }
\item{rg}{ range of slices to process (default option \code{rg=c(23,23)}); \code{rg=NULL} processes all slices. }
\item{swap}{ toggle radiological/neurological orientation  (default: \code{FALSE}). }
\item{mddratio}{ GQI parameter (default: 1.2). }
\item{depth}{ sampling density on the hemisphere used in simulation (default N=321; depth=3). }
\item{btoption}{ b-table selection between \file{dsi203_bmax4000.txt} (\code{btoption=1}), and the 3D-DSI grid b-table extracted from the diffusion data set (\file{data.bvec} and \file{data.bval}). By default, the 3D-DSI grid b-table is used (\code{btoption=2}). }
\item{threshold}{ thresholding generalized fractional anisotropy (GFA) value at each voxel (default: 0.4). }
\item{kdir}{ maximum number of fibre directions to map (default: 4). }
\item{zfactor}{ parameter controlling z-value in relief overlay maps (default: 5). }
\item{showglyph}{ logical variable controlling visualization of voxel glyphs (default: \code{FALSE}). }
\item{snapshot}{ logical variable controlling line-maps snapshots (default: \code{FALSE}). }
\item{showimage}{ object controlling visualization of line-maps (default: \code{"linesgfa"}). Alternative options are \code{c("none", "gfa", "lines", "linesgfa", "linesrgbmap", "linesdata")} (see Details). }
\item{bview}{ MRI slice view selection in \{\code{axial, coronal, sagittal}\} (default: \code{"axial"}). }
\item{savedir}{ directory for saving/loading processed results (default: \code{tempdir()}. }
\item{pngfig}{ base-name of the PNG file when \code{snapshot} is \code{TRUE} (default \code{"odfvmf"}). }
\item{bg}{ map background colour (default \code{"white"}) }
\item{texturefile}{ name of the PNG file to be used as RGB map overlay in some 'showimage' options (default \code{NA} - no texture). }
}

\details{
Generalized q-Sampling Imaging (GQI) is used for orientation distribution function (ODF) reconstruction.
The function \code{gqi.odfvmflines} implements a mixture-model approach to clustering orientation distribution functions (ODFs) based on the von Mises-Fisher distributions.
The method focus on clustering data on the unit sphere, where complexity arises from representing ODF profiles as directional data.

Starting with the raw diffusion signal acquired on a grid of q-space, the ODF profile is estimated at each voxel, considering a sampling density of unit vectors on a unit S2 grid.
When a threshold is applied to the estimated ODF at each voxel, the non-thresholded unit vectors provide directional statistics information about the estimated  ODF profile.
The main ODF orientations at each voxel relevant for fibre tracking may be estimated by clustering the non-thresholded unit vectors.

The main diffusion data set used in the examples is a DICOM data set provided by the "Advanced Biomedical MRI Lab, National Taiwan University Hospital",  which is included in the "DSI Studio" package,  publicly available from the NITRC repository (\url{http://www.nitrc.org}).
Two b-tables defining the acquisition setup are specified.
One is a b-table for a S2-like grid denoted by \file{dsi203\_bmax4000.txt}. 
The other is the b-table for the 3D-DSI sampling scheme used in the DICOM data acquisition.
This b-table has 203 points uniformly distributed on a 3D grid limited to the volume of the unit sphere.
In both tables, the b-values range from 0 to 4000.

Slice map display and overlay selection is controlled by specifying one the arguments \code{c("none", "gfa", "lines", "linesgfa", "linesrgbmap", "linesdata")} for \code{showimages}. Meanings are as follows: \code{"none"} - no visualization; \code{"gfa"} - GFA map only; \code{"lines"} - line map only; \code{"linesgfa"} - GFA overlayed on line map; \code{"linesrgbmap"} - lines overlayed on RGB map (if available); \code{"linesdata"} - \file{data_brain.nii.gz} is overlayed on line map.
}

\value{
\code{gqi.odfvmflines} produces line-maps of ODF profiles for diffusion data slices.
The line-maps may be overlayed with generalized fractional anisotropy (GFA) relief maps, diffusion data maps or ROI maps.
The file \file{V1list.RData} containing the first main orientation directions for all processed voxels is output for further posterior orientation processing.
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

Zeileis, A., Hornik, K., and Murrell, P. Escaping RGBland: Selecting colors for statistical graphics. Computational Statistics & Data Analysis 53 (2009), 3259-3270.
}


\author{ Adelino Ferreira da Silva, Universidade Nova de Lisboa,
  Faculdade de Ciencias e Tecnologia, Portugal, \email{afs at fct.unl.pt}
}

\seealso{
  \code{\link{gqi.odfvmf}},
  \code{\link{gqi.odfpeaks}},
  \code{\link{gqi.odfvmflines}},
  \code{\link{gqi.odfvxgrid}},
  \code{\link{plotglyph}},
  \code{\link{rgbvolmap}},
  \code{\link{s2tessel.zorder}},
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
  ##-------------
  ## von Mises Fisher fibre orientation mapping
  ## for a range of slices
  gqi.odfvmflines(run=TRUE, rg=c(1,1), depth=2,
    showimage="linesdata", threshold=0.5)
  ## display line-maps only
  gqi.odfvmflines(run=FALSE, depth=2, showimage="lines")
  ## using GFA overlay
  gqi.odfvmflines(run=FALSE, depth=2, showimage="linesgfa")
  ##-------------
  ## Show reconstructed glyphs in ODF processing 
  ## for principal direction determination
  gqi.odfvmflines(run=TRUE, depth=3,
    showimage="linesdata", showglyph=TRUE, threshold=0.5)
  ##-------------
  ## using a ROI overlay
  gqi.odfvmflines(run=TRUE, depth=3, roi="slfcst")
  ##-------------
  ## coronal view with texture for a single slice
  rgbvolmap(texture=TRUE, transparent=TRUE)
  gqi.odfvmflines(threshold=0.5, showimage="linesrgbmap",
    texturefile=paste(tempdir(),"/rgbmap.png", sep=""))
}
}

\keyword{ methods }
\keyword{ clustering }
\keyword{ dplot }

