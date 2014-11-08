\name{gqi.odfvmflines}
\alias{gqi.odfvmflines}
\concept{Diffusion Magnetic Resonance}
\concept{GQI Reconstruction}
\concept{von Mises Distributions}
\concept{Fibre Orientation Mapping}

\title{Fibre Orientation Mapping Based on von Mises-Fisher Distributions with GQI Reconstruction}
\description{
In order to enable mapping complex white matter fibres in the brain, \code{gqi.odfvmflines}
 implements a new methodology based on directional statistics to estimate fibre profiles
 from high angular resolution diffusion imaging data.
Statistical orientation estimation in \code{gqi.odfvmf} and \code{gqi.odfvmflines} is based
 on von Mises-Fisher clustering procedures provided by the R-package \pkg{movMF},
 by Kurt Hornik and Bettina Gruen.
}

\usage{
gqi.odfvmflines(gdi="gqi", fbase=NULL, roi=NULL,  rg=c(1,1), swap=FALSE,
 lambda=NULL, depth=3, btoption=2, threshold=0.4, kdir=6, zfactor=5,
 showglyph=FALSE, showimage="linesgfa", bview="coronal", bg="white",
 texture=NULL, clusterthr=0.6, aniso=NULL, ...)
}

\arguments{
\item{gdi}{ method of ODF reconstruction to use \code{c("gqi", "gqi2")} (default: \code{"gqi"}). }
\item{fbase}{ Directory where the required input data files are located. } 
\item{roi}{ Region of Interest (ROI) to use as mask; default mask (\code{roi=NULL}) uses an all brain mask for the supplied data set. }
\item{rg}{ range of slices to process (default option \code{rg=c(1,1)}); \code{rg=NULL} processes all slices. }
\item{swap}{ toggle radiological/neurological orientation  (default: \code{FALSE}). }
\item{lambda}{ diffusion sampling length in \code{gdi="gqi"} and \code{gdi="gqi2"}.
 By default the following default values are used when \code{lambda=NULL} is specified: 1.25 in "gqi", 3 in "gqi2". }
\item{depth}{ sampling density on the hemisphere used in simulation (default N=321; depth=3). }
\item{btoption}{ b-table selection between \file{btable.txt} (\code{btoption=1}), and the 3D-DSI grid b-table extracted from the diffusion data set (\file{data.bvec} and \file{data.bval}). By default, the 3D-DSI grid b-table is used (\code{btoption=2}). }
\item{threshold}{ thresholding generalized fractional anisotropy (GFA) value at each voxel (default: 0.4). }
\item{kdir}{ maximum number of fibre directions to map (default: 6). }
\item{zfactor}{ parameter controlling z-value in relief overlay maps (default: 5). }
\item{showglyph}{ logical variable controlling visualization of voxel glyphs (default: \code{FALSE}). }
\item{showimage}{ object controlling visualization of line-maps (default: \code{"linesgfa"}). \cr
  Alternative options are:\cr
  \code{c("none", "gfa", "lines", "linesgfa", "linesrgbmap", "linesdata")} \cr
 (see Details). }
\item{bview}{ MRI slice view selection in \{\code{axial, coronal, sagittal}\} (default: \code{"coronal"}). }
\item{bg}{ map background colour (default \code{"white"}). }
\item{texture}{ name of the PNG file to be used as RGB map overlay in some 'showimage' options (default \code{NULL} - no texture). }
\item{clusterthr}{ thresholding orientations based on ODF values at each voxel for directional clustering (default: 0.6). }
\item{aniso}{ anisotropic parameter in the range "[0,1)" or \code{NULL} to use in ODF pos-processing default: \code{NULL}. }
\item{...}{ additional material properties for geometry appearance as specified in \code{rgl.material}, or specification of non-default control parameters as detailed in \code{movMF}. }
}

\details{
The function \code{gqi.odfvmflines} implements a mixture-model approach to clustering orientation distribution functions (ODFs) based on the von Mises-Fisher distributions.
The method focus on clustering data on the unit sphere, where complexity arises from representing ODF profiles as directional data.
GQI (Yeh et.al. 2010) or GQI2 (Garyfallidis 2012) may be used for ODF reconstruction.

Starting with the raw diffusion signal acquired on a grid of q-space, the ODF profile is estimated at each voxel, considering a sampling density of unit vectors on a unit S2 grid.
When a threshold is applied to the estimated ODF at each voxel, the non-thresholded unit vectors provide directional statistics information about the estimated  ODF profile.
The main ODF orientations at each voxel relevant for fibre tracking may be estimated by clustering the non-thresholded unit vectors.

The main diffusion data set used in the examples is a DICOM data set provided by the "Advanced Biomedical MRI Lab, National Taiwan University Hospital",  which is included in the "DSI Studio" package,  publicly available from the NITRC repository (\url{http://www.nitrc.org}).
Two b-tables defining the acquisition setup are specified.
One is a b-table for a S2-like grid denoted by \file{btable.txt}. 
The other is the b-table for the 3D-DSI sampling scheme used in the DICOM data acquisition.
This b-table has 203 points uniformly distributed on a 3D grid limited to the volume of the unit sphere.
In both tables, the b-values range from 0 to 4000.

Slice map display and overlay selection is controlled by specifying one the arguments \cr
 \code{c("none", "gfa", "lines", "linesgfa", "linesrgbmap", "linesdata")} \cr
for \code{showimages}.
Meanings are as follows: \code{"none"} - no visualization; \code{"gfa"} - GFA map only; \code{"lines"} - line map only; \code{"linesgfa"} - GFA overlayed on line map; \code{"linesrgbmap"} - lines overlayed on RGB map (if available); \code{"linesdata"} - \file{data_brain.nii.gz} is overlayed on line map.
}

\value{
\code{gqi.odfvmflines} produces line-maps of ODF profiles for diffusion data slices.
The line-maps may be overlayed with generalized fractional anisotropy (GFA) relief maps, diffusion data maps or ROI maps.
The file \file{V1list.RData} containing the first main orientation directions for all processed voxels is output for further posterior orientation processing.
}

\references{
Ferreira da Silva, A. R., Modeling White-matter Fiber-orientation Uncertainty for Improved Probabilistic Tractography. In \emph{Proc. of the  2nd Int. Congr. on Neurotechnology, Electronics and Informatics, NEUROTECHNIX 2014}, Rome, Italy, (Oct., 2014), A. Londral and P. Encarnacao, Eds., INSTICC, SCITEPRESS, pp 71-78.

Ferreira da Silva, A. R. Computational Representation of White Matter Fiber Orientations, \emph{International Journal of Biomedical Imaging}, Vol. 2013, Article ID 232143, Hindawi Publishing Corporation \url{http://dx.doi.org/10.1155/2013/232143}.

Ferreira da Silva, A. R. Facing the Challenge of Estimating Human Brain White Matter Pathways. In \emph{Proc. of the 4th International Joint Conference on Computational Intelligence} (Oct. 2012), K. Madani, J. Kacprzyk, and J. Filipe, Eds., SciTePress, pp. 709-714.

Hornik, K., and Gruen, B. \pkg{movMF}: Mixtures of von Mises-Fisher Distributions, 2014. R package version 0.2-0.

Yeh, F.-C., Wedeen, V. J., and Tseng, W.-Y. I. Generalized q-Sampling Imaging. \emph{IEEE Transactions on Medical Imaging} 29, 9 (2010), 1626-1635.

Garyfallidis E., Towards an Accurate Brain Tractography, 2012, PhD Thesis, University of Cambridge.
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
  \code{\link{data.bval}},
  \code{\link{data.bvec}},
  \code{\link{btable}}

}

\examples{
\dontrun{
  ##-------------
  ## von Mises-Fisher fibre orientation mapping
  ## for a range of slices
  gqi.odfvmflines(gdi="gqi", rg=c(1,1), depth=2,
    showimage="linesdata", threshold=0.5)
  gqi.odfvmflines(gdi="gqi2", rg=c(1,1), depth=2,
    showimage="linesdata", threshold=0.5)
  ## display line-maps only
  gqi.odfvmflines(depth=2, showimage="lines")
  ## using GFA overlay
  gqi.odfvmflines(depth=2, showimage="linesgfa")
  ##-------------
  ## Show reconstructed glyphs in ODF processing 
  ## for principal direction determination
  gqi.odfvmflines(depth=3, showimage="linesdata",
    showglyph=TRUE, threshold=0.5)
  ## show glyphs with using 'aniso' parameter 
  gqi.odfvmflines(depth=3, showimage="linesdata",
    showglyph=TRUE, threshold=0.5, aniso=0.3)
  ##-------------
  ## using a ROI overlay
  gqi.odfvmflines(depth=3, roi="slfcst.nii.gz")
  ##-------------
  ## coronal view with texture for a single slice
  texturefname <- file.path(tempdir(),"rgbmap.png")
  rgbvolmap(texture=texturefname, bg="transparent")
  gqi.odfvmflines(threshold=0.5, showimage="linesrgbmap",
    texture=texturefname)
  ##-------------
  ## speeded up approximations: hardmax and common/numeric kappa
  gqi.odfvmflines(gdi="gqi", rg=c(1,1), depth=2,
    showimage="linesdata", threshold=0.5,
    E="hardmax", kappa=list(common = TRUE))
  gqi.odfvmflines(gdi="gqi", rg=c(1,1), depth=2,
    showimage="linesdata", threshold=0.5, E="hardmax", kappa=20)
}
}

\keyword{ methods }
\keyword{ clustering }
\keyword{ dplot }

