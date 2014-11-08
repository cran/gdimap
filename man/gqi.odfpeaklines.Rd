\name{gqi.odfpeaklines}
\alias{gqi.odfpeaklines}
\concept{GQI Reconstruction}
\concept{ODF Line-maps}
\concept{Peak Detection}

\title{Fibre Orientation Mapping Based on Local Peak Detection}
\description{
\code{gqi.odfpeaklines} produces line-maps of ODF profiles for diffusion data slices.
}

\usage{
gqi.odfpeaklines(gdi="gqi", fbase=NULL, roi=NULL, rg=c(1,1),
 swap=FALSE, lambda=NULL, depth=3, btoption=2, threshold=0.4, kdir=2,
 zfactor=5, showglyph=FALSE, showimage="linesgfa", bview="coronal",
 savedir=tempdir(), bg="white", texture=NULL, aniso=NULL, ...)
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
\item{btoption}{ b-table selection between \file{btable.txt} (\code{btoption=1}), and the 3D-DSI grid b-table
  extracted from the diffusion data set (\file{data.bvec} and \file{data.bval}). By default, the 3D-DSI grid b-table is used (\code{btoption=2}). }
\item{threshold}{ thresholding generalized fractional anisotropy (GFA) value at each voxel (default: 0.4). }
\item{kdir}{ maximum number of fibre directions to map (default: 2). }
\item{zfactor}{ parameter controlling z-value in relief overlay maps (default: 5). }
\item{showglyph}{ logical variable controlling visualization of voxel glyphs (default: \code{FALSE}). }
\item{showimage}{ object controlling visualization of line-maps (default: \code{"linesgfa"}). \cr
  Alternative options are: \cr
  \code{c("none", "gfa", "lines", "linesgfa", "linesrgbmap", "linesdata")} \cr
  (see Details). }
\item{bview}{ MRI slice view selection in \{\code{axial, coronal, sagittal}\} (default: \code{"coronal"}). }
\item{savedir}{ directory for saving/loading processed results (default: \code{tempdir()}. }
\item{bg}{ map background colour (default \code{"white"}) }
\item{texture}{ name of the PNG file to be used as RGB map overlay in some 'showimage' options (default \code{NULL} - no texture). }
\item{aniso}{ anisotropic parameter in the range "[0,1)" or \code{NULL} to use in ODF pos-processing default: \code{NULL}. }
\item{...}{ additional material properties for geometry appearance as specified in \code{rgl.material}. } 
}

\details{
The identification of fibre directions is performed by extracting the local maxima of the reconstructed ODF, where this function surpasses a certain threshold.
GQI (Yeh et.al. 2010) or GQI2 (Garyfallidis 2012) may be used for ODF reconstruction.
\code{gqi.odfpeaklines} implements the standard method of fibre orientation detection.
Local maxima of the reconstructed ODF are located simply by selecting a large number of sampled points on the sphere and searching within a fixed radius neighbourhood. 
For a single main fibre orientation the method performs well.
However, for crossing fibres and other complex fibre configurations the peaks of the ODF profiles identified by the methods do not necessarily match the orientations of the distinct fibre populations.
A more robust method is implemented in \code{gqi.odfvmflines}.

Starting with the raw high angular resolution diffusion signal acquired on a grid of q-space, the ODF profile is reconstructed at each voxel, considering a sampling density of unit vectors on a unit S2 grid.
Generalized q-Sampling Imaging (GQI) is used for orientation distribution function (ODF) reconstruction.
Two b-tables defining the acquisition setup are specified.
One is a b-table for a S2-like grid denoted by \file{btable.txt}. 
The other is the b-table for the 3D-DSI sampling scheme used in the DICOM data acquisition.
This b-table has 203 points uniformly distributed on a 3D grid limited to the volume of the unit sphere.
In both tables, the b-values range from 0 to 4000.

Slice map display and overlay selection is controlled by specifying one the arguments\cr
 \code{c("none", "gfa", "lines", "linesgfa", "linesrgbmap", "linesdata")} for \code{showimages}.
 Meanings are as follows:
 \code{"none"} - no visualization;
 \code{"gfa"} - GFA map only;
 \code{"lines"} - line map only;
 \code{"linesgfa"} - GFA overlayed on line map;
 \code{"linesrgbmap"} - lines overlayed on RGB map (if available);
 \code{"linesdata"} - \file{data_brain.nii.gz} is overlayed on line map;
}

\value{
\code{gqi.odfpeaklines} produces line-maps of ODF profiles for diffusion data slices.
The line-maps may be overlayed with generalized fractional anisotropy (GFA) relief maps, diffusion data maps or ROI maps.
The file \file{V1list.RData} containing the first main orientation directions for all processed voxels is output for further posterior orientation processing.
}

\references{
Ferreira da Silva, A. R. Computational Representation of White Matter Fiber Orientations, \emph{International Journal of Biomedical Imaging}, Vol. 2013, Article ID 232143, Hindawi Publishing Corporation \url{http://dx.doi.org/10.1155/2013/232143}.

Ferreira da Silva, A. R. Facing the Challenge of Estimating Human Brain White Matter Pathways. In \emph{Proc. of the 4th International Joint Conference on Computational Intelligence} (Oct. 2012), K. Madani, J. Kacprzyk, and J. Filipe, Eds., SciTePress, pp. 709-714.

Yeh, F.-C., Wedeen, V. J., and Tseng, W.-Y. I. Generalized q-Sampling Imaging. \emph{IEEE Transactions on Medical Imaging} 29, 9 (2010), 1626-1635.

Garyfallidis E., Towards an Accurate Brain Tractography, 2012, PhD Thesis, University of Cambridge.
}

\author{ Adelino Ferreira da Silva, Universidade Nova de Lisboa,
  Faculdade de Ciencias e Tecnologia, Portugal, \email{afs at fct.unl.pt}
}

\seealso{
  \code{\link{gqi.odfpeaks}},
  \code{\link{gqi.odfvmf}},
  \code{\link{gqi.odfvmflines}},
  \code{\link{gqi.odfvxgrid}},
  \code{\link{s2tessel.zorder}},
  \code{\link{plotglyph}},
  \code{\link{rgbvolmap}},
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
  ## Line map using ODF peak detection 
  gqi.odfpeaklines(showimage="lines")
  ## idem with "gqi2"
  gqi.odfpeaklines(gdi="gqi2", showimage="lines")
  ## display line-map overlayed on GFA map 
  gqi.odfpeaklines(showimage="linesgfa")
  ##-------------
  ## generate slice texture first from default data file
  texturefname <- file.path(tempdir(),"rgbmap.png")
  rgbvolmap(texture=texturefname, bg="transparent")
  ## Line map with RGB map overlay
  gqi.odfpeaklines(showimage="linesrgbmap",
    texture=texturefname)
  ##-------------
  ## Show examples of reconstructed glyphs in ODF processing 
  gqi.odfpeaklines(showimage="lines", showglyph=TRUE)
  ##------------
  ## using a ROI overlay
  gqi.odfpeaklines(roi="slfcst.nii.gz", showimage="linesgfa")
  ## using data overlay
  gqi.odfpeaklines(showimage="linesdata")
}
}

\keyword{ methods }
\keyword{ dplot }
\keyword{ utilities }

