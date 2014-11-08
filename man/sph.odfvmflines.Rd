\name{sph.odfvmflines}
\alias{sph.odfvmflines}
\concept{Diffusion Magnetic Resonance}
\concept{QBI Reconstruction}
\concept{von Mises Distributions}
\concept{Fibre Orientation Mapping}

\title{Fibre Orientation Mapping Based on von Mises-Fisher Distributions with QBI reconstruction}
\description{
In order to enable mapping complex white matter fibres in the brain, \code{sph.odfvmflines} implements a new methodology based on directional statistics to estimate fibre profiles from high angular resolution diffusion imaging data.
Statistical orientation estimation in \code{sph.odfvmf} is based on von Mises-Fisher clustering procedures provided by the R-package \pkg{movMF}, by Kurt Hornik and Bettina Gruen.
}

\usage{
sph.odfvmflines(fbase=NULL, roi=NULL,  rg=c(1,1), swap=FALSE,
 btoption=2, threshold=0.4, kdir=4, zfactor=5, showglyph=FALSE, showimage="linesgfa",
 bview="coronal", bg="white", order=4, texture=NULL, clusterthr=0.6, aniso=NULL, ...)
}

\arguments{
\item{fbase}{ Directory where the required input data files are located. } 
\item{roi}{ Region of Interest (ROI) to use as mask; default mask (\code{roi=NULL}) uses an all brain mask for the supplied data set. }
\item{rg}{ range of slices to process (default option \code{rg=c(1,1)}); \code{rg=NULL} processes all slices. }
\item{swap}{ toggle radiological/neurological orientation  (default: \code{FALSE}). }
\item{btoption}{ b-table selection between \file{btable.txt} (\code{btoption=1}), and the 3D-DSI grid b-table extracted from the diffusion data set (\file{data.bvec} and \file{data.bval}). By default, the 3D-DSI grid b-table is used (\code{btoption=2}). }
\item{threshold}{ thresholding generalized fractional anisotropy (GFA) value at each voxel (default: 0.4). }
\item{kdir}{ maximum number of fibre directions to map (default: 4). }
\item{zfactor}{ parameter controlling z-value in relief overlay maps (default: 5). }
\item{showglyph}{ logical variable controlling visualization of voxel glyphs (default: \code{FALSE}). }
\item{showimage}{ object controlling visualization of line-maps (default: \code{"linesgfa"}). \cr
  Alternative options are:\cr
  \code{c("none", "gfa", "lines", "linesgfa", "linesrgbmap", "linesdata")} \cr
 (see Details). }
\item{bview}{ MRI slice view selection in \{\code{axial, coronal, sagittal}\} (default: \code{"coronal"}). }
\item{bg}{ map background colour (default \code{"white"}) }
\item{order}{ parameter associated with the order of the spherical harmonics approximation (default: 4). }
\item{texture}{ name of the PNG file to be used as RGB map overlay in some 'showimage' options (default \code{NULL} - no texture). }
\item{clusterthr}{ thresholding orientations based on ODF values at each voxel for directional clustering (default: 0.6). }
\item{aniso}{ anisotropic parameter in the range "[0,1)" or \code{NULL} to use in ODF pos-processing default: \code{NULL}. }
\item{...}{ additional material properties for geometry appearance as specified in \code{rgl.material}. }
}

\details{
The function \code{sph.odfvmflines} implements a mixture-model approach to clustering orientation distribution functions (ODFs) based on the von Mises-Fisher distributions.
The method focus on clustering data on the unit sphere, where complexity arises from representing ODF profiles as directional data.

Starting with the raw high angular resolution diffusion signal acquired on a S2-shell of q-space, the ODF profile is reconstructed at each voxel, considering a sampling density of unit vectors on a unit S2 shell.
Q-ball imaging (QBI) is used for orientation distribution function (ODF) reconstruction.
For comparison with GQI, the b-table \file{byable.txt} has been used in the examples.  This b-table has 203 points distributed on a S2-shell.

Slice map display and overlay selection is controlled by specifying one the arguments \cr
 \code{c("none", "gfa", "lines", "linesgfa", "linesrgbmap", "linesdata")} \cr
for \code{showimages}.
Meanings are as follows: \code{"none"} - no visualization; \code{"gfa"} - GFA map only; \code{"lines"} - line map only; \code{"linesgfa"} - GFA overlayed on line map; \code{"linesrgbmap"} - lines overlayed on RGB map (if available); \code{"linesdata"} - \file{data_brain.nii.gz} is overlayed on line map.
}

\value{
\code{sph.odfvmflines} produces line-maps of ODF profiles for diffusion data slices.
The line-maps may be overlayed with generalized fractional anisotropy (GFA) relief maps, diffusion data maps or ROI maps.
The file \file{V1list.RData} containing the first main orientation directions for all processed voxels is output for further posterior orientation processing.
}


\references{
Ferreira da Silva, A. R. Computational Representation of White Matter Fiber Orientations, \emph{International Journal of Biomedical Imaging}, Vol. 2013, Article ID 232143, Hindawi Publishing Corporation \url{http://dx.doi.org/10.1155/2013/232143}.

Ferreira da Silva, A. R. Facing the Challenge of Estimating Human Brain White Matter Pathways. In \emph{Proc. of the 4th International Joint Conference on Computational Intelligence} (Oct. 2012), K. Madani, J. Kacprzyk, and J. Filipe, Eds., SciTePress, pp. 709-714.

Hornik, K., and Gruen, B. \pkg{movMF}: Mixtures of von Mises-Fisher Distributions, 2012. R package version 0.1-0.

Tuch D. S., Q-Ball Imaging, \emph{Magnetic Resonance in Medicine} 52 (2004), 1358-1372. 

Tabelow K., Polzehl J.: \pkg{dti}: DTI/DWI Analysis, 2012. R package version 1.1-0.
}


\author{ Adelino Ferreira da Silva, Universidade Nova de Lisboa,
  Faculdade de Ciencias e Tecnologia, Portugal, \email{afs at fct.unl.pt}
}

\seealso{
  \code{\link{sph.odfvmf}},
  \code{\link{sph.odfpeaks}},
  \code{\link{sph.odfvmflines}},
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
  sph.odfvmflines(rg=c(1,1), showimage="linesdata",
    threshold=0.5)
  ## display line-maps only
  sph.odfvmflines(showimage="lines")
  ## using GFA overlay
  sph.odfvmflines(showimage="linesgfa")
  ##-------------
  ## Show reconstructed glyphs in ODF processing 
  ## for principal direction determination
  sph.odfvmflines(showimage="linesdata",
    showglyph=TRUE, threshold=0.5)
  ##-------------
  ## using a ROI overlay
  sph.odfvmflines(roi="slfcst.nii.gz")
  ##-------------
  ## speeded up approximations: hardmax and common/numeric kappa
  sph.odfvmflines(rg=c(1,1), showimage="linesdata",
    E="hardmax", kappa=list(common = TRUE))
  sph.odfvmflines(rg=c(1,1), showimage="linesdata",
    E="hardmax", kappa=20)
}
}

\keyword{ methods }
\keyword{ clustering }
\keyword{ dplot }

