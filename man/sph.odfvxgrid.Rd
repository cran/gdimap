\name{sph.odfvxgrid}
\alias{sph.odfvxgrid}
\concept{GQI Volume Processing}
\concept{Glyph maps}
\title{Glyph Maps}

\description{
\code{sph.odfvxgrid} produces glyph-map fields of reconstructed ODF profiles for voxels in slices.
}

\usage{
sph.odfvxgrid(fbase=NULL, rg=c(1,1), swap=FALSE, btoption=1, threshold=0.4,
 kdir=4, zfactor=5, snapshot=FALSE, showimage="glyphgfa", bview="coronal",
 savedir=tempdir(), pngfig="odfglyphs", bg="white", order=4, texture=NULL, ...)
}

\arguments{
\item{fbase}{ Directory where the required input data files are located. } 
\item{rg}{ range of slices to process (default option \code{rg=c(1,1)}); \code{rg=NULL} processes all slices. }
\item{swap}{ toggle radiological/neurological orientation  (default: \code{FALSE}). }
\item{btoption}{ b-table selection using the included data file \file{btable.txt} (\code{btoption=1})
 or user defined files. With \code{btoption=2}, diffusion data files \file{data.bvec} and \file{data.bval}
 are selected for input, as used in single-shell QBI acquisition. In this case, the files must be provided by the user.
 (Default \code{btoption=1}). }
\item{threshold}{ thresholding generalized fractional anisotropy (GFA) value at each voxel (default: 0.4). }
\item{kdir}{ maximum number of fibre directions to map (default: 4). }
\item{zfactor}{ parameter controlling z-value in relief overlay maps (default: 5). }
\item{snapshot}{ logical variable controlling line-maps snapshots (default: \code{FALSE}). }
\item{showimage}{ object controlling visualization of line-maps (default: \code{"glyphgfa"}). \cr
  Alternative options are: \cr
 \code{c("none", "gfa", "glyph", "glyphgfa", "glyphrgbmap", "glyphdata")} \cr 
 (see Details). }
\item{bview}{ MRI slice view selection in \{\code{axial, coronal, sagittal}\} (default: \code{"coronal"}). }
\item{savedir}{ directory for saving/loading processed results (default: \code{tempdir()}. }
\item{pngfig}{ base-name of the PNG file when \code{snapshot} is \code{TRUE} (default \code{"odfglyphs"}). }
\item{bg}{ map background colour (default \code{"white"}) }
\item{order}{ parameter associated with the order of the spherical harmonics approximation (default: 4). }
\item{texture}{ name of the PNG file to be used as RGB map overlay in some 'showimage' options (default \code{NULL} - no texture). }
\item{...}{ additional material properties for geometry appearance as specified in \code{rgl.material}. } 
}

\details{
Slice glyph map display and overlay selection is controlled by specifying one the arguments \cr
\code{c("none", "gfa", "glyph", "glyphgfa", "glyphrgbmap", "glyphdata")} \cr
 for \code{showimages}.
Meanings are as follows: \code{"none"} - no overlay; \code{"gfa"} - GFA map only; \code{"glyph"} - glyph map only; \code{"glyphgfa"} - GFA overlayed on glyph map; \code{"glyphrgbmap"} - glyphs overlayed on RGB map (if available); \code{"glyphdata"} - \file{data_brain.nii.gz} is overlayed on line map.

Q-ball imaging (QBI) is used for orientation distribution function (ODF) reconstruction.
One QBI-compatible b-table defining the acquisition setup on a S2-shell is included in the package. The b-table \file{btable.txt} has 203 points distributed on a S2-shell.
}

\references{
Ferreira da Silva, A. R. Computational Representation of White Matter Fiber Orientations, \emph{International Journal of Biomedical Imaging}, Vol. 2013, Article ID 232143, Hindawi Publishing Corporation \url{http://dx.doi.org/10.1155/2013/232143}.

Ferreira da Silva, A. R. Facing the Challenge of Estimating Human Brain White Matter Pathways. In \emph{Proc. of the 4th International Joint Conference on Computational Intelligence} (Oct. 2012), K. Madani, J. Kacprzyk, and J. Filipe, Eds., SciTePress, pp. 709-714.

Tuch D. S., Q-Ball Imaging, \emph{Magnetic Resonance in Medicine} 52 (2004), 1358-1372. 
Tabelow K., Polzehl J.: \pkg{dti}: DTI/DWI Analysis, 2012. R package version 1.1-0.
} 

\author{ Adelino Ferreira da Silva, Universidade Nova de Lisboa,
  Faculdade de Ciencias e Tecnologia, Portugal, \email{afs at fct.unl.pt}
}

\seealso{
  \code{\link{sph.odfpeaks}},
  \code{\link{sph.odfvmf}},
  \code{\link{sph.odfvmflines}},
  \code{\link{s2tessel.zorder}},
  \code{\link{plotglyph}},
  \code{\link{rgbvolmap}},
  \code{\link{simulglyph.vmf}},
  \code{\link{simul.fandtasia}},
  \code{\link{simul.simplefield}}
}

\examples{
\dontrun{
  ## show glyph-map for selected slice
  sph.odfvxgrid(threshold=0.5, showimage="glyphgfa")
  ##----
  ## generate slice RGB map as texture
  rgbvolmap(texture=TRUE, transp=TRUE)
  ## coronal glyphs + RGBmap
  sph.odfvxgrid(threshold=0.5, showimage="glyphrgbmap",
    texture=file.path(tempdir(),"rgbmap.png"))
}
}

\keyword{ methods }
\keyword{ dplot }

