\name{gqi.odfvxgrid}
\alias{gqi.odfvxgrid}
\concept{GQI Volume Processing}
\concept{Glyph maps}
\title{Glyph Maps}

\description{
\code{gqi.odfvxgrid} produces glyph-map fields of reconstructed ODF profiles for voxels in slices.
}

\usage{
gqi.odfvxgrid(gdi="gqi", fbase=NULL, rg=c(1,1), swap=FALSE, lambda=NULL,
 depth=3, btoption=2, threshold=0.4, kdir=4, zfactor=5,
 snapshot=FALSE, showimage="glyphgfa", bview="coronal",
 savedir=tempdir(), pngfig="odfglyphs", bg="white", texture=NULL, ...)
}

\arguments{
\item{gdi}{ method of ODF reconstruction to use \code{c("gqi", "gqi2")} (default: \code{"gqi"}). }
\item{fbase}{ Directory where the required input data files are located. } 
\item{rg}{ range of slices to process (default option \code{rg=c(1,1)}); \code{rg=NULL} processes all slices. }
\item{swap}{ toggle radiological/neurological orientation  (default: \code{FALSE}). }
\item{lambda}{ diffusion sampling length in \code{gdi="gqi"} and \code{gdi="gqi2"}.
 By default the following default values are used when \code{lambda=NULL} is specified: 1.24 in "gqi", 3 in "gqi2". }
\item{depth}{ sampling density on the hemisphere used in simulation (default N=321; depth=3). }
\item{btoption}{ b-table selection between \file{btable.txt} (\code{btoption=1}), and the 3D-DSI grid b-table extracted from the diffusion data set (\file{data.bvec} and \file{data.bval}). By default, the 3D-DSI grid b-table is used (\code{btoption=2}). }
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
\item{texture}{ name of the PNG file to be used as RGB map overlay in some 'showimage' options (default \code{NULL} - no texture). }
\item{...}{ additional material properties for geometry appearance as specified in \code{rgl.material}. } 
}

\details{
GQI (Yeh et.al. 2010) or GQI2 (Garyfallidis 2012) may be used for ODF reconstruction.
Slice glyph map display and overlay selection is controlled by specifying one the arguments \cr
\code{c("none", "gfa", "glyph", "glyphgfa", "glyphrgbmap", "glyphdata")} \cr
 for \code{showimages}.
Meanings are as follows: \code{"none"} - no overlay; \code{"gfa"} - GFA map only; \code{"glyph"} - glyph map only; \code{"glyphgfa"} - GFA overlayed on glyph map; \code{"glyphrgbmap"} - glyphs overlayed on RGB map (if available); \code{"glyphdata"} - \file{data_brain.nii.gz} is overlayed on line map.
}

\references{
Adler, D., and Murdoch, D. \pkg{rgl}: 3D visualization device system (OpenGL), 2012. R package version 0.92.880.

Auguie, B. \pkg{gridExtra}: functions in Grid graphics, 2012. R package version 0.9.1.

Barber, C. B., Habel, K., Grasman, R., Gramacy, R. B., Stahel, A., and Sterratt, D. C. \pkg{geometry}: Mesh generation and surface tesselation, 2012. R package version 0.3-2.

Zeileis, A., Hornik, K., and Murrell, P. Escaping RGBland: Selecting colors for statistical graphics. Computational Statistics & Data Analysis 53 (2009), 3259-3270.
} 

\author{ Adelino Ferreira da Silva, Universidade Nova de Lisboa,
  Faculdade de Ciencias e Tecnologia, Portugal, \email{afs at fct.unl.pt}
}

\seealso{
  \code{\link{gqi.odfpeaks}},
  \code{\link{gqi.odfvmf}},
  \code{\link{gqi.odfvmflines}},
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
  gqi.odfvxgrid(threshold=0.5, depth=2, showimage="glyphgfa")
  gqi.odfvxgrid(gdi="gqi", threshold=0.5, depth=2, showimage="glyphgfa")
  ##----
  ## generate slice RGB map as texture
  texturefname <- file.path(tempdir(),"rgbmap.png")
  rgbvolmap(texture=texturefname, bg="transparent")
  ## coronal glyphs + RGBmap
  gqi.odfvxgrid(threshold=0.5, depth=2, showimage="glyphrgbmap",
    texture=texturefname)
}
}

\keyword{ methods }
\keyword{ dplot }

