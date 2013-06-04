\name{sph.odfpeaklines}
\alias{sph.odfpeaklines}
\concept{QBI Reconstruction}
\concept{ODF Line-maps}
\concept{Peak Detection}

\title{Fibre Orientation Mapping Based on Local Peak Detection with QBI Reconstruction}
\description{
\code{sph.odfpeaklines} produces line-maps of ODF profiles for diffusion data slices using
 a regularized spheric harmonics method for ODF reconstruction.
}

\usage{
sph.odfpeaklines(run=TRUE, fbase=NULL, roi=NULL, rg=c(1,1), btoption=1,
 swap=FALSE, threshold=0.4, kdir=2,  zfactor=5, showglyph=FALSE,
 snapshot=FALSE, showimage="linesgfa", bview="coronal", savedir=tempdir(),
 pngfig="odfpeak", bg="white", texturefile=NA, order=4)
}

\arguments{
\item{run}{ logical variable enabling loading previously processed data (default: \code{TRUE}). } 
\item{fbase}{ Directory where the required input data files are located. } 
\item{roi}{ Region of Interest (ROI) to use as mask; default mask (\code{roi=NULL}) uses an all brain mask
 for the supplied data set. }
\item{rg}{ range of slices to process (default option \code{rg=c(1,1)}); \code{rg=NULL} processes all slices. }
\item{btoption}{ b-table selection using the included data file \file{btable.txt} (\code{btoption=1})
 or user defined files. With \code{btoption=2}, diffusion data files \file{data.bvec} and \file{data.bval}
 are selected for input, as used in single-shell QBI acquisition. In this case, the files must be provided by the user.
 (Default \code{btoption=1}). }
\item{swap}{ toggle radiological/neurological orientation  (default: \code{FALSE}). }
\item{threshold}{ thresholding generalized fractional anisotropy (GFA) value at each voxel (default: 0.4). }
\item{kdir}{ maximum number of fibre directions to map (default: 2). }
\item{zfactor}{ parameter controlling z-value in relief overlay maps (default: 5). }
\item{showglyph}{ logical variable controlling visualization of voxel glyphs (default: \code{FALSE}). }
\item{snapshot}{ logical variable controlling line-map snapshots (default: \code{FALSE}). }
\item{showimage}{ object controlling visualization of line-maps (default: \code{"linesgfa"}). \cr
  Alternative options are: \cr
  \code{c("none", "gfa", "lines", "linesgfa", "linesrgbmap", "linesdata")} \cr
  (see Details). }
\item{bview}{ MRI slice view selection in \{\code{axial, coronal, sagittal}\} (default: \code{"coronal"}). }
\item{savedir}{ directory for saving/loading processed results (default: \code{tempdir()}). }
\item{pngfig}{ base-name of the PNG file when \code{snapshot} is \code{TRUE} (default \code{"odfpeak"}). }
\item{bg}{ map background colour (default \code{"white"}) }
\item{texturefile}{ name of the PNG file to be used as RGB map overlay in some 'showimage' options (default \code{NA} - no texture). }
\item{order}{ parameter associated with the order of the spherical harmonics approximation (default: 4). }
}

\details{
The identification of fibre directions is performed by extracting the local maxima of the reconstructed ODF, where this function surpasses a certain threshold.
The Q-ball reconstruction method with Aganj regularization as implemented in \pkg{dti} (Tabelow and Polzehl) is used for orientation distribution function (ODF) reconstruction.
\code{sph.odfpeaklines} implements the standard method of fibre orientation detection.
Local maxima of the reconstructed ODF are located simply by selecting a large number of sampled points on the sphere and searching within a fixed radius neighbourhood. 
For a single main fibre orientation the method performs well.
However, for crossing fibres and other complex fibre configurations the peaks of the ODF profiles identified by the methods do not necessarily match the orientations of the distinct fibre populations.
A more robust method is implemented in \code{sph.odfvmflines}.

Starting with the raw high angular resolution diffusion signal acquired on a S2-shell of q-space, the ODF profile is reconstructed at each voxel, considering a sampling density of unit vectors on a unit S2 shell.
Q-ball imaging (QBI) is used for orientation distribution function (ODF) reconstruction.
For comparison with GQI, the b-table \file{btable.txt} has been used in the examples. This b-table has 203 points distributed on a S2-shell.

Slice map display and overlay selection is controlled by specifying one the arguments\cr
 \code{c("none", "gfa", "lines", "linesgfa", "linesrgbmap", "linesdata")} for \code{showimages}.
 Meanings are as follows:
 \code{"none"} - no visualization;
 \code{"gfa"} - GFA map only;
 \code{"lines"} - line map only;
 \code{"linesgfa"} - GFA overlayed on line map;
 \code{"linesrgbmap"} - lines overlayed on RGB map (if available);
 \code{"linesdata"} - \file{data_brain.nii.gz} is overlayed on line map.
}

\value{
\code{sph.odfpeaklines} produces line-maps of ODF profiles for diffusion data slices.
The line-maps may be overlayed with generalized fractional anisotropy (GFA) relief maps, diffusion data maps or ROI maps.
The file \file{V1list.RData} containing the first main orientation directions for all processed voxels is output for further posterior orientation processing.
}

\references{
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
  sph.odfpeaklines(run=TRUE, showimage="lines")
  ## display line-map overlayed on GFA map 
  sph.odfpeaklines(run=FALSE, showimage="linesgfa")
  ##-------------
  ## Show examples of reconstructed glyphs in ODF processing 
  sph.odfpeaklines(showimage="lines", showglyph=TRUE)
  ##------------
  ## using a ROI overlay
  sph.odfpeaklines(roi="slfcst", showimage="linesgfa")
  ## using data overlay
  sph.odfpeaklines(showimage="linesdata")
}
}

\keyword{ methods }
\keyword{ dplot }
\keyword{ utilities }

