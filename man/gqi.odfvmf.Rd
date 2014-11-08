\name{gqi.odfvmf}
\alias{gqi.odfvmf}
\concept{Diffusion Magnetic Resonance}
\concept{GQI Reconstruction}
\concept{von Mises distributions}
\concept{Orientation Distribution Function}
\concept{RGB maps}


\title{Fibre Orientation Estimation Based on von Mises Distributions with GQI Reconstruction}
\description{
Fibre orientations in multiple fibre voxels are estimated using a mixture of von Mises-Fisher (vMF) distributions.
This statistical estimation procedure is used to resolve crossing fibre mappings.
}

\usage{
gqi.odfvmf(gdi="gqi", run=TRUE, fbase=NULL, savedir=tempdir(), rg=NULL, swap=FALSE,
 lambda=NULL, depth=3, btoption=2, threshold=0.4, showglyph=FALSE, bview="coronal",
 clusterthr=0.6, aniso=NULL, ...)
}

\arguments{
\item{gdi}{ method of ODF reconstruction to use \code{c("gqi", "gqi2")} (default: \code{"gqi"}). }
\item{run}{ logical variable enabling loading previously processed data (default: \code{TRUE}). } 
\item{fbase}{ Directory where the required input data files are located. } 
\item{savedir}{ directory for saving/loading processed results (default: \code{tempdir()}. }
\item{rg}{ range of slices to process; default option \code{rg=NULL} processes all slices. }
\item{swap}{ toggle radiological/neurological orientation  (default: \code{FALSE}). }
\item{lambda}{ diffusion sampling length in \code{gdi="gqi"} and \code{gdi="gqi2"}.
 By default the following default values are used when \code{lambda=NULL} is specified: 1.25 in \dQuote{gqi}, 3 in \dQuote{gqi2}. }
\item{depth}{ sampling density on the hemisphere used in simulation (default N=321; depth=3). }
\item{btoption}{ b-table selection between \file{btable.txt} (\code{btoption=1}), and the 3D-DSI grid b-table extracted from the diffusion data set (\file{data.bvec} and \file{data.bval}). By default, the 3D-DSI grid b-table is used (\code{btoption=2}). }
\item{threshold}{ thresholding generalized fractional anisotropy (GFA) value at each voxel (default: 0.4). }
\item{bview}{ MRI slice view selection in \{\code{axial, coronal, sagittal}\} (default: \code{"coronal"}). }
\item{showglyph}{ logical variable controlling visualization of voxel glyphs (default: \code{FALSE}). }
\item{clusterthr}{ thresholding orientations based on ODF values at each voxel for directional clustering (default: 0.6). }
\item{aniso}{ anisotropic parameter in the range "[0,1)" or \code{NULL} to use in ODF pos-processing default: \code{NULL}. }
\item{...}{ optional specification of non-default control parameters as detailed in \code{movMF}. }
}

\details{
GQI methods specify an operational sampling scheme in q-space from which the ODF can be estimated.
GQI (Yeh et.al. 2010) or GQI2 (Garyfallidis 2012) may be used for ODF reconstruction.
For directional clustering estimation \code{gqi.odfvmf} uses a mixture of 2 and 4 von Mises-Fisher (vMF) distributions that serves as a model for directional ODF profile data, corresponding to multiple fibre orientations.
Statistical orientation estimation in \code{gqi.odfvmf} is based on von Mises clustering procedures provided by the R-package \pkg{movMF}, by Kurt Hornik and Bettina Gruen.

Starting with the raw diffusion signal acquired on a grid of q-space, the ODF profile is estimated at each voxel, considering a sampling density of unit vectors on a unit S2 grid.
When a threshold is applied to the estimated ODF at each voxel, the non-thresholded unit vectors provide directional statistics information about the estimated  ODF profile.
The main ODF orientations at each voxel relevant for fibre tracking may be estimated by clustering the non-thresholded unit vectors.

The main diffusion data set used in the examples is a DICOM data set provided by the "Advanced Biomedical MRI Lab, National Taiwan University Hospital",  which is included in the "DSI Studio" package,  publicly available from the NITRC repository (\url{http://www.nitrc.org}).
Two b-tables defining the acquisition setup are specified.
One is a b-table for a S2-like grid denoted by \file{btable.txt}. 
The other is the b-table for the 3D-DSI sampling scheme used in the DICOM data acquisition.
This b-table has 203 points uniformly distributed on a 3D grid limited to the volume of the unit sphere.
In both tables, the b-values range from 0 to 4000.
Sampling densities of N=81 (depth=2) and N=321 (depth=3) on the hemisphere are often used in ODF profile reconstruction from diffusion acquisitions. 

The output files \file{data_V1_gqi.nii.gz}, \file{data_V2_gqi.nii.gz}, \file{data_V3_gqi.nii.gz}, and \file{data_gfa_gqi.nii.gz} may be used for probabilistic white matter tractography. These principal diffusion direction (PDD) files retain information about the  'theta' and 'alpha' parameters of the von Mises-Fisher mixture at each voxel. The file  \file{data_V123_gqi.nii.gz} joins all three PDD files in a single NIfTI file.  For visualization purposes via a external tool such as "FSL/fslview" the voxel PDDs must be normalized to the unit sphere beforehand by using \code{niinorm}.
}

\value{
\code{gqi.odfvmf} outputs three data files in NIfTI format named \cr
 \file{data_V1_gqi.nii.gz}, \file{data_V2_gqi.nii.gz}, and
 \file{data_gfa_gqi.nii.gz}.
 The first and second main fibre directions per voxel are contained in \file{data_V1_gqi.nii.gz}, \file{data_V2_gqi.nii.gz}, respectively. The file \file{data_gfa_gqi.nii.gz} contains the GFA metric per voxel.
}

\references{
Ferreira da Silva, A. R., Modeling White-matter Fiber-orientation Uncertainty for Improved Probabilistic Tractography. In \emph{Proc. of the  2nd Int. Congr. on Neurotechnology, Electronics and Informatics, NEUROTECHNIX 2014}, Rome, Italy, (Oct., 2014), A. Londral and P. Encarnacao, Eds., INSTICC, SCITEPRESS, pp 71-78.

Ferreira da Silva, A. R. Computational Representation of White Matter Fiber Orientations, \emph{International Journal of Biomedical Imaging}, Vol. 2013, Article ID 232143, Hindawi Publishing Corporation \url{http://dx.doi.org/10.1155/2013/232143}.

Ferreira da Silva, A. R. Facing the Challenge of Estimating Human Brain White Matter Pathways. In \emph{Proc. of the 4th International Joint Conference on Computational Intelligence} (Oct. 2012), K. Madani, J. Kacprzyk, and J. Filipe, Eds., SciTePress, pp. 709-714.

Hornik, K., and Gruen, B. \pkg{movMF}: Mixtures of von Mises-Fisher Distributions, 2014. R package version 0.2-0.

Yeh, F.-C., Wedeen, V. J., and Tseng, W.-Y. I. Generalized q-Sampling Imaging. \emph{IEEE Transactions on Medical Imaging} 29, 9 (2010), 1626-1635.

Garyfallidis E., Towards an Accurate Brain Tractography, 2012, PhD Thesis, University of Cambridge.

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
  \code{\link{data.bval}},
  \code{\link{data.bvec}},
  \code{\link{btable}}
}

\examples{
\dontrun{
    ## Generate ODF volumes (GQI volume processing)
    ## for a range of slices using von Mises-Fisher clustering
    gqi.odfvmf(depth=2, showglyph=FALSE, threshold=0.5, savedir=tempdir())
    ## RGB maps for range of slices processed by gqi.odfvmf()
    rgbvolmap(fbase=tempdir(), rg=c(1,4), bview="coronal")
    ##-------------
    ## Show reconstructed glyphs in ODF processing 
    ## for first and second main fibre direction determination
    gqi.odfvmf(gdi="gqi", rg=c(1,1), bview="coronal", depth=3,
      showglyph=TRUE, threshold=0.5)
    gqi.odfvmf(gdi="gqi2", rg=c(1,1), bview="coronal", depth=3,
      showglyph=TRUE, threshold=0.5)
    ##-------------
    ## speeded up approximations: hardmax and numeric kappa
    gqi.odfvmf(depth=2, showglyph=FALSE, threshold=0.5, savedir=tempdir(),
      E="hardmax", kappa=20)
    rgbvolmap(fbase=tempdir(), rg=c(1,4), bview="coronal")
}
}

\keyword{ methods }
\keyword{ clustering }


