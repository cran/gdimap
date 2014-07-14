\name{sph.odfvmf}
\alias{sph.odfvmf}
\concept{Diffusion Magnetic Resonance}
\concept{QBI Reconstruction}
\concept{von Mises distributions}
\concept{Orientation Distribution Function}
\concept{RGB maps}

\title{Fibre Orientation Estimation Based on von Mises Distributions with Q-ball Reconstruction}
\description{
Fibre orientations in multiple fibre voxels are estimated using a mixture of von Mises-Fisher (vMF) distributions.
This statistical estimation procedure is used to resolve crossing fibre mappings.
}

\usage{
sph.odfvmf(run=TRUE, fbase=NULL, savedir=tempdir(), rg=NULL, swap=FALSE,
 btoption=2, threshold=0.4, showglyph=FALSE, bview="coronal", order=4,
 clusterthr=0.6, aniso=NULL, ...)
}

\arguments{
\item{run}{ logical variable enabling loading previously processed data (default: \code{TRUE}). } 
\item{fbase}{ Directory where the required input data files are located. } 
\item{savedir}{ directory for saving/loading processed results (default: \code{tempdir()}. }
\item{rg}{ range of slices to process; default option \code{rg=NULL} processes all slices. }
\item{swap}{ toggle radiological/neurological orientation  (default: \code{FALSE}). }
\item{btoption}{ b-table selection between \file{btable.txt} (\code{btoption=1}), and the 3D-DSI grid b-table extracted from the diffusion data set (\file{data.bvec} and \file{data.bval}). By default, the 3D-DSI grid b-table is used (\code{btoption=2}). }
\item{threshold}{ thresholding generalized fractional anisotropy (GFA) value at each voxel (default: 0.4). }
\item{bview}{ MRI slice view selection in \{\code{axial, coronal, sagittal}\} (default: \code{"coronal"}). }
\item{showglyph}{ logical variable controlling visualization of voxel glyphs (default: \code{FALSE}). }
\item{order}{ parameter associated with the order of the spherical harmonics approximation (default: 4). }
\item{clusterthr}{ thresholding orientations based on ODF values at each voxel for directional clustering (default: 0.6). }
\item{aniso}{ anisotropic parameter in the range "[0,1)" or \code{NULL} to use in ODF pos-processing default: \code{NULL}. }
\item{...}{ optional specification of non-default control parameters as detailed in \code{movMF}. }
}

\details{
For directional clustering estimation \code{sph.odfvmf} uses a mixture of 2 and 4 von Mises-Fisher (vMF) distributions that serves as a model for directional ODF profile data, corresponding to multiple fibre orientations.
Statistical orientation estimation in \code{sph.odfvmf} is based on von Mises clustering procedures provided by the R-package \pkg{movMF}, by Kurt Hornik and Bettina Gruen.

Starting with the raw diffusion signal acquired on a grid of q-space, the ODF profile is estimated at each voxel, considering a sampling density of unit vectors on a unit S2 grid.
When a threshold is applied to the estimated ODF at each voxel, the non-thresholded unit vectors provide directional statistics information about the estimated  ODF profile.
The main ODF orientations at each voxel relevant for fibre tracking may be estimated by clustering the non-thresholded unit vectors.
The Q-ball reconstruction method with Aganj regularization as implemented in \pkg{dti} (Tabelow and Polzehl) is used for orientation distribution function (ODF) reconstruction.

The main diffusion data set used in the examples is a DICOM data set provided by the "Advanced Biomedical MRI Lab, National Taiwan University Hospital",  which is included in the "DSI Studio" package,  publicly available from the NITRC repository (\url{http://www.nitrc.org}).
One QBI-compatible b-table defining the acquisition setup on a S2-shell is included in the package. 
The b-table for the shell acquisition used in 
For comparison with GQI, the b-table \file{btable.txt} has been used in the examples. This b-table has 203 points distributed on a S2-shell.

The output files \file{data_V1_sph.nii.gz}, \file{data_V2_sph.nii.gz} and \file{data_gfa_sph.nii.gz} are used to visualize RGB maps through \code{rgbvolmap()} or via the "FSL/fslview" tool.
These files may be used to perform white matter fibre tractography. 
}

\value{
\code{sph.odfvmf} outputs three data files in NIfTI format named \cr
 \file{data_V1_sph.nii.gz}, \file{data_V2_sph.nii.gz}, and
 \file{data_gfa_sph.nii.gz}.
 The first and second main fibre directions per voxel are contained in \file{data_V1_sph.nii.gz}, \file{data_V2_sph.nii.gz}, respectively. The file \file{data_gfa_sph.nii.gz} contains the GFA metric per voxel.
}

\references{
Ferreira da Silva, A. R. Computational Representation of White Matter Fiber Orientations, \emph{International Journal of Biomedical Imaging}, Vol. 2013, Article ID 232143, Hindawi Publishing Corporation \url{http://dx.doi.org/10.1155/2013/232143}.

Ferreira da Silva, A. R. Facing the Challenge of Estimating Human Brain White Matter Pathways. In \emph{Proc. of the 4th International Joint Conference on Computational Intelligence} (Oct. 2012), K. Madani, J. Kacprzyk, and J. Filipe, Eds., SciTePress, pp. 709-714.

Hornik, K., and Gruen, B. \pkg{movMF}: Mixtures of von Mises-Fisher Distributions, 2012. R package version 0.1-0.

Jenkinson, M., Beckmann, C. F., Behrens, T. E., Woolrich, M. W., and Smith, S. M. FSL. \emph{NeuroImage} 62, 2 (2012), 782-790.

Tuch D. S., Q-Ball Imaging, \emph{Magnetic Resonance in Medicine} 52 (2004), 1358-1372. 

Tabelow K., Polzehl J.: \pkg{dti}: DTI/DWI Analysis, 2012. R package version 1.1-0.
}


\author{ Adelino Ferreira da Silva, Universidade Nova de Lisboa,
  Faculdade de Ciencias e Tecnologia, Portugal, \email{afs at fct.unl.pt}
}

\seealso{
  \code{\link{sph.odfvmflines}},
  \code{\link{sph.odfpeaklines}},
  \code{\link{gqi.odfvxgrid}},
  \code{\link{rgbvolmap}},
  \code{\link{sph.odfpeaks}},
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
    ## Generate ODF volumes (QBI volume processing)
    ## for a range of slices using von Mises-Fisher clustering
    sph.odfvmf(showglyph=FALSE, threshold=0.5, savedir=tempdir())
    ## RGB maps for range of slices processed by sph.odfvmf()
    rgbvolmap(fbase=tempdir(), rg=c(1,4), bview="coronal")
    ##-------------
    ## Show reconstructed glyphs in ODF processing 
    ## for first and second main fibre direction determination
    sph.odfvmf(rg=c(1,1), bview="coronal", showglyph=TRUE, threshold=0.5)
}
}

\keyword{ methods }
\keyword{ clustering }


