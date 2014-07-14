\name{sph.odfpeaks}
\alias{sph.odfpeaks}

\concept{Diffusion Magnetic Resonance}
\concept{GQI Reconstruction}
\concept{Orientation Distribution Function}
\concept{RGB maps}

\title{Main Fibre Orientation Determination via Peak Detection with Q-ball Reconstruction}

\description{
The Q-ball reconstruction method with Aganj regularization as implemented in \pkg{dti} (Tabelow and Polzehl) is used for orientation distribution function (ODF) reconstruction.
For identifying voxel fibre directions, the local maxima of the reconstructed ODF are identified,
 where this function surpasses a certain threshold.
}

\usage{
sph.odfpeaks(fbase=NULL, rg=NULL, swap=FALSE, btoption=2,
 threshold=0.4, showglyph=FALSE, bview="coronal", savedir=tempdir(), order=4)
}

\arguments{
\item{fbase}{ Directory where the required input data files are located. } 
\item{rg}{ range of slices to process; default option \code{rg=NULL} processes all slices. }
\item{swap}{ toggle radiological/neurological orientation  (default: \code{FALSE}). }
\item{btoption}{ b-table selection between \file{btable.txt} (\code{btoption=1}),
 and 3D-DSI grid b-table supplied with the diffusion data set (default \code{btoption=2}). }
\item{threshold}{ thresholding generalized fractional anisotropy (GFA) value at each voxel (default: 0.4). }
\item{showglyph}{ logical variable controlling visualization of voxel glyphs (default: \code{FALSE}). }
\item{bview}{ MRI slice view selection in \{\code{axial, coronal, sagittal}\}
 (default: \code{"coronal"}). }
\item{savedir}{ directory for saving/loading processed results (default: \code{tempdir()}). }
\item{order}{ parameter associated with the order of the spherical harmonics approximation (default: 4). }
}

\details{
The most commonly used approach for identifying fibre directions is to extract the local maxima of the reconstructed ODF, where this function surpasses a certain threshold.
The application of the approach is based on the assumption that the principal directions extracted from the ODF can be interpreted as principal directions of the underlying fibre architecture.
Thresholding avoids selecting smaller ODF peaks that may appear due to noise.
Local maxima of the reconstructed ODF are located simply by selecting a large number of sampled points on the sphere and searching within a fixed radius neighbourhood. 

Starting with the raw high angular resolution diffusion signal acquired on a S2-shell of q-space, the ODF profile is reconstructed at each voxel, considering a sampling density of unit vectors on a unit S2 shell.
Q-ball imaging (QBI) is used for orientation distribution function (ODF) reconstruction.
For comparison with GQI, the b-table \file{btable.txt} has been used in the examples. This b-table has 203 points distributed on a S2-shell.

The output files \file{data_V1_sph.nii.gz} and \file{data_gfa_sph.nii.gz} are used to visualize RGB maps through \code{rgbvolmap()} or via the "FSL/fslview" tool.
These files may also be used to perform white matter fibre tractography.
}

\value{
\code{sph.odfpeaks} outputs two data files in NIfTI format named \cr
 \file{data_V1_sph.nii.gz} and \file{data_gfa_sph.nii.gz}.
The first main fibre directions per voxel are contained in \file{data_V1_sph.nii.gz}.
The file \file{data_gfa_sph.nii.gz} contains the GFA values per voxel.
}

\references{
Ferreira da Silva, A. R. Computational Representation of White Matter Fiber Orientations, \emph{International Journal of Biomedical Imaging}, Vol. 2013, Article ID 232143, Hindawi Publishing Corporation \url{http://dx.doi.org/10.1155/2013/232143}.

Ferreira da Silva, A. R. Facing the Challenge of Estimating Human Brain White Matter Pathways. In \emph{Proc. of the 4th International Joint Conference on Computational Intelligence} (Oct. 2012), K. Madani, J. Kacprzyk, and J. Filipe, Eds., SciTePress, pp. 709-714.

Jenkinson, M., Beckmann, C. F., Behrens, T. E., Woolrich, M. W., and Smith, S. M. FSL. \emph{NeuroImage} 62, 2 (2012), 782-790.

Tuch D. S., Q-Ball Imaging, \emph{Magnetic Resonance in Medicine} 52 (2004), 1358-1372. 

Tabelow K., Polzehl J.: \pkg{dti}: DTI/DWI Analysis, 2012. R package version 1.1-0.
}
 
\author{ Adelino Ferreira da Silva, Universidade Nova de Lisboa,
  Faculdade de Ciencias e Tecnologia, Portugal, \email{afs at fct.unl.pt}
}


\seealso{
  \code{\link{sph.odfpeaklines}},
  \code{\link{plotglyph}},
  \code{\link{sph.odfvmf}},
  \code{\link{sph.odfvmflines}},
  \code{\link{gqi.odfvxgrid}},
  \code{\link{s2tessel.zorder}},
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
    ## Generate ODF volumes (GQI volume processing) for
    ## range of slices using deterministic find peaks algorithm
    sph.odfpeaks(rg=c(1,1), bview="coronal", showglyph=FALSE)
    ## Show RGB maps for range of slices processed by sph.odfpeaks()
    rgbvolmap(fbase=tempdir())
    ##-------------
    ## Process whole volume: generate NIfTI files
    ##   data_gfa_sph.nii.gz and data_V1_sph.nii.gz
    sph.odfpeaks(rg=NULL)
    rgbvolmap(fbase=tempdir(), rg=NULL, bview="coronal")
    ##-------------
    ## Show reconstructed glyphs in ODF processing 
    ## for principal direction determination
    sph.odfpeaks(rg=c(1,1), bview="coronal", showglyph=TRUE)
}
}

\keyword{ methods }
\keyword{ utilities }

