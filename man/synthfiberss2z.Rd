\name{synthfiberss2z}
\alias{synthfiberss2z}
\concept{Simulation}
\concept{Diffusion signal simulation}
\concept{Glyph mapping}

\title{Voxel Diffusion Profiles for Multiple Fibre Simulation}
\description{
\code{synthfiberss2z} simulates apparent diffusion coefficient (ADC) profiles in multi-direction, diffusion-weighted MR data, for testing ODF reconstruction and fibre orientation estimation.
}

\usage{
synthfiberss2z(g0, angles=c(20,100), b=3000, S0=1, sigma=NULL,
 logplot=TRUE, pos=c(0,0,0), showglyph=FALSE, new=TRUE, wi=NULL)
}

\arguments{
\item{g0}{ matrix of 3D points on the S2 shell used in simulation. }
\item{angles}{ angles in degrees of fibres to be used in simulation (default: two fibres with angles \code{c(20,100)}). }
\item{b}{ strength of the magnetic diffusion gradient (default b-value=3000). }
\item{S0}{signal intensity without the diffusion weighting (default: 1). }
\item{sigma}{ Rician noise level used in simulation (default \code{NULL}). }
\item{logplot}{ logical variable for selecting log-scale (default \code{TRUE}). }
\item{pos}{ 3D positional coordinate (default \code{c(0,0,0)}). }
\item{showglyph}{ logical variable controlling visualization of voxel glyph (default: \code{TRUE}). }
\item{new}{ starts a new figure if \code{TRUE} (default \code{new=TRUE}). }
\item{wi}{ weight given to fiber's volume fraction. Example for two fibers with different weights \code{wi=c(0.7,0.3)} (default \code{NULL} gives equal weigth to all fibers.) }
}

\value{
\code{synthfiberss2z} plots the diffusion profile and returns the synthesized diffusion signal.
}

\details{
The simulation models the profile of the ADC over the sphere.
Prolate diffusion tensor (DT) white matter profiles are estimated with eigenvalues \{1700, 200, 200\}(x 10^(-6) mm2/s) (see D.C. Alexander, 2002).
Diffusion profiles for crossing fibres are simulated from prolate DTs in equal proportions, where each fibre is represented by a prolate DT. 
Noisy profiles may be simulated by adding Rician noise to the simulated diffusion profile, with a user defined standard deviation level specified as \eqn{\sigma} (SNR=1/\eqn{\sigma}). 
Typically, noise values of SNR~30 are used in simulated dMRI.
}

\references{
Barber, C. B., Habel, K., Grasman, R., Gramacy, R. B., Stahel, A., and Sterratt, D. C. \pkg{geometry}: Mesh generation and surface tessellation, 2012. R package version 0.3-2.

Adler, D., and Murdoch, D. \pkg{rgl}: 3D visualization device system (OpenGL), 2012. R package version 0.92.880.

Alexander, D. C., Barker, G. J., and Arridge, S. R. Detection and
Modeling of Non-Gaussian Apparent Diffusion Coefficient Profiles in Human
Brain Data. Magnetic Resonance in Medicine 48 (2002), 331-340.
} 

\author{ Adelino Ferreira da Silva, Universidade Nova de Lisboa,
  Faculdade de Ciencias e Tecnologia, Portugal, \email{afs at fct.unl.pt}
}

\seealso{
  \code{\link{simulglyph.vmf}},
  \code{\link{plotglyph}}
  \code{\link{gqi.odfvmflines}},
  \code{\link{rgbvolmap}},
  \code{\link{gqi.odfpeaks}},
  \code{\link{gqi.odfpeaklines}},
  \code{\link{gqi.odfvxgrid}},
  \code{\link{simulglyph.vmf}},
  \code{\link{simul.fandtasia}},
  \code{\link{simul.simplefield}}
}

\examples{
\dontrun{
    ## S2 grid
    s2 <- s2tessel.zorder(depth=3)
    g0 <- s2$pc
    ## synthetize diffusion signal (two crossing fibres)
    open3d()
    angles=c(20,100); b=3000
    S <- synthfiberss2z(g0=g0, angles=angles, b=b)
    ## synthetize signal with different volume fractions
    S <- synthfiberss2z(g0=g0, angles=angles, b=b, wi=c(0.7,0.3))
    ## synthesize diffusion signal (three crossing fibres)
    angles <- c(0,60,120); b <- 3000
    S <- synthfiberss2z(g0=g0, angles=angles, b=b)
}
}

\keyword{ utilities }
\keyword{ dplot }

