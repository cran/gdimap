\name{simul.fandtasiaSignal}
\alias{simul.fandtasiaSignal}
\concept{Diffusion Field Simulation}
\concept{ODF visualization}
\concept{Glyph field mapping}
\title{Simulation of Crossing-Fibre Diffusion Profiles}
\description{
\code{simul.fandtasiaSignal} generates a synthesized field of diffusion profiles following Barmpoutis's algorithm. 
}

\usage{
simul.fandtasiaSignal(g, gridsz=32, b=4000, sigma=NULL, savedir=tempdir())
}

\arguments{
\item{g}{ matrix of 3D points on the S2 shell used in simulation. }
\item{gridsz}{ dimension of squared grid to use in simulation (default 32). }
\item{b}{ strength of the magnetic diffusion gradient (default b-value=4000). }
\item{sigma}{ Rician noise level used in simulation; (default \code{NULL}). }
\item{savedir}{ directory for saving/loading processed results (default: \code{tempdir()}). }
}

\details{
The diffusion field simulates a field of dimension \code{gridsz x gridsz} of diffusion profiles for testing crossing fibre orientation mapping.
The procedure \code{simul.fandtasiaSignal} is an adaptation of Barmpoutis' MATLAB code to generate synthetic tensor diffusion-weighted MRI fields.
}

\value{
\code{simul.fandtasiaSignal} returns a field of 32x32 diffusion profiles in NIfTI format.
}

\references{
Ferreira da Silva, A. R. Facing the Challenge of Estimating Human Brain White Matter Pathways. In \emph{Proc. of the 4th International Joint Conference on Computational Intelligence} (Oct. 2012), K. Madani, J. Kacprzyk, and J. Filipe, Eds., SciTePress, pp. 709-714.

Hornik, K., and Gruen, B. \pkg{movMF}: Mixtures of von Mises-Fisher Distributions, 2012. R package version 0.0-2.

Barmpoutis, A. Tutorial on Diffusion Tensor MRI using Matlab. Electronic Edition, University of Florida, 2010,\cr
\url{http://www.mathworks.com/matlabcentral/fileexchange/file_infos/26997-fandtasia-toolbox}.
} 

\author{ Adelino Ferreira da Silva, Universidade Nova de Lisboa,
  Faculdade de Ciencias e Tecnologia, Portugal, \email{afs at fct.unl.pt}
}

\seealso{
  \code{\link{simul.fandtasia}},
  \code{\link{simul.simplefield}}
}

\examples{
\dontrun{
    ## shell grid
    s2 <- s2tessel.zorder(depth=3) 
    simul.fandtasiaSignal(g=s2$pc, gridsz=16, sigma=0.01)
}
}

\keyword{ utilities }

