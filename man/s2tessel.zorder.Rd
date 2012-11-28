\name{s2tessel.zorder}
\alias{s2tessel.zorder}
\concept{Tessellation}
\title{3D Shell Grid Tessellation}
\description{
3D shell grid tessellation of the icosahedron by Delaunay triangulation.
Point are sorted in ZZ to enable hemi-sphere sampling.
}

\usage{
s2tessel.zorder(depth=3, viewgrid=FALSE, saveg=FALSE)
}

\arguments{
\item{depth}{ sampling density on the hemisphere used in simulation (default N=321; depth=3). }
\item{viewgrid}{ logical variable controlling shell grid visualization (default: \code{FALSE}). }
\item{saveg}{ logical variable controlling grid data storage (default: \code{FALSE}). }
}

\details{
Sampling densities of N=81 (depth=2) and N=321 (depth=3) on the hemisphere are often used in ODF profile reconstruction from raw diffusion acquisitions. 
}

\references{
Adler, D., and Murdoch, D. \pkg{rgl}: 3D visualization device system (OpenGL), 2012. R package version 0.92.880.

Barber, C. B., Habel, K., Grasman, R., Gramacy, R. B., Stahel, A., and Sterratt, D. C. \pkg{geometry}: Mesh generation and surface tessellation, 2012. R package version 0.3-2.
} 

\author{ Adelino Ferreira da Silva, Universidade Nova de Lisboa,
  Faculdade de Ciencias e Tecnologia, Portugal, \email{afs at fct.unl.pt}
}

\seealso{
  \code{\link{gqi.odfvmflines}},
  \code{\link{rgbvolmap}},
  \code{\link{gqi.odfpeaks}},
  \code{\link{gqi.odfpeaklines}},
  \code{\link{gqi.odfvxgrid}},
  \code{\link{plotglyph}},
  \code{\link{simulglyph.vmf}},
  \code{\link{simul.fandtasia}},
  \code{\link{simul.simplefield}}
}

\examples{
\dontrun{
  s2tessel.zorder(depth=2, viewgrid=TRUE) # 162 shell points
  s2tessel.zorder(depth=3, viewgrid=TRUE) # 642 shell points
}
}

\keyword{ utilities }

