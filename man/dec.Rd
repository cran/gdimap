\name{dec}
\alias{dec}
\concept{Directional color code}

\title{Directionally-Encoded Color (DEC) Representation}

\description{
\code{dec} illustrates the RGB encoding used in DEC representations
}

\usage{
dec(depth=4, new=TRUE)
}

\arguments{
\item{depth}{ code{depth} controls the sampling density on the hemisphere (default \code{depth=4}). }
\item{new}{ start a new plot with \code{new=TRUE} (default). }
}

\details{
The DEC encoding is the standard colouring scheme to represent the orientation of anisotropic tissue, as originally proposed by Pajevic, et.al.
}

\references{
Pajevic, S., and Pierpaoli, C. Color Schemes to Represent the Orien-
tation of Anisotropic Tissues From Diffusion Tensor Data: Application to
White Matter Fiber Tract Mapping in the Human Brain. Magnetic Resonance in Medicine 42 (1999), 526-540.
}


\author{ Adelino Ferreira da Silva, Universidade Nova de Lisboa,
  Faculdade de Ciencias e Tecnologia, Portugal, \email{afs at fct.unl.pt}
}

\examples{
\dontrun{
    dec()
}
}

\keyword{ dplot }

