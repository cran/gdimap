\name{niisetup}
\alias{niisetup}
\concept{Utility}
\title{Initialization utility for writing data arrays as NIfTI files using Rniftilib}

\description{
\code{niisetup} is a utility function used as structure initialization in writing an array as a NIfTI file using \pkg{Rniftilib}
}

\usage{
niisetup(savedir=tempdir(), filename, dim)
}

\arguments{
\item{savedir}{ directory where the file is to be stored (default: \code{tempdir()}.) }
\item{filename}{ name of the NIfTI output file without extension specifications. }
\item{dim}{ dimension of the array to be used for writing (default: \code{tempdir()}.) }
}

\value{
\code{niisetup} retuns a NIfTI object corresponding to the basic initialization step of a nifti.image structure as specified in \pkg{Rniftilib}.
}

\author{ Adelino Ferreira da Silva, Universidade Nova de Lisboa,
  Faculdade de Ciencias e Tecnologia, Portugal, \email{afs at fct.unl.pt}
}

\keyword{ utilities }

