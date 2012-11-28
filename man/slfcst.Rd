\name{slfcst}
\alias{slfcst}
\docType{data}
\title{A Region-of-Interest (ROI) File for Diffusion MRI Analysis}

\description{
The \file{slfcst.nii.gz} is a region-of-interest (ROI) file that was formed by extracting the superior longitudinal fasciculus (SLF) and corticospinal tract (CST) regions from the main data diffusion set \file{data.nii.gz}.

For reasons of limited CRAN storage space, a trimmed down data set with just 4 coronal slices was selected. Coronal slices {29-32} were extracted from the original data set. Hence, the argument \code{bview="coronal"} is used in the examples, and the range of coronal slices range from 1 to 4. User specified complete data sets may be used without restrictions. 
}

\format{
	The file \file{slfcst.nii.gz} is in gzipped NIfTI format.
	The R-package \pkg{oro.nifti} is required to read NIfTI files.
}

\details{
The extracted regions were registered to the DSI data set using the FSL/FLIRT tool.
The procedure uses the \file{JHU-ICBM-labels-2mm.nii.gz} atlas included in the FSL toolbox.

The ROI file \file{slfcst.nii.gz} depicts brain regions where anatomic white matter fibre crossings are known to exist, forming multiple pathway bundles connected to the cerebral cortex.
}

\references{
Jenkinson, M., Beckmann, C. F., Behrens, T. E., Woolrich, M. W., and Smith, S. M. Fsl. \emph{NeuroImage} 62, 2 (2012), 782-790.

Whitcher, B., Schmid, V. J., and Thornton, A. Working with the DICOM and NIfTI data standards in R. \emph{Journal of Statistical Software} 44, 6 (2011), 1-28.
}

\keyword{ data }
