\name{prostate}
\Rdversion{1.1}
\alias{prostate}
\docType{data}
\title{
Prostate cancer serum 2D DIGE data
}
\description{
The 2D DIGE experiment was carried out on 32 serum samples from patients (18 prostate cancer (Gleason 5) and 14 benign prostatic hyperplasia (BPH)).
300 spots were selected and logrithms form was taken. The spots coordinates are included 
as first two columns for the spots correlation study. It followed by spots expression data.  
}
\usage{data(prostate)}  

\format{

  \describe{
    \item{\code{x}}{x coordinate}
    \item{\code{y}}{y coordinate}
    \item{\code{BPH_N}}{BPH patient no. N}
    \item{\code{G5_N}}{Gleason 5 patient no.N}
  }
}

\details{

The prostate data is included in prostate.txt which can be uploaded into digeR.
}
\keyword{datasets}
