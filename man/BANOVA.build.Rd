\name{BANOVA.build}
\alias{BANOVA.build}
\title{Build BANOVA models}
\description{
\code{BANOVA.build} builds(compiles) BANOVA models.
}
\usage{
BANOVA.build(BANOVA_model)
}
\arguments{
  \item{BANOVA_model}{an object of class \code{"BANOVA.model"}}
}

\value{
  \code{BANOVA.build} returns an object of class \code{"BANOVA.build"}. The returned object is a list containing:
  
  \item{stanmodel}{the compiled Stan model}
  \item{model_name}{the model name}
  \item{single_level}{if the model is a single level model}

}

\examples{
if(0){
model <- BANOVA.model('Poisson', single_level = FALSE) 
Poisson_model <- BANOVA.build(model)
}
}



