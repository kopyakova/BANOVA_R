\name{BANOVA.mediation}
\alias{BANOVA.mediation}
\alias{print.BANOVA.mediation}
\title{Mediation analysis based on BANOVA models}
\description{
\code{BANOVA.mediation} conducts mediation/moderated mediation analysis based on various BANOVA models.
}
\usage{
BANOVA.mediation(sol_1, sol_2, xvar, mediator)
\method{print}{BANOVA.mediation}(x, ...)
}
\arguments{
  \item{sol_1}{a BANOVA.* model based on an outcome variable, a causal variable, a mediator and possible moderators}

  \item{sol_2}{a BANOVA.Normal model for the mediator which inlcudes the causal variable and moderators}

  \item{xvar}{the causal variable}
  
  \item{mediator}{the mediator variable}
  
  \item{x}{a BANOVA.mediation object}

  \item{\dots}{additional arguments, currently ignored}
}

\details{
A mediation or moderated mediation analysis (Baron and Kenny 1986; Zao, Lynch and Chen 2010; Zhang, Wedel and Pieters 2008) based on BANOVA models is conducted, in which posterior distributions of the direct effect and indirect effect are calculated based on posterior samples. Means and 95\% credible intervals are reported.
  
}
\value{
  \code{BANOVA.mediation} returns an object of class \code{"BANOVA.mediation"}. The returned object is a list containing:
  \item{dir_effects}{tables of the direct effect}
  \item{m1_effects}{tables of the effct of the mediator on the outcome}
  \item{m2_effects}{tables of the effct of the causal variable on the mediator}
  \item{indir_effects}{tables of the indirect effect}
  \item{xvar}{the name of the causal variable}
  \item{mediator}{the name of the mediator}

}

\references{ 

Baron, R. M. and Kenny, D. A. (1986) \emph{Moderator Mediator Variables Distinction in Social Psychological Research: Conceptual, Strategic, and Statistical Considerations}, Journal of Personality and Social Psychology, Vol. 51, No. 6, pp. 1173-82.

Zhang, J., Wedel,M. and Pieters, R. G.M.  (2009) \emph{Sales Effects of Attention to Feature Advertisements: A Bayesian Mediation Analysis}, Journal of Marketing Research, Vol.46, No.5, pp. 669-681.

Ying, Y. and MacKinnon,D. P.  (2009) \emph{Bayesian Mediation Analysis}, Psychological Methods, Vol. 14, No.4, pp. 301-322.

Zhao, X., John G. L. and Chen, Q. (2010) \emph{Reconsidering Baron and Kenny: Myths and Truths About Mediation Analysis}, Journal of Consumer Research, Vol.37, No.2, pp. 197-206.

Wedel, M. and Dong, C. (2016) \emph{BANOVA: Bayesian Analysis of Variance for Consumer Research}. Submitted.

}

\examples{
data(condstudy)
if(0){
out2 <- BANOVA.Normal(att~cond+pict, ~type, condstudy, condstudy$id, 
burnin = 5000, sample = 1000, thin = 20)
out3 <- BANOVA.Normal(pict~cond, ~type, condstudy, condstudy$id, 
burnin = 5000, sample = 1000, thin = 20)
# (moderated) mediation
BANOVA.mediation(out2, out3, xvar='cond', mediator='pict') 
}

}

