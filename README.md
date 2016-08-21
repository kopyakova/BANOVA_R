# BANOVA R package
V0.6 improvements:
  1. add checking missing data for independent variables
  2. Fix typos
  3. add warnings for mean centering of numeric variables
  4. print full convergence diag. for the Heidelberg and Welch diagnostic.
  5. Table of means can now print any level of interactions

V0.7 improvements:
  1. exclude numeric variables for predictions 
  2. format change of convergence diag 
  3. table of means and prediction change (exp(mu + sigma^2/2)) for the Poisson model
  4. table of means included in the summary
  5. TODO: median -> mean 
  6. TODO: predictors that have 4 levels or fewer will be considered as factors
