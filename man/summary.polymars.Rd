\name{summary.polymars}
\alias{summary.polymars}
\title{
summary.polymars - summarizes a polymars model
}
\usage{
summary.polymars(mars.model)
}
\arguments{
\item{mars.model}{
a model returned from a call to polymars
}}
\value{
The returned object contains information about the fitting steps and the model
selected. The first data frame contains a row for each step of the fitting
procedure. In the columns are: a 1 for an addition step or a 0 for a deletion
step, the size of the model at each step, residual sums of squares (RSS) and 
the generalized cross validation value (GCV), testset residual sums of squares
or testset misclassification, whatever was used for the model selection.




The second data frame, model, contains a row for each basis function of the
model. Each row corresponds to one basis
function (with two possible components). The pred1 column contains the indices
of the first predictor of the basis function. Column knot1 is a possible knot 
in this predictor. If this column is NA, the first component is linear. If 
any of the basis functions of the model is categorical then there will be a 
level1 column. Column pred2 is the possible second predictor involved (if 
it is NA the basis function only depends on one 
predictor). Column knot2 contains the possible knot for the predictor pred2, 
and it is NA when this component is linear. This is a similar format 
to the startmodel argument together with an additional first row corresponding to the 
intercept but the startmodel doesn't use a separate column to specify levels of a 
categorical variable . If any predictor in pred2 is categorical then there will be a level2 
column. The column "coefs" (more than one column in the case of multiple response 
regression) contains the coefficients.


}
\description{
is.polymars and summary.polymars are available. print.polymars defaults to 
summary.polymars.


}
\seealso{
plot.polymars, polymars.fun, summary.polymars, print.polymars


}
\examples{
fit <- polymars(mydata[,1], mydata[,c(2,3,4,5,6)])


summary(fit)
}
% Converted by Sd2Rd version 0.2-a4.
