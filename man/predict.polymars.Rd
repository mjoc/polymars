



\name{predict.polymars}
\alias{predict.polymars}
\title{
predict.polymars - fitted values from a `polymars' model
}
\usage{
predict.polymars(polymars.model, x, classify=F, intercept)


}
\arguments{
\item{polymars.model}{
an object of the class \code{polymars}, usually the result of
a call of the function polymars.


}
\item{x}{
the predictor values at which the fitted values will be computed.  The
predictor values can be in a number of formats. It can take the form of a
vector of length equal to the number of predictors in the original data set
or it can be shortened to the length of only those predictors that occur in
the model, in the same order as they appear in the original data set. 


Similarly, x can take the form of a matrix with the number of columns equal to
the number of predictors in the original data set, or shortened to the
number of predictors in the model.




}
\item{classify}{
if the original call to polymars was for a classification problem and you would
like the classifications for the x-values, set this option = T. Otherwise the
function returns a response column for each class (the highest values in each
row is its class for the case when classify = T).


}
\item{intercept}{
is by default = T so the models intercept coefficient is included in the 
evaluation. Setting intercept = F evaluates without the intercept. The 
intercept may also be given any numerical value which overrides the fitted 
coefficient from the model.




}}
\value{
the returned object is a matrix of fitted values, fitted according to the
model produced by the function polymars. The number of columns in the
returned matrix equals the number of responses in the original call to the 
polymars function.




}
\description{


}
\seealso{
polymars, plot.polymars, summary.polymars
}
\examples{


state.pm<-polymars(state.region,state.x77,knots=15,classify=T,additive=T)
predict(state.pm,x=state.x77,classify=T)


}
% Converted by Sd2Rd version 0.2-a4.




