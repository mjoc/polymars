\name{plot.polymars}
\alias{plot.polymars}
\alias{polymars.persp}
\title{
plot.polymars - plots for polymars models
}
\description{
  Plot objects of class \code{"polymars"}
}
\usage{
plot.polymars(polymars.model, predictor, response, predictor2, x,
add=F, n, xyz=F, contour.polymars=F, intercept,\dots)
}
\arguments{
\item{polymars.model}{
a model produced by the function \code{polymars}.
}
\item{predictor}{
the index of a predictor that was used when the polymars model was fit.
For the two dimensional plots, this variable is plotted along the X-axis.
}
\item{predictor2}{
the index of a predictor that was used when the polymars model was fit. For the three 
dimensional plots, this variable is plotted along the Y-axis.  
See xyz.
}
\item{response}{
if the model was fitted to multiple response data the response index
should be specified. By default response = 1.
}
\item{x}{
should be a vector of length equal to the number of predictors in the 
original data set. The values should be in the same order as in the original 
dataset. By default the function uses the median values of the data that was 
used to fit the model.  Although the values for predictor and predictor2 are 
not used, they should still be provided.
}
\item{add}{
should the plot be added to a previously created plot? Works only for two
dimensional plots.
}
\item{n}{
number of plotting points (2 dimensional plot) or plotting points along each
axis (3 dimensional plot). The default is n = 100 for 2 dimensional plots and
n = 33 for 3 dimensional plots.
}
\item{xyz}{
is the plot being made a 3 dimensional plot? If it is xyz should be set = T.
If there is only one response it need not be set, if two numerical values 
accompany the model in the call they will be understood as two predictors
for a 3-d plot. By default a 3-d plot uses the S persp function.
Note: categorical predictors cannot be used for 3 dimensional plots.
}
\item{contour.ploymars}{
if the plot being made a 3 dimensional plot and contour = T a plot using the S
function \code{contour} is being made. Note: categorical
predictors cannot be used for 3 dimensional plots.
}
\item{intercept}{
is by default set = T so the fitted intercepts value in the model is included in 
the plots evaluation. Setting intercept = F evaluates a plot without the model's 
intercept. The intercept may also be given any numerical value which then 
overrides the fitted coefficient for the model.
}
\item{...}{
Inherits arguments from the S functions \code{plot} and \code{persp}.
}}
\value{
This function produces a 2-d plot of 1 predictor and response of a polymars fit
at n equally spaced points or a 3-d plot of two predictors and response of a
polymars fit.  The range of the plot is by default equal to the range of the
particular predictor(s) in the original data, but this can be changed by 
xlim=c(from,to) (and ylim=c(from,to))
}
\seealso{
  \code{\link{polymars}}, \code{\link{predict.polymars}},
  \code{\link{summary.polymars}}
}
\examples{
data(state)
state.pm <- polymars(state.x77[,5], state.x77[,-5], gcv = 1)
plot(state.pm,4,7)
}
% Converted by Sd2Rd version 0.2-a4.




