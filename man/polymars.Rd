

\name{polymars}
\alias{polymars}
\title{


Multi-response Multivariate Adaptive Regression Splines


}
\usage{


polymars(responses, predictors, weights, maxsize, gcv=4, additive=F,
startmodel, no.interact, knots, knot.space=3, ts.resp, ts.pred, 
ts.weights, classify, factors, tolerance, verbose=F)






}
\arguments{
\item{responses}{


vector of responses, or a matrix for multiple response regression.
In the case of a matrix each column corresponds to a response and each
row corresponds to an observation. Missing values are not allowed.


}
\item{predictors}{


matrix of predictor variables for the regression. Each column corresponds to a
predictor and each row corresponds to an observation in the same order as
they appear in the response argument. Missing values are not allowed.




}
\item{weights}{


vector of observation weights; if supplied, the algorithm fits to minimize the
sum of the weights multiplied by the squared residuals.  The length of
weights must be the same as the number of observations. The weights must
be nonnegative.


}
\item{maxsize}{
the maximum number of basis functions that the model is allowed to grow to in
the stepwise addition procedure. The default is
min(6*(number.of.observations^(1/3)),number.of.observations/4,100).


}
\item{gcv}{
parameter used to find the overall best model from a sequence of fitted models.
The residual sum of squares of a model is penalized by dividing by the square of  
1-(gcv x model size)/cases.  
A larger gcv value would tend to produce a smaller model. By default gcv is set 
to 4.


}
\item{additive}{
Should the fitted model be additive in the predictors?


}
\item{startmodel}{
the first model that is to be fit by the polymars function. It is either an
object of the class polymars or a model specified by the user.  In that case,
it takes the form of a 4 x n matrix (n = number of basis functions in the starting model excluding the intercept term). Each row corresponds to one basis
function (with two possible components). Column 1 is the index of the first
predictor involved. Column 2 is a possible knot in this predictor. If column
2 is NA, the first component is linear. Column 3 is the possible second
predictor involved (if column 3 is NA the basis function only depends on one
predictor). Column 4 contains the possible knot for the predictor in column 3,
and it is NA when this component is linear.  Example: if a row reads
3 NA 2 4.7, the corresponding basis function is [X3 * (X2-4.7)_+]; if a row
reads 2 4.3 NA NA the corresponding basis function is [(X2-4.3)_+].


A fifth column can be added with 1's and 0's, The 1's specify which basis
functions of the startmodel must be in each model. Thus, these
functions stay in the model during the whole stepwise fitting procedure.


If no startmodel is specified polymars starts with a model that only contains 
the intercept.


 
}
\item{no.interact}{
a matrix used if certain predictor interactions are not allowed in the model.
It is given as a matrix (2 * m), with predictor indices as entries.  The two
predictors of any row cannot have interaction terms with each other.


}
\item{knots}{


defines how the function is to find potential knots for the spline basis
functions.  This can be set to the maximum number of knots you would
like to be considered for each predictor.
Usually, to avoid the design matrix becoming singular the actual number of
knots produced is constrained to at most every third order statistic in any
predictor. This constraint can be adjusted using the \code{knot.space} argument.
It can also
be a vector with the number of potential knots for each predictor.
Again the actual number of knots produced is constrained to be at most every
third order statistic any predictor. 
A third possibility is to provide a matrix where each columns corresponds
to the ordered knots you would like to have considered for that predictor.
This matrix should be filled out to a rectangular data structure with NAs.
The default is min(20,round(ncases/4)) knots per predictor.


When specifying knots as a vector an entry of -1 indicates that the predictor 
is a categorical variable and each unique entry in it's column is treated as a 
level.


 
When specifying knots as a single number or a matrix and there are categorical 
variables these are specified separately as such using the factor argument.




}
\item{knot.space}{
is an integer describing the minimum number of order statistics apart that
two knots can be. Knots should not be too close to insure numerical stability.
The default is 3.


}
\item{ts.resp}{


testset responses for model selection. Should have the same number of columns
as the training set response. A testset can be used for the model selection.
Depending on the value of classify, either the model with the smallest testset
residual sum of squares or the smallest testset classification error is
provided.  Overrides gcv.


}
\item{ts.pred}{


testset predictors. Should have the same number of columns
as the training set predictors.


}
\item{ts.weights}{


testset observation weights. A vector of length equal to the number of cases
of the testset. All weights must be non-negative.


}
\item{classify}{


when the response is discrete (categorical), polymars can be used for
classification. In particular, when classify = T, a discrete response
with K levels is replaced by K indicator variables as response. Model
selection is still being carried out using gcv, except when a testset is
provided, in which case testset misclassification is used to select the best
model.




}
\item{factors}{


used to indicate that certain variables in the predictor set are categorical 
variables. Specified as a vector containing the appropriate predictor 
indices (column numbers of categorical variables in predictors matrix). Factors 
can also be set when the "knots" argument is given as a vector, with -1 as 
the appropriate entries for factors.


}
\item{tolerance}{


for each possible candidate to be added/deleted the resulting residual sums 
of squares of the model, with/without this candidate, must be calculated. 
The inversion of of the "X-transpose by X" matrix, X being the design matrix, 
is done by an updating procedure c.f. C.R. Rao - Linear Statistical Inference 
and Its Applications, 2nd. edition, page 33.  
In the inversion the size of the bottom right-hand entry of this matrix is 
critical. If it's value is near zero or the value of it's inverse is  
almost zero then the 
inversion procedure becomes somewhat inaccurate. The lower the tolerance value the 
more careful the procedure is in selecting candidates for addition to the model 
but it may exclude too conservatively. And the other hand if the tolerance is set
too high a spurious result with a singular or otherwise sub-optimal model may 
occur. By default tolerance is set to 1.0e-5.


}
\item{verbose}{


when set =T, the function will print out a line for each addition or deletion 
stage. For example, " + 8 : 5 3.25 2 NA" means adding interaction basis function 
of predictor 5 with knot at 3.25 and predictor 2 (linear),
to make a model of size 8, including intercept.


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


The returned object also contains the fitted values and residuals of the data used
in fitting the model.




}
\description{
is.polymars and summary.polymars are available. print.polymars defaults to 
summary.polymars.




}
\references{


Friedman, J. H. (1991). Multivariate adaptive regression splines
(with discussion).  The Annals of Statistics, 19, 1-141.


Kooperberg, C., Bose, S. and Stone, C. J. (1997). Polychotomous Regression.
Journal of the American Statistical Association, 92,  - .




}
\seealso{


plot.polymars, predict.polymars, summary.polymars


 
}
\examples{


state.pm<-polymars(state.region,state.x77,knots=15,classify=T)


data<-na.omit(auto.stats)
auto.pm<-polymars(data[1:40,1],data[1:40,-1],factor=c(2,3),ts.pred=data[41:66,-1],ts.resp=data[41:66,1],maxsize=40)


}
% Converted by Sd2Rd version 0.2-a4.


