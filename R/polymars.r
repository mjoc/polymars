# Time-stamp: <97/07/08 16:02:17 martin>  

#******************************************************************************
# (c) Charles Kooperberg and Martin O'Connor 1997                             *
# This function is part of an implementation of the Multivariate Adaptive     *
# Regression Splines (MARS) methodology, first proposed by J.H. Friedman(1991)*
# The Annals of Statistics, 19, 1 - 141.                                      *
# The program is a modified version of MARS (PolyMARS) described in           *
# Kooperberg, C., Bose, S. and Stone C.J. (1997) ``Polychotomous Regression'',*
# 92, 117 - 127.                                                              *
# Journal of the American Statistical Association                             *
# You are free to use this program, for non-commercial purposes only,         *
# under the condition that                                                    *
#                          this note is not to be removed.                    *
#                                                                             *
# The program is not formally maintained, but we are interested in hearing    *
# from people who have problems with it, although we may not be able to solve *
# them.                                                                       *
# Email clk@fhcrc.org or martin@stat.washington.edu                           *
#                       Charles Kooperberg and Martin O'Connor, May 20, 1997  *
#*****************************************************************************/

polymars<-function(responses,predictors,maxsize,gcv=4.0,additive=F,startmodel,weights,no.interact,knots,knot.space=3,ts.resp,ts.pred,ts.weights,classify,factors,tolerance=1e-5,verbose=F)
{
 #responses  - a vector (or matrix) of responses. (Can be a a vector of characters for classification)
 #predictors - a matrix of predictors with same number of cases as response. Columns are predictors.
 #OPTIONAL ARGUEMENTS
 #maxsize    - maximum number of basis function the model can contain 
 #gcv        - parameter for overall best model seletion
 #additive   - boolean, is the model to be additive
 #startmodel - either a matrix (m*4 or m*5) or a polymars object from a previous call to polymars 
 #             an initial model the procedure should start with in model selection
 #weights    - a vector of length equal to the number of cases
 #no.interact- a 2*l matrix of columns numbers of the predictor matrix( each row pair cannot 
 #              have interaction terms)
 #knots      - a vector specifying many knots per predictor are wanted (with -1 for categorical 
 #             variables) ncol(predictors)==length(knots), or a matrix with ncol(predictors) == 
 #             ncol(knots) with actual knot specified and filled out with NA's.
 #             Can also be a single number - "knots" number of knots per predictor                    
 #knot.space - minimum number of order statistics between knots
 #ts.resp    - testset reponses, same format as responses
 #ts.pred    - testset predictors, same format as predictors
 #ts.weights - testset weights, same format as weights
 #classify   - whether classification is to be done, set = T if the response vector is integer, if 
 #             if character classify is automatically true
 #factors    - a vector of column numbers of the predictor matrix of categorical variables
 #tolerance  - a numerical parameter which may need to be made smaller if the program crashes

 #store the call to the mars function
 call <- match.call()

 responses<-as.matrix(responses)
 predictors<-data.matrix(predictors)
 nresponses<-ncol(responses)
 npredictors<-ncol(predictors)
 ncases<-nrow(predictors)
 

 if(missing(classify))classify<-F
 if(mode(responses)=="character" || classify == T)
  {
   if(ncol(responses) > 1)
    {
     stop("When using character responses  or classify = T only 1 response per case is allowed\n")
     
    }
   char.responses<-responses
   int.responses<-as.integer(as.factor(responses))
   
   nresponses<-length(unique(responses))
   responses<-matrix(ncol=nresponses,nrow=ncases,data=int.responses)
   for(i in 1:nresponses)
    {
     responses[,i]<-(responses[,i] == (unique(int.responses)[i]))
    }
   conversion <- matrix(ncol=2,nrow=nresponses,c(unique(char.responses),unique(int.responses)))

   classify<-T
   
   if(!missing(ts.resp))
    {
    
     char.responses.test<-ts.resp

     ts.resp<-matrix(ncol=nresponses,nrow=length(char.responses.test),data=0)

     for(i in 1:nresponses)
      {
       

       ts.resp[,i]<- as.integer(char.responses.test == conversion[i,1])
      }
     
    }


  }
 else
  {
   conversion <- F
   classify<-F
  }

 #maxsize that the model can grow to
 if(missing(maxsize))maxsize<-ceiling(min(6*(ncases**(1/3)),ncases/4,100))
 
 #if a testset is to be used in model selection
 if(!missing(ts.resp) || !missing(ts.pred))
  {
   if(missing(ts.resp) || missing(ts.pred))
    {
     stop("Both ts.resp (testsets responses) and ts.pred (testset predictors) should be specified\n")
    }
   if(!is.matrix(ts.resp))ts.resp<-as.matrix(ts.resp)
   if(!is.matrix(ts.pred))ts.pred<-as.matrix(ts.pred)
   if(ncol(ts.resp) != nresponses)
    {
     stop("Testset should have the same number of responses as the training set\n ")
     
    }
   if(ncol(ts.pred) != npredictors)
    {
     stop("Testset should have the same number of predictors as the training set\n ")
     
    }
   if(nrow(ts.resp) != nrow(ts.pred))
    {
     stop("Testset ts.pred and ts.resp should have the same number of cases (rows)");
    }
   testsetmatrix<-cbind(ts.resp,ts.pred)
   testsetcases<-nrow(testsetmatrix)
   testset<-T

   if(!missing(ts.weights))
    {
     if(length(ts.weights) != testsetcases)
      {
       stop("length of testset weights misspecified\n")
       
      }
     testset.weighted<-T
    }
   else
    {
     testset.weighted<-F
     ts.weights<-0
    }
  }
 else
  {
   testsetmatrix<-0
   testsetcases<-0
   testset<-F
   testset.weighted<-F
   ts.weights<-0
  }
   

 #If the mesh is specified by the knots arguement this will be changed to
 #true later
 mesh.specified<-F
 mesh.vector<-0
 if(nrow(responses) != nrow(predictors))
  {
   
   stop("The number of rows (cases) of the response and predictor matricies should be the same")
  }
 
 if(!missing(knots) && !is.matrix(knots) && length(knots) != npredictors && length(knots) !=1)
  {
   
   stop("Length of vector of `knots per predictor\nshould be equal to number of predictors or 1\n")
   
  }


 if(!missing(knots))
  {
   if(!is.matrix(knots))
    {
     #if knots is specified as a single number it is expanded to a vector 
     #length npredictors
     if(length(knots) ==1){knots<-rep(knots,npredictors)
     if(!missing(factors))
      {
       for(i in 1:length(factors))
        {
         if(! is.vector(factors))
          {
           stop("`factors' should be a vector whose elements are indicies of predictors that are factors\n")
	   
          } 
          # in knots the number of knots(per predictor) is specified
          # or -1 if the predictor is a factor and all it values are levels  
          knots[factors[i]] <- -1
         }
       }
     }
    }
   else
    {
     mesh<-knots
     mesh.vector<-vector(length=ncol(mesh)*nrow(mesh),mode="double")
     knots<-vector(length=npredictors,mode="integer")
     k<-0
     for(i in 1:npredictors)
      {
       knots[i]<-length(unique(mesh[is.na(mesh[,i])==F,i]))
       for(j in 1:knots[i])
        {
         k<-k+1 
         mesh.vector[k]<-unique(mesh[!is.na(mesh[,i]),i])[j]
        }
      }
     if(!missing(factors))
      {
       for(i in 1:length(factors))
        {
         if(! is.vector(factors))
          {
           stop("`factors' should be a vector whose elements are indicies of predictors that are factors\n")
	   
          } 
         # in knots the number of knots(per predictor) is specified
         # or -1 if the predictor is a factor and all it values are levels  
         knots[factors[i]] <- -1
        }
      }
     mesh.specified <-T
    }
   } 

  if(missing(knots))
   {
   knots<-rep(min(20,round(ncases/4)), npredictors)
   if(!missing(factors))
    {
     for(i in 1:length(factors))
      {
       if(! is.vector(factors))
        {
         stop("`factors' should be a vector whose elements are indicies of predictors that are factors\n")
	 
        } 
         # in knots the number of knots(per predictor) is specified
         # or -1 if the predictor is a factor and all it values are levels  
         knots[factors[i]] <- -1
      }
    }  
  }
  startmodelsize<-1
  #A starting model must be specified as a object of class polymars
  #or a matrix with 4 or 5 columns
  no.remove<-0
  no.remove.size<-0


  if(!missing(startmodel))
   {
    
    if(is.vector(startmodel))startmodel<-t(as.matrix(startmodel))
    
    if(!(is.matrix(startmodel) || is.polymars(startmodel)) || 
      (is.matrix(startmodel) && (ncol(startmodel) != 4 &&  
      (ncol(startmodel) != 5))))
     {
      stop(
       paste("startmodel should be a matrix with each row corresponding to\n",
           "a function with number of columns = 4 (or 5 for extra boolean)\n",
           "column specifying predictors which cannot be removed \n"))
      
     }
   
    if(is.matrix(startmodel))
     {
      #Fifth column denotes which basis functions must remain in the model at 
      #all times
      if(ncol(startmodel) == 5)
       {
        no.remove<-vector(length=(nrow(startmodel)))
        j<-0; 
        for(i in 1:nrow(startmodel))
         {
          if(startmodel[i,5]==T)
           {
            j<-j+1
            no.remove[j]<-i
           }
         }
        no.remove.size<-j
       }
      
     
      #The startknots are taken from the startmodel and put into a vector
      #The startmodel becomes a 4*n matrix with a "1" in the 2nd and 4th 
      #columns where knots appear
      startknots<-as.vector(t(cbind(startmodel[,2],startmodel[,4])))
      
      startknots[is.na(startknots)]<-0.0
      
      startmodel<-matrix(startmodel[,1:4],ncol=4)

      startmodel[!is.na(startmodel[,2]),2]<-1

     startmodel[is.na(startmodel[,2]),2]<-0

    startmodel[is.na(startmodel[,3]),3]<-0
    startmodel[startmodel[,3]==0,4]<-0
      for(i in 1:nrow(startmodel))
       {
        if((!is.na(startmodel[i,4])) && startmodel[i,3]!=0)
        startmodel[i,4]<-1
       }
      startmodel[is.na(startmodel[,4]),4]<-0
      startmodelsize<-nrow(startmodel)+1

     }
    else
     {
      startmodelsize<-startmodel$model.size
      
     startmodel<-startmodel$model[-1,]
     startknots1<-startmodel$knot1
      startknots2<-startmodel$knot2
     
      L1<-F
     

      if(!is.null(startmodel$level1))
       {
        L1<-T
        level1<-startmodel$level1
       }
     

      
      if(L1)
       {
       
        startmodel$knot1[!is.na(level1)]<-1
       
        startknots1[!is.na(level1)]<-level1[!is.na(level1)]
       }
    
     

      startknots<-cbind(startknots1,startknots2)
      
      
      startknots<-as.vector(t(startknots))
      startknots[is.na(startknots)]<-0.0
      startmodel<-cbind(startmodel[,"pred1"],
                        startmodel[,"knot1"],
                        startmodel[,"pred2"],
                        startmodel[,"knot2"])
      startmodel[,2]<-!is.na(startmodel[,2])
      startmodel[,4]<-!is.na(startmodel[,4])
     }

    }
   else
    {
     startmodel<-0
     startknots<-0
    }


 if(!missing(weights)) 
  {
   if(length(weights) != ncases)
    {
     stop("Number of weights not equal to the numnber of cases\n")
     
    }
   weighted<-T
  } 
 else
  {
   weighted<-F
   weights<-0
  }
 
 datamatrix<-cbind(responses,predictors)

 #Predictors which cannot interact together in the model are specified 
 #by a 2*n matrix of predictor indicies
 if(!missing(no.interact))
  {
   if(!is.matrix(no.interact) || ncol(no.interact)!=2) 
    {
     stop("list of interactions disallowed has been misspecified,must be a 2*n matrix")
    }
   no.interact<-t(no.interact)
   no.interact.size<-ncol(no.interact)
  } 
 else 
  {
   no.interact.size<-0
   no.interact<-0
  }
 
 if(startmodelsize > maxsize)
  {
   stop("start model should not be of greater size than the max model size\n")
   
  }

 #Some error checking on the startmodel
 
 if(startmodelsize != 1)
  {
   for(i in 1:(startmodelsize-1))
    {
     if(startmodel[i,1] ==0)
      {
       stop("first column of startmodel cannot be zero\n")
       
      } 
     
     if(startmodel[i,2] == 1)
      {
     	if(startknots[(i*2)-1] < min(predictors[,startmodel[i,1]]) 
          || startknots[(i*2)-1] > max(predictors[,startmodel[i,1]]))
        {
         stop("Knot out of range of its predictor \n")
        }	
      }
    
    
     if(startmodel[i,4] == 1)
      {
       
       if(startknots[(i*2)] <= min(predictors[,startmodel[i,3]]) 
           || startknots[(i*2)] >= max(predictors[,startmodel[i,3]]))
        {
         
         stop("Knot out of range of its predictor\n")
         
        }
      }
    }
  
   if(max(startmodel[,c(1,3)]>npredictors))
    {
     stop("Initial model misspecified on input\n")
     
    }
  }
 
 
 
 startmodel<-t(startmodel) 
 resultmodelsize<-0
 end.state<-0
 step.count <-0
 
 z <- .C("polymars",
         as.integer(npredictors),
         as.integer(nresponses),
         as.integer(ncases),
	 as.double(datamatrix),
         as.integer(knots),
         as.double(mesh.vector),
         as.integer(mesh.specified),
         as.integer(maxsize),
         as.double(gcv),
	 as.integer(additive),
         as.integer(startmodelsize),
         start.model=as.integer(startmodel),
         start.knots=as.double(startknots),
         as.integer(weighted),
         as.double(weights),
         as.integer(no.interact.size),
         as.integer(no.interact),
         as.integer(no.remove.size),
         as.integer(no.remove),
         as.integer(knot.space),
         as.integer(testset),
         as.double(testsetmatrix),
         as.integer(testsetcases),
         as.integer(testset.weighted),
         as.double(ts.weights),
         as.integer(classify),
         as.double(tolerance),
         as.integer(verbose),
         best.model = as.integer(matrix(nrow=maxsize,
                                        ncol=4,
                                        data=rep(0,maxsize*4))),
         coefficients = as.double(matrix(nrow=maxsize,
                                        ncol=nresponses,
                                        data=rep(0.0,maxsize*nresponses))),
         steps = as.integer(matrix(nrow=maxsize*2,
                                   ncol=2,
                                   data=rep(0,maxsize*4))),
         rss.gcv = as.double(matrix(nrow=maxsize*2,
                                    ncol=nresponses+1,
                                    data=rep(0.0,maxsize*2*(nresponses+1)))),
         modelsize=as.integer(resultmodelsize),
	 modelknots = as.double(matrix(nrow=maxsize,
                                       ncol=2,
                                       data=rep(0.0,maxsize*2))),
         coefficient.se.term = as.double(rep(0.0,maxsize)),
	 end.state = as.integer(end.state),
         step.count = as.integer(step.count))
 
 #The C function returns information about how it ended
 
 if(z$end.state != 0 && z$end.state !=5)
  {
   
   switch(z$end.state,
   stop("Mis-specification of initial model\n"),
   stop("Initial model with non-linear function must contain the corresponding linear function\n"),
   stop("Initial model contains two-predictor functions that require prerequisite functions\n"))
  
  }
 else
  {
   
   
   
   model<-matrix(z$best.model[1:((z$modelsize-1)*4)],ncol=4,byrow=T)
   knot.values<-matrix(z$modelknots[1:((z$modelsize-1)*2)],ncol=2,byrow=T)
   
   for(i in 1:nrow(model))
    {
     if(model[i,2] != 0){model[i,2]<-knot.values[i,1]}else{model[i,2]<-NA}
     if(model[i,4] != 0){model[i,4]<-knot.values[i,2]}else{model[i,4]<-NA}
    }
  

  
  if(length(knots[model[,1]]) !=0 && min(knots[model[,1]])<0)
   {
    
    factor1<-T
    levels1<-rep(NA,z$modelsize-1)
    
    factor.variables<-unique(model[knots[model[,1]]<0,1])
    
    
    for( i in 1:length(factor.variables))
     {
      for( j in 1:length(model[,1]))
       {	
        if(model[j,1] == factor.variables[i]){levels1[j]<-model[j,2];}
       }
      model[model[,1] == factor.variables[i],2]<-NA
     }
    levels1<-c(NA,levels1)
    
   } 
  else 
   {factor1<-F}
  

  
  

  coefs<-matrix(z$coefficients[1:(z$modelsize*nresponses)],ncol=nresponses)
   #The model that the C-function returns does not explicitly contain an intercept
   #so in formatting the output one is added
   
   
   if(z$modelsize > 1)
    { 
     
     
     if(factor1==F)
      {
       model<-rbind(c(0,NA,0,NA),model)
       model<-data.frame(model,coefs)
       if(nresponses == 1)
        {
         dimnames(model)<-list(1:z$modelsize,c("pred1","knot1","pred2",
                               "knot2","coefs"))
        }
       else
        {
         dimnames(model)<-list(1:z$modelsize,c("pred1","knot1","pred2","knot2",
                               paste("Coefs",1:nresponses)))
        }
      }
    
      
     if(factor1==T)
      {
              
       
       model[(knots[model[,1]]<0),2]<-NA
       model<-rbind(c(0,NA,0,NA),model)
       model<-data.frame(model[,1:2],levels1,model[,3:4],coefs)
       if(nresponses == 1)
        {
         dimnames(model)<-list(1:z$modelsize,c("pred1","knot1","level1",
                                               "pred2","knot2","coefs"))
        }
       else
        {
         dimnames(model)<-list(1:z$modelsize,c("pred1","knot1","level1",
                               "pred2","knot2",paste("Coefs",1:nresponses)))
        }
      }
      
       
    }
   else
    {
     
     
     model<-data.frame(0,NA,0,NA,coefs)
     
     if(nresponses == 1)
      {
       dimnames(model)<-list(1:z$modelsize,c("pred1","knot1","pred2",
                                           "knot2","coefs"))
      }
     else
      {
       dimnames(model)<-list(1:z$modelsize,c("pred1","knot1","pred2","knot2",
                                           paste("Coefs",1:nresponses)))
      }
    }
 
   #for later plotting the ranges and medians of the predictors are stored
   ranges.and.medians<-matrix(ncol=npredictors,nrow=3,data=0)
   
   
   
   for(i in 1:npredictors)
    {ranges.and.medians[1,i]<-min(predictors[,i])}
   for(i in 1:npredictors)
    {ranges.and.medians[2,i]<-max(predictors[,i])}
   for(i in 1:npredictors)
    {ranges.and.medians[3,i]<-median(predictors[,i])}

   # A table with information from the fitting is formatted here
   steps<-matrix(z$steps[1:(2*(z$step.count+1))],ncol=2,byrow=T)
   rss.gcv<-matrix(z$rss.gcv[1:((nresponses+1)*(z$step.count+1))],
                   ncol=nresponses+1,
                   byrow=T)
   fitting<-data.frame(steps,rss.gcv)
  
   if(testset == F)
    {
     if(nresponses == 1)
      {
       dimnames(fitting) <- list(1:(nrow(fitting)), 
                                  c("0/1","size","RSS","GCV"))
      }
     else
      {
       dimnames(fitting) <- list(1:nrow(fitting), 
                                 c("0/1",
                                 "size", 
                                 paste("RSS", 1:nresponses), 
                                 "GCV"))
      }
    }
   else
    {
     if(classify ==F)
      {
       if(nresponses == 1)
        {
         dimnames(fitting) <- list(1:(nrow(fitting)), 
                                  c("0/1","size","RSS","T.S. RSS"))
        }
       else
        {
         dimnames(fitting) <- list(1:nrow(fitting), 
                                   c("0/1",
                                   "size", 
                                   paste("RSS", 1:nresponses), 
                                   "T.S. RSS"))
        } 
      }
    else
     {
      if(nresponses == 1)
       {
        dimnames(fitting) <- list(1:(nrow(fitting)), 
                                  c("0/1","size","RSS","T.S.M.C."))
        }
       else
        {
         dimnames(fitting) <- list(1:nrow(fitting), 
                                   c("0/1",
                                   "size", 
                                   paste("RSS", 1:nresponses), 
                                   "T.S.M.C."))
        } 
      }
    }
 

   #calculates fitted values and residual of the data according to the
   #model returned 
   if(z$modelsize >1)
    {
     temp<-list(model=model,
                model.size = z$modelsize,
                ranges.and.medians=ranges.and.medians,
                responses = nresponses
                )
     class(temp)<-"polymars"
     fitted<-matrix(ncol=nresponses,
                    nrow=ncases,
                    data=rep(0,nresponses*ncases))
     residuals<-matrix(ncol=nresponses,
                       nrow=ncases,
                       data=rep(0,nresponses*ncases))
     
       fitted<-predict.polymars(temp,x=predictors)
       residuals<-responses-fitted
    #model<-model[,-c(nrow(model)-nreponses,nrow(model))]
    }
   else
    {
     fitted<-matrix(ncol=nresponses,nrow=ncases,data=coefs[1,1])
     residuals<-matrix(ncol=nresponses,nrow=ncases,data=responses-coefs[1,1])
    }


  # if their are factors present in the model the factors must be stored for use during plotting
  if(factor1 == T)
   {
    model2<-model[-1,]
    factors.in.model<-unique(model2[knots[model2[,1]]<0,1])
    
    maxfactors<-0
    for(i in 1:length(factors.in.model))
     {
      maxfactors<-max(maxfactors,length(unique(predictors[,factors.in.model[i]])))
     
     }
    

   
   
   
    factor.matrix<-matrix(ncol=length(factors.in.model),
                        nrow=maxfactors+2,data=NA)
    for(i in 1:length(factors.in.model))
     {
      factor.matrix[1,i]<-factors.in.model[i]
      factor.matrix[2,i]<-length(unique(predictors[,factors.in.model[i]]))
      for(j in 3:(length(unique(predictors[,factors.in.model[i]]))+2))
       {
        factor.matrix[j,i]<-unique(predictors[,factors.in.model[i]])[j-2]
       }
     }
    
   }
  else
   {
    factor.matrix<-0
   }
   
   
   
   
   
   if(nresponses ==1)
    {
   
     SE <- round(sqrt((sum(residuals^2)/(ncases-z$modelsize))*z$coefficient.se.term[1:z$modelsize]),4)
     model<-cbind(model,SE)
   
     dimnames(model)[[2]][length(dimnames(model)[[2]])]<-"SE"
   
    }
   else
    {
     
     for(i in 1:nresponses)
      {
       SE <- round(sqrt((sum(residuals[,i]^2)/(ncases-z$modelsize))*z$coefficient.se.term[1:z$modelsize]),4)
        
       model<-cbind(model,SE)
       
       dimnames(model)[[2]][length(dimnames(model)[[2]])]<-paste("SE",i)
      }

    } 
  
   
   
   
   
   
   
   if(nresponses ==1)
    {
     
     Rsquared<-1-(sum(residuals^2)/sum((responses-mean(responses))^2))
    }
   else
    {
     Rsquared<-NULL
    } 
   result<-list(model=model,
                fitting=fitting,
                model.size = z$modelsize,
                fitted=fitted,
                responses=nresponses,
                residuals=residuals,
	        ranges.and.medians=ranges.and.medians,
                call=call,
                conversion=conversion,
                factor.matrix=factor.matrix,
                Rsquared=Rsquared)
   class(result)<-"polymars"
   return(result)
  }
}


################################################################################################

predict.polymars<-function(mars.model,x,classify=F,intercept)
{
 # Produces predicted values for a polymars object
 # mars.model  an object returned from a call to polymars
 # x           a matrix with number of columns equal to number of columns of predictor matrix in 
 #             original call to polymars and predictor values in the corresponding columns. Can 
 #             also be a matrix with number of column equal to the number of predictors in the 
 #             model, in the order of the original dataset.
 # classify    If the original call to polymars was for classification setting  classify=T will 
 #             the new data otherwise it will return the multi-response fitted values.
 # intercept   By default T. The full intercept is included; or when F the intercept is left out.
 #             Can also be givebn a numerical value
 
 if(missing(intercept))
  {
   intercept<-T
  }
 # some error checking
 if(!(is.polymars(mars.model)))
  {
   stop("The first arguement must be an object returned from a call to `mars'\n")
   
  }
 # The x matrix number of columns can be of length equal to the number of 
 # predictors in the original model or shorten to the number of predictors in 
 # the model in `mars.model'
 if(!(is.matrix(x))) 
  {
   if(length(unique(mars.model$model[, "pred1"]))== 1 ||  ncol(mars.model$ranges.and.medians)== 1  )
    {
	x<-matrix(data=x,ncol=1)
    }
  }
 if((is.matrix(x) && ncol(x) 
    != length(unique(mars.model$model[,"pred1"]))))
  {
   if(ncol(x) != ncol(mars.model$ranges.and.medians))
    {
     
     
     stop("Input should be a matrix with number of columns equal to either number of original predictors or number of predictors in model\n")
     
    }
  }
 
 # If the number of columns of the matrix is not length equal to number of 
 # predictors it is expanded to that size.
 if(is.matrix(x) && ncol(x) == length(unique(mars.model$model[, "pred1"])) && ncol(x) != ncol(mars.model$ranges.and.medians))
  {
   tempmatrix<-x
   
   x<-matrix(nrow=nrow(tempmatrix),ncol=ncol(mars.model$ranges.and.medians),data = 0)
   for(i in 1:length(unique(mars.model$model[, "pred1"]))) 
    {
     for(j in 1:nrow(tempmatrix))
      {
       x[j,sort(unique(mars.model$model[,"pred1"]))[i]]<-x[j]
      }
    }	
  }
 # If x is a vector put it into matrix form expanding it if it is of 
 # length equal to only the number of predictors in the model in `mars.model'
 if(!(is.matrix(x))) 
  {
   if(!(length(x) == ncol(mars.model$ranges.and.medians) || length(x) == unique(mars.model$model[, "pred1"]))) 
    {
     stop("The vector of values must be equal in length to either the number of original predictors or predictors in the model\n")
    
    }
   if(length(x) == unique(mars.model$model[, "pred1"]) && length(x) != ncol(mars.model$ranges.and.medians)) 
    {
     x <- rep(0, ncol(mars.model$ranges.and.medians))
     for(i in 1:length(unique(mars.model$model[, "pred1"]))) 
      {
       x[sort(unique(mars.model$model[, "pred1"]))[i]]<-x[i]
      }
    }
   x <- t(as.matrix(x))
  }
 
 # Checking to see if there are factor variables in the model
 if(dimnames(mars.model$model)[[2]][3] == "level1")
  {
   level1<-T
   mars.model$model<-mars.model$model[,c(1:(5+mars.model$responses))]


   #if(dimnames(mars.model$model)[[2]][6] == "level2"){level2<-T}else{level2<-F}
  }
 else
  {
   level1<-F
   mars.model$model<-mars.model$model[,c(1:(4+mars.model$responses))]
   #if(dimnames(mars.model$model)[[2]][5] == "level2")
   # {level2<-T}else{level2<-F}
  }
 # Setting up the fitted responses matrix
 responses<-mars.model$responses
 Y <- matrix(ncol = responses, nrow = nrow(x), data = rep(0, nrow(x)))
 Y1 <- matrix(ncol = 1, nrow = nrow(x), data = rep(0, nrow(x)))
 Y2 <- matrix(ncol = 1, nrow = nrow(x), data = rep(0, nrow(x)))
 if(is.logical(intercept))
  {
   if(intercept==T)
    {
     for(i in 1:responses)Y[,i] <- mars.model$model[1,ncol(mars.model$model)-responses+i]
     }
   else
     {
      if(intercept==F)
       {
        for(i in 1:responses)Y[,i] <- 0.0
        
       }
     }
  }
 else
  {
   if(is.numeric(intercept))
    {
     if(length(intercept)==responses)
      {
       for(i in 1:responses)Y[,i] <- intercept[i]
      }
     else
      {
       if(length(intercept) != 1)
        {
       	 stop("Intercept arguement mispecified \n")
       	 
        }
       for(i in 1:responses)Y[,i] <- intercept
      }

    }
  }
 # Computing fitted values
 if(mars.model$model.size>1)
  {
   for(i in 2:mars.model$model.size) 
    {
  
     Y2[] <- 1
   
     Y1[] <- x[,mars.model$model[i, "pred1"]]
     if(!is.na(mars.model$model[i, "knot1"])) 
      {
       Y1 <- Y1 - mars.model$model[i,"knot1"]
       Y1[Y1 < 0,] <- 0
      }
     if(level1)
      {

     if(!is.na(mars.model$model[i, "level1"]))
      {
       Y1<- (Y1 == mars.model$model[i, "level1"])
      }
    }
   if(!is.na(mars.model$model[i, "pred2"]) & mars.model$model[i, "pred2"] != 0) 
    {
     Y2[] <- x[,mars.model$model[i,"pred2"]]
     if(!is.na(mars.model$model[i,"knot2" ])) 
      {
       Y2 <- Y2 - mars.model$model[i,"knot2"]
       Y2[Y2 < 0,] <- 0
      }
     #if(level2)
     #{
     #  if(!is.na(mars.model$model[i, "level2"]))
     #	{
     #   Y2<- (Y2 == mars.model$model[i, "level2"])
     #  }
     #}
    }
   
    for(j in 1:responses){Y[,j]<-Y[,j]+(Y1 * Y2 * mars.model$model[i,ncol(mars.model$model)-responses+j])}
   }
  }
 # If classification is to be used the original polymars fitting expanded the 
 # response into a vector of indicator variables. The largest of the responses 
 # correspondes to the fitted class for each case.
 if(classify == T)
  {
   for(i in 1:nrow(Y))
    {
     Y[i,]<-Y[i,]==max(Y[i,])
    }
   if(is.matrix(mars.model$conversion))
   Z<-Y
   
   Y<-matrix(ncol=1,nrow=nrow(Z))
   for(i in 1:nrow(Y))
    {
     for(j in 1:ncol(Z))
      {
       
       if(Z[i,j] == 1) Y[i,] <- mars.model$conversion[j]
      }     
    }
   }
 else
  {
   # if classification was used but the full multiple response fitted response 
   # matrix is requested the response names (corresponding to the classes) are 
   # added.
   if(is.matrix(mars.model$conversion))
    {
     dimnames(Y)[[2]]<-list(mars.model$conversion[,1])

    }
  }
	
 return(Y)
}
################################################################################################
is.polymars<-function(object)
{
if(is.na(class(object)[1]=="polymars")){return(F)}
if(class(object) == "polymars") {return(T)}
else{return(F)}

}
################################################################################################
print.polymars<-function(mars.model)
{
        UseMethod("summary")
}
################################################################################################
summary.polymars<-function(mars.model)
{
        cat("Call:\n")
        print(mars.model$call)
        cat("\nModel fitting\n\n")
        print(mars.model$fitting)
        cat("\n\nModel produced\n\n")
        print(mars.model$model)
        if(mars.model$responses != 1)
                cat("\nRESPONSES :", mars.model$responses, "\n")
        if(!is.null(mars.model$Rsquared))
                cat("\nRsquared :",round(mars.model$Rsquared,3),"\n")
        invisible()
}


plot.polymars<-function(mars.model,predictor1,response,predictor2,x,add=F,n,xyz=F,contour.ploymars=F,xlim,ylim,main,intercept,...)
{
 # mars.model      a polymars object
 # predictor1      the column number in the original predictor matrix of the predictor of interest
 # response        with multi-response polymars the  column number in the original response matrix.
 #                 Default is 1
 # predictor2      the second predictor for a contour or persp plot. For single response data 
 #                 plot(pmars1,2,6) is understood as 3-d plot of predictors 2 and 6.
 # x               Values for the other predictors can be given, using the same format as for the 
 #                 predict fuhnction. By default median values are used.
 # add             should the plot be added to another. (for 2-d plots only)
 # n               For 2-d plot the number of points the function is interploted over. For 3-d plots
 #                 the a n*n mesh is interploted over. Default 2-d: 100, 3-d 33.
 # xyz             sometimes a call can be ambiguous: plot(pmars1,6,2) a 2-d plot with 2nd response 
 #                 or a 3-d plot. Use xyz=T for 3-d.
 #contour.ploymars By default a 3-d if a `persp' plot. contour.ploymars=T asks for a `contour' plot.
 #intercept        same as for predict function. =T intercepr is included =F it is left out, or can
 #                 be given a numerical value.

 if(missing(x))x<-mars.model$ranges.and.medians[3,]
 if(length(x) != ncol(mars.model$ranges.and.medians))
  {
   stop("x should be of length equal to the number of predictors in original data\n")
  
  }
 if(missing(predictor2) && (!missing(response)) && mars.model$responses == 1)
  {
   if(missing(predictor1) && xyz == T)
    {

     stop("You must specify 2 predictor numbers")
     

    }
   xyz<-T
   predictor2<-response
   response<-1
   
  } 
 if(contour.ploymars == T)
  {
   xyz<-T
  }
 
 if(missing(intercept))
  {
   intercept<-T
  }
 
 if(xyz==T)
  {
   
   
   if(missing(n))n<-33
   if(missing(response))
    {
     
     if(missing(xlim))
      {
       polymars.persp(mars.model,
                  predictor1,
                  predictor2,
                  n=n,
                  
                  contour.ploymars=contour.ploymars,
                  intercept=intercept,
                  ...)
       
      }
     else
      {
      
       polymars.persp(mars.model,
                  predictor1,
                  predictor2,
                  n=n,
                  xlim=xlim,
                  
                  contour.ploymars=contour.ploymars,
                  intercept=intercept,
                  ...)
      }
    }
   else
    {
     
     if(missing(xlim))
      {
       polymars.persp(mars.model,
                  predictor1,
                  predictor2,
                  response,
                  n=n,
                  
                  contour.ploymars=contour.ploymars,
                  intercept=intercept,
                  ...)
      }
     else
      {
      polymars.persp(mars.model,
                 predictor1,
                 predictor2,
                 response,
                 n=n,
                 xlim=xlim,
                 contour.ploymars=contour.ploymars,
                 intercept=intercept,
                 ...)
      }
    }
   invisible(return())
   }
 else
  {
   
   if(missing(predictor1))
    {
     cat("predictor should be specified \n")
    }
  if(mars.model$responses != 1 && missing(response)&& missing(predictor2))
   {
    cat("Response should be specified  (default: response =1)\n")
   }
  
	
  #check to see that the predictor is in the model
  inmodel<-F
  for(i in 2:mars.model$model.size)
   {
    if(mars.model$model[i,"pred1"] == predictor1)inmodel<-T
   }
  if(inmodel == F)
  {
   stop("The predictor specified is not in the model\n")
  
  }
  #check to see if the predictor is a factor

  if(is.matrix(mars.model$factor.matrix))
   {
    if(length(mars.model$factor.matrix[1,mars.model$factor.matrix[1,]== predictor1]) != 0)
     {
      isfactor<-T
     }
    else
     { 
      isfactor<-F
     }
   }
  else
   {
    isfactor<-F
   }

  if(isfactor == T)
   { 
   

    pred.values <- matrix(nrow = mars.model$factor.matrix[2,mars.model$factor.matrix[1,]==predictor1], ncol = ncol(mars.model$ranges.and.medians),data = x, byrow = T)
    factors<-mars.model$factor.matrix[-c(1,2),mars.model$factor.matrix[1,]==predictor1]
    pred.values[,predictor1]<- factors[!is.na(factors)]

   
   
   mesh<-factors[!is.na(factors)]	
   }
  else
   {
    if(missing(n))n<-100
    if(missing(xlim))xlim<-c(mars.model$ranges.and.medians[1,predictor1],mars.model$ranges.and.medians[2,predictor1])
    
    pred.values <- matrix(nrow = n, ncol = ncol(mars.model$ranges.and.medians), 
    data = x, byrow = T)
    
    mesh <- matrix(seq(xlim[1],xlim[2],(xlim[2]-xlim[1])/(n-1)),nrow=1)
    pred.values[,predictor1]<-mesh
   }
  if(missing(response) && missing(predictor2))response<-1
  if(response > mars.model$responses || response < 0)
   {
    stop("response arguement = ",response,"is out of range\n")
    
   }

  model<-mars.model$model

 
 Y<-predict.polymars(mars.model,pred.values,intercept=intercept)
 
 if(isfactor == F)
  {
   if(add == F)
    {
     if(mars.model$responses == 1)
      {
       plot(mesh,Y,...,type="l",xlab=paste("Predictor ",predictor1),ylab="Response")
	
      }
     else
      {
	
       plot(mesh,
            Y[,response],
            type="l",
            xlab=paste("Predictor ",predictor1),
            ylab=paste("Response ",response),
            ...)
      }
    }
  else
   {
	
    points(mesh,
           Y,
           type="l")
   }
  }

 if(isfactor == T)
  {
	 
   if(add == F)
    {
     if(mars.model$responses == 1)
      {
       plot(mesh,Y,...,xlab=paste("Predictor ",predictor1),ylab="Response")
     }  
    else
     { 
       plot(mesh,
           Y[,response],
           type="l",
           xlab=paste("Predictor ",predictor1),
           ylab=paste("Response ",response),
           ...)
     }
    } 
   else
    {
     points(mesh,
            Y,
            type="l")
    }
  }





	
  invisible()
 }
}


################################################################################################
polymars.persp<-function(mars.model, predictor1, predictor2, response, n= 33,xlim,ylim,x,contour.ploymars,main,intercept,...)
{
 # used by the plot.polymars function
 # not designed for stand alone use.

 if(missing(x))x<-mars.model$ranges.and.medians[3,]
 if(missing(xlim))xlim<-c(mars.model$ranges.and.medians[1,predictor1],mars.model$ranges.and.medians[2,predictor1])
 if(missing(ylim))ylim<-c(mars.model$ranges.and.medians[1,predictor2],mars.model$ranges.and.medians[2,predictor2])
 if(missing(predictor1) || missing(predictor2)) 
  {
   stop("You must specify 2 predictor numbers\n")
   
  }
 if(mars.model$responses != 1 && missing(response)) 
  {
   cat("Response should be specified  (default: response =1)\n")
  }
 if(missing(response))response <- 1
 if(response > mars.model$responses || response < 0)
  {
   stop("response arguement = ",response,"is out of range\n")
   
  }

 if(sum(as.integer(predictor1==mars.model$model[,1])) == 0)
  {
   stop("Predictor 1 not in model\n")
   
  }
 if(sum(as.integer(predictor2==mars.model$model[,1])) == 0)
  {
   stop("Predictor 2 not in model\n")
  
  }
 

 X <- seq(xlim[1],xlim[2],(xlim[2] - xlim[1])/(n-1))

 y <- seq(ylim[1],ylim[2],(ylim[2] - ylim[1])/(n-1)) 
 meshX <- rep(X, n)
 meshY <- rep(y, n)
 meshY <- sort(meshY)
 
 pred.values <- matrix(nrow = n^2, ncol = ncol(mars.model$ranges.and.medians), 
 data = x, byrow = T)
 
 for(i in 1:(n^2))pred.values[i, predictor1] <- meshX[i]
 for(i in 1:(n^2))pred.values[i, predictor2] <- meshY[i]
 Z <- predict.polymars(mars.model, pred.values,intercept=intercept)[,response]
 Z <- matrix(Z, ncol = n, byrow = F)
 
 
 xtitle<-paste("Predictor", predictor1)
 ytitle<-paste("Predictor", predictor2)
 if(mars.model$responses > 1) 
  {
 
   if(missing(main) && (!contour.ploymars))
    {
     ztitle <- paste("Response", response)
    }
     
   if(missing(main) && (contour.ploymars))
     {
      ztitle <- paste("Contour of response",response)
     }
  }
 else
  {
   if(missing(main) && (!contour.ploymars))ztitle <- "Response"
   if(missing(main) && contour.ploymars)ztitle <- paste("Contour of response")
  }
 if(!contour.ploymars)
  {
   image(X, y, Z, xlab = xtitle, ylab= ytitle, zlab = ztitle)
  }
 else
  {
   
     contour(X, y, Z)
   
   }
 invisible()
}
################################################################################################

