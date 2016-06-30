#' Stepped wedge design model 
#'
#' @description Performs the mean values of a given Stepped wedge design (SWD)
#' @param I number of clusters (design parameter)
#' @param TP number of timepoints (design parameter)
#' @param mu baseline mean (model parameter)
#' @param theta treatment effect (model parameter)
#' @param beta.j vector of time trents (model parameter)
#' @param sigma.alpha between cluster variability as standard deviation (model parameter)
#' @param X.i.j.0 assumed treatment model matrix for a SWD study (model parameter)
#' @param X.i.j data model matrix of real intervention implementation (model parameter), by default (NULL) the same as the X.i.j.0
#' @return Data frame with mean intensities corresponding to the SWD model and full model parameter information
#' @examples
#' noCl<-10
#' noT<-6
#' switches<-2
#' DM<-designMatrix.SWD(noCl,noT,switches)
#' buildmodeldesign(I=noCl,TP=noT, mu=0,theta=1,beta.j=rep(1,noT),sigma.alpha=0.5, X.i.j.0=DM)
#' @export
buildmodeldesign<-function(I,TP, mu,theta,beta.j,sigma.alpha, X.i.j.0, X.i.j=NULL){
  
  if(is.null(X.i.j)){X.i.j<-X.i.j.0}
  
  if(length(beta.j)!=TP){
    stop("number of time point effects not equal to number of time points")
  }else{
  # Cluster effects
  alpha.i <- rnorm(I,mean = 0,sd = sigma.alpha)  
  
  # Mean of cluster i at time j 
  mu.i.j <- rep(0, I*TP)
  dim(mu.i.j) <- c(I, TP)
  
  for(i in 1:I){
    for(j in 1:TP){
      # Mean of cluster i at time j 
      mu.i.j[i, j] <- mu + alpha.i[i] + beta.j[j] + theta * X.i.j[i, j] 
    }
  }
  
  
  # Names TPime points
  colnames(mu.i.j) <- paste("t",1:TP, sep = ".")
  colnames(X.i.j) <- paste("xt",1:TP, sep = ".")
  colnames(X.i.j.0) <- paste("xt",1:TP, sep = ".")
  
  # Name of clusters
  #cluster  <- factor(paste("cluster", 1:I, sep = "."))
  cluster  <- 1:I
  
  # Convert to data frame
  mu.i.j <- data.frame(mu.i.j, cluster=cluster)
  x.i.j <- data.frame(X.i.j, cluster = cluster)
  x.i.j.0 <- data.frame(X.i.j.0, cluster = cluster)
  
  tmp1 <- reshape(mu.i.j, timevar = "time",
                  varying = paste("t",1:TP, sep = "."), 
                  idvar = "cluster", 
                  direction = "long"          
  )
  names(tmp1) <- c("cluster", "time", "mu")
  
  # Convert design matrix of data to data frame
  tmp2 <- reshape(x.i.j, timevar = "time",
                  varying = paste("xt",1:TP, sep = "."), 
                  idvar = "cluster", 
                  direction = "long"          
  )
  
  # Convert design matrix to data frame
  tmp2.0 <- reshape(x.i.j.0, timevar = "time",
                    varying = paste("xt",1:TP, sep = "."), 
                    idvar = "cluster", 
                    direction = "long"          
  )
  
  # Data frame with full information at the level of the cluster
  dat <-  cbind(tmp1, xt = tmp2[, "xt"], xt.0 = tmp2.0[, "xt"])
  names(dat) <- c("cluster", "time", "mu", "xt", "xt.0")
  
  return(dat)
  
  }#end ifelse number beta.j
}

#' Sampling Response of individuals within a SWD model
#'
#' @description Sample data (response) for given numbers of individuals by given group means within a SWD model 
#' @param I number of clusters (design parameter)
#' @param TP number of timepoints (design parameter)
#' @param mu baseline mean (model parameter)
#' @param theta treatment effect (model parameter)
#' @param beta.j vector of time trents (model parameter)
#' @param sigma.alpha between cluster variability as standard deviation (model parameter)
#' @param X.i.j.0 assumed treatment model matrix for a SWD study (model parameter)
#' @param X.i.j data model matrix of real intervention implementation (model parameter), by default (NULL) the same as the X.i.j.0
#' @param N number of individuals (fixed) for all clusters and timepoints
#' @param sigma.e random error variability as standard deviation (model parameter)
#' @param sigma.ind individual variability as standard deviation (model parameter), if it is an longitudinal model, by default (NULL) it is an cross-sectional model
#' @return Data frame with individuals intensities corresponds to the SWD model and full model parameter information
#' @examples
#' noCl<-10
#' noT<-6
#' switches<-2
#' DM<-designMatrix.SWD(noCl,noT,switches)
#' sampleData(I=noCl,TP=noT, mu=0,theta=1,beta.j=rep(1,noT),sigma.alpha=0.5, X.i.j.0=DM,N=10,sigma.e=1)
#' @export
sampleData<-function(I,TP, mu,theta,beta.j,sigma.alpha, X.i.j.0, X.i.j=NULL,N,sigma.e, sigma.ind=NULL){  
  
  #create SWD model given the design parameters
  model<-buildmodeldesign(I,TP, mu,theta,beta.j,sigma.alpha, X.i.j.0)
  
  #N individuals pro I Cluster haben innersubjectvarianz über alle TP Timpoints hinweg
  if(!is.null(sigma.ind)){
    #print("indiv var")
    error.inner.i.n<-replicate(N,with(model,rnorm(I, mean = 0, sd = sigma.ind)[as.numeric(cluster)]))
    #print(error.inner.i.n)
  }
  
  y <- rep(0, N*TP*I)
  dim(y) <- c(TP*I, N)
  blocks <- paste("block", 1:N, sep=".")
  colnames(y) <- blocks
  
  for(i in 1:N){y[, i] <- with(model, rnorm(TP*I, mean = mu, sd = rep(sigma.e, length(mu))))  }
  
  if(!is.null(sigma.ind)){y<-y+error.inner.i.n}
  y <- data.frame(y)
  y <- stack(y)
  
  datall <- NULL
  for(i in 1:N)
  {
    tm <- cbind(model, y[y$ind == blocks[i], ])  
    datall <- rbind(datall, tm)
  }
  
  datall$id <- with(datall, paste(cluster, ind, sep = "."))
  datall$id <- factor(datall$id)
  datall$time <- factor(datall$time)
  
  return(datall)
  
}

#' For given information of loss sampling response of individuals within a SWD model 
#'
#' @description Data for given numbers of individuals by given group means within a SWD model and derivations with loss of data
#' @param data.all Sampled data (response) for given numbers of individuals by given group means within a SWD model 
#' @param I number of clusters (design parameter)
#' @param TP number of timepoints (design parameter)
#' @param B timepoint of cluster loss with 4 possibilities: "0": default - no cluster at no timepoint get lost, "1" - Cluster missing at random from timepoint 2 untill TP, "2" - Cluster is missing at beginning (1/3 of timepoints after the first), "3" - Cluster is missing at end (1/3 of the last timepoints). 
#' @param C number of cluster loss, by default zero. If a cluster get lost from time point i, all indiviual responses of that cluster will be deleted from timepoint i until timpeoint TP (end).
#' @param D number of individuals loss, by default zero. If not zero, then individual responses to delete are selected at random from timepoints and clusters.
#' @return Data frame with individuals intensities corresponds to the SWD model and full model parameter information and derivation information
#' @examples
#' noCl<-10
#' noT<-6
#' switches<-2
#' DM<-designMatrix.SWD(noCl,noT,switches)
#' #cross-sectional SWD (10 cluster and 6 time points)
#' #no derivation from perfect 100 percent effectiveness pattern
#' #no data loss (no missing)
#' data<-SWD.datasampling(I=noCl,TP=noT, mu=0,theta=1,beta.j=rep(1,noT),sigma.alpha=0.5, X.i.j.0=DM,N=10,sigma.e=1)
#' #no missing in data
#' Data.loss.SWD(data.all=data, I=noCl,TP=noT)
#' #missing individuals
#' Data.loss.SWD(data.all=data, I=noCl,TP=noT, D=5)
#' #missing 2 cluster at random
#' Data.loss.SWD(data.all=data, I=noCl,TP=noT, B="1", C=2)
#' @export
Data.loss.SWD<-function(data.all, I,TP,B="0",C=0, D=0){
  
  ########  data deletion depends on:  B,C, D      ##############
  
  if(C==0){data.delete<-data.all
  }else{  
    
    #number of cluster loss
    clusternr<-sample(1:I, size=C, replace=FALSE)
    
    #timepoint of cluster loss
    switch(B,
           #B0: no Cluster missing
           "0"={timepoint<-0},
           #B1: Cluster missing at random
           "1"={timepoint<-sample(2:TP, size=C, replace=TRUE)  },
           #B2: Cluster is missing at beginning
           "2"={timepoint<-sample((1:round(TP/4))+1, size=C, replace=TRUE)  },
           #B3: Cluster is missing at end
           "3"={timepoint<-sample((TP-round(TP/4)+1):TP, size=C, replace=TRUE)  }           
    )
    
    dc<-NULL
    for(i in 1:C){
      
      which.tp<-timepoint[i]:TP
      dc<-rbind(dc, cbind(rep(clusternr[i],length(which.tp)), timepoint[i]:TP))
    }  
    colnames(dc)<-c("cluster", "time")
    
    ###given vector of  timepoints and corresponding vector of cluster
    #collect rows whcih has to delete
    collect.all<-NULL
    for(i in 1:dim(dc)[1]){   
      collect<-unlist(sapply(1:dim(data.all)[1], function(j){
        
        if(identical(unname(unlist(data.all[j,c("cluster", "time")]),force=TRUE),unname(unlist(dc[i,]),force=TRUE))){return(j)}
      }))
      collect.all<-c(collect.all,collect )
    }
    #now delete rows in data 
    data.delete<-data.all[-collect.all, ] 
    
  }
  if(D!=0){#additionally individiual loss
    
    loss.ind<-sample(1:dim(data.delete)[1], size=D, replace=FALSE)
    data.delete<-data.delete[-loss.ind, ] 
  }
  
  return(data.delete)
}  

#' Sampling Response of individuals within a SWD model
#'
#' @description Sample data (response) for given numbers of individuals by given group means within a SWD model and derivations
#' @param I number of clusters (design parameter)
#' @param TP number of timepoints (design parameter)
#' @param mu baseline mean (model parameter)
#' @param theta treatment effect (model parameter)
#' @param beta.j vector of time trents (model parameter)
#' @param sigma.alpha between cluster variability as standard deviation (model parameter)
#' @param X.i.j.0 assumed treatment model matrix for a SWD study (model parameter)
#' @param X.i.j data model matrix of real intervention implementation (model parameter), by default (NULL) the same as the X.i.j.0
#' @param N number of individuals (fixed) for all clusters and timepoints
#' @param sigma.e random error variability as standard deviation (model parameter)
#' @param sigma.ind individual variability as standard deviation (model parameter), if it is an longitudinal model, by default (NULL) it is an cross-sectional model
#' @param A derivation from perfect 100 percent effectiveness pattern
#' @param B timepoint of cluster loss with 4 possibilities: "0": default - no cluster at no timepoint get lost, "1" - Cluster missing at random from timepoint 2 untill TP, "2" - Cluster is missing at beginning (1/3 of timepoints after the first), "3" - Cluster is missing at end (1/3 of the last timepoints). 
#' @param C number of cluster loss, by default zero. If a cluster get lost from time point i, all indiviual responses of that cluster will be deleted from timepoint i until timpeoint TP (end).
#' @param D number of individuals loss, by default zero. If not zero, then individual responses to delete are selected at random from timepoints and clusters.

#' @return Data frame with individuals intensities corresponds to the SWD model and full model parameter information and derivation information
#' @examples
#' noCl<-10
#' noT<-6
#' switches<-2
#' DM<-designMatrix.SWD(noCl,noT,switches)
#' #cross-sectional SWD (10 cluster and 6 time points)
#' #no derivation from perfect 100 percent effectiveness pattern
#' #no data loss (no missing)
#' #SWD.datasampling(I=noCl,TP=noT, mu=0,theta=1,beta.j=rep(1,noT),sigma.alpha=0.5, X.i.j.0=DM,N=10,sigma.e=1)
#' #cross-sectional SWD (10 cluster and 6 time points)
#' #no derivation from perfect 100 percent effectiveness pattern
#' #missing individuals
#' #SWD.datasampling(I=noCl,TP=noT, mu=0,theta=1,beta.j=rep(1,noT),sigma.alpha=0.5, X.i.j.0=DM,N=10,sigma.e=1, D=5)
#' #cross-sectional SWD (10 cluster and 6 time points)
#' #no derivation from perfect 100 percent effectiveness pattern
#' #missing 2 cluster at random
#' #SWD.datasampling(I=noCl,TP=noT, mu=0,theta=1,beta.j=rep(1,noT),sigma.alpha=0.5, X.i.j.0=DM,N=10,sigma.e=1 ,B="1", C=2)
#' #longitudinal SWD (10 cluster and 6 time points)
#' #no derivation from perfect 100 percent effectiveness pattern
#' #no data loss (no missing)
#' #SWD.datasampling(I=noCl,TP=noT, mu=0,theta=1,beta.j=rep(1,noT),sigma.alpha=0.5, X.i.j.0=DM,N=10,sigma.e=1, sigma.ind=0.5)
#' @export
SWD.datasampling<-function(I,TP, mu,theta,beta.j,sigma.alpha, X.i.j.0, N,sigma.e, sigma.ind=NULL, A=NULL, B="0", C=0, D=0){
  
  ########  model matrix depends on         ############## 
  #######   Factor A: Learning Effect       ##############
  ########################################################
  
  if (is.null(A)){
    X.i.j.B <-X.i.j.0
  }else{
    X.i.j.B <-A
  }
  
  if(length(beta.j)!=TP){
    stop("number of time point effects not equal to number of time points")
  }else{
    
    
    #sample data for given individuals parameter within a model gien the model parameters
    data.all<-sampleData(I,TP, mu,theta,beta.j,sigma.alpha, X.i.j.0, X.i.j,N,sigma.e, sigma.ind)
    #data.delete<-data.all
    ## if given data deletion
    data.delete<-Data.loss.SWD(data.all, I,TP,B,C, D)
    
    return(data.delete)
  }#end ifelse error
}

