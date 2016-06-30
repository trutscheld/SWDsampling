
#' Simulate data and estimate lmm within stepped wedge design model
#'
#' @description For a given one scenario build mean intensities of SWD model, sample data and estimate treatment parameter using a linear mixed model
#' TP number of timepoints, I number of cluster. The design matrix has to be coded by zeros and ones.
#' @param I number of clusters (design parameter)
#' @param TP number of timepoints (design parameter)
#' @param mu baseline mean (model parameter)
#' @param theta treatment effect (model parameter)
#' @param beta.j vector of time trents (model parameter)
#' @param sigma.alpha between cluster variability as standard deviation (model parameter)
#' @param X.i.j.0 assumed treatment model matrix for a SWD study (model parameter)
#' @param N number of individuals (fixed) for all clusters and timepoints
#' @param sigma.e random error variability as standard deviation (model parameter)
#' @param sigma.ind individual variability as standard deviation (model parameter), if it is an longitudinal model, by default (NULL) it is an cross-sectional model
#' @param A derivation from perfect 100 percent effectiveness pattern (simulation parameter)
#' @param B timepoint of cluster loss (simulation parameter) with 4 possibilities: "0": default - no cluster at no timepoint get lost, "1" - Cluster missing at random from timepoint 2 untill TP, "2" - Cluster is missing at beginning (1/3 of timepoints after the first), "3" - Cluster is missing at end (1/3 of the last timepoints). 
#' @param C number of cluster loss (simulation parameter), by default zero. If a cluster get lost from time point i, all indiviual responses of that cluster will be deleted from timepoint i until timpeoint TP (end).
#' @param D number of individuals loss (simulation parameter), by default zero. If not zero, then individual responses to delete are selected at random from timepoints and clusters.
#' @return treatment and time effects estimated with linear mixed model  
#' #@examples
#' @export
simABC<-function(I,TP, mu,theta,beta.j,sigma.alpha,  X.i.j.0, N,sigma.e, sigma.ind=NULL, A=NULL,B="0",C=0, D=0){
  
  #data sampling given model and data deletion of the simulation
  data.delete<- SWD.datasampling(I,TP, mu,theta,beta.j,sigma.alpha,  X.i.j.0, N,sigma.e,sigma.ind, A, B, C, D)      
  
  ####################    Estimation by linear mixed model ########################
  #timpoint 1 <- 0
  data.delete$time<-as.factor((as.numeric(data.delete$time)-1))
  
  if(is.null(sigma.ind)){#cross-sectional model
    m0.i <- lmer(values ~ 1 + xt.0 + time +  (1 | cluster),data = data.delete)
    ##nullmodel, if there is no treatment
    mbase.i <- update(m0.i, values ~ 1 + time +  (1 | cluster))
  }else{#longitudinal model
    m0.i <- lmer(values ~ 1 + xt.0 + time +  (1 | cluster)+  (1 | id),data = data.delete)
    ##nullmodel, if there is no treatment
    mbase.i <- update(m0.i, values ~ 1 + time +  (1 | cluster)+  (1 | id))
  }
  
  summ<-summary(m0.i)
  #print(summ)
  #return(coef(summ))
  #return(coef(summ)[c("xt.0",paste("time",1:(TP-1),sep="")),"Estimate"])
  return(c(coef(summ)[c("xt.0",paste("time",1:(TP-1),sep="")),"Estimate"], 
           SE=coef(summ)[c("xt.0",paste("time",1:(TP-1),sep="")), "Std. Error"],
           anova.p=anova(m0.i,mbase.i)[8][2,1]))
  
}

#' Simulate data and estimate lmm within stepped wedge design model
#'
#' @description For a given one scenario build mean intensities of SWD model, sample data and estimate treatment parameter using a linear mixed model
#' TP number of timepoints, I number of cluster. The design matrix has to be coded by zeros and ones.
#' @param I number of clusters (design parameter)
#' @param TP number of timepoints (design parameter)
#' @param mu baseline mean (model parameter)
#' @param theta treatment effect (model parameter)
#' @param beta.j vector of time trents (model parameter)
#' @param sigma.alpha between cluster variability as standard deviation (model parameter)
#' @param X.i.j.0 assumed treatment model matrix for a SWD study (model parameter)
#' @param N number of individuals (fixed) for all clusters and timepoints
#' @param sigma.e random error variability as standard deviation (model parameter)
#' @param sigma.ind individual variability as standard deviation (model parameter), if it is an longitudinal model, by default (NULL) it is an cross-sectional model
#' @param A derivation from perfect 100 percent effectiveness pattern (simulation parameter)
#' @param B timepoint of cluster loss (simulation parameter) with 4 possibilities: "0": default - no cluster at no timepoint get lost, "1" - Cluster missing at random from timepoint 2 untill TP, "2" - Cluster is missing at beginning (1/3 of timepoints after the first), "3" - Cluster is missing at end (1/3 of the last timepoints). 
#' @param C number of cluster loss (simulation parameter), by default zero. If a cluster get lost from time point i, all indiviual responses of that cluster will be deleted from timepoint i until timpeoint TP (end).
#' @param D number of individuals loss (simulation parameter), by default zero. If not zero, then individual responses to delete are selected at random from timepoints and clusters.
#' @return linear mixed model  
#' #@examples
#' @export
simABC.model<-function(I,TP, mu,theta,beta.j,sigma.alpha,  X.i.j.0, N,sigma.e, sigma.ind=NULL, A=NULL,B="0",C=0, D=0){
  
  #data sampling given model and data deletion of the simulation
  data.delete<- SWD.datasampling(I,TP, mu,theta,beta.j,sigma.alpha,  X.i.j.0, N,sigma.e,sigma.ind, A, B, C, D)      
  
  ####################    Estimation by linear mixed model ########################
  #timpoint 1 <- 0
  data.delete$time<-as.factor((as.numeric(data.delete$time)-1))
  
  if(is.null(sigma.ind)){#cross-sectional model
    m0.i <- lmer(values ~ 1 + xt.0 + time +  (1 | cluster),data = data.delete)
    ##nullmodel, if there is no treatment
    mbase.i <- update(m0.i, values ~ 1 + time +  (1 | cluster))
  }else{#longitudinal model
    m0.i <- lmer(values ~ 1 + xt.0 + time +  (1 | cluster)+  (1 | id),data = data.delete)
    ##nullmodel, if there is no treatment
    mbase.i <- update(m0.i, values ~ 1 + time +  (1 | cluster)+  (1 | id))
  }
  
  summ<-summary(m0.i)
  #print(summ)
  #return(coef(summ))
  #return(coef(summ)[c("xt.0",paste("time",1:(TP-1),sep="")),"Estimate"])
  return(m0.i)
  
}
