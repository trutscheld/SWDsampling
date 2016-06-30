#' Power calculation within stepped wedge design model by Hussey et.al 
#'
#' @description Calculation of power for a lmm with cluster as random effect, fixed timepoint effects, but set to null, TP number of timepoints, I number of cluster. The design matrix has to be coded by zeros and ones.
#' @param ThetaEst expected treatment effect
#' @param alpha singificance level (by default 0.05)
#' @param Design design matrix for a given SWD model
#' @param sigmaq within cluster variance
#' @param tauq between cluster variance
#' @param time a logical (FALSE, if no time trends are expected, otherwise TRUE)
#' @return Aproximated power of two tailed test, although the design matrix is fractionated, then power is not valid
#' #@references \cite{Hussey2007}
#' @examples
#' noCl<-10
#' noT<-6
#' switches<-2
#' DM<-designMatrix.SWD(noCl,noT,switches)
#' sigma.e <- 2
#' sigma.alpha <- 2    
#' calcPower.Hussey_SWD(ThetaEst=1,Design=DM, sigmaq=sigma.e^2, tauq=sigma.alpha^2, time=FALSE)
#' calcPower.Hussey_SWD(ThetaEst=1,Design=v, sigmaq=sigma.e^2, tauq=sigma.alpha^2, time=TRUE)
#' @export

calcPower.Hussey_SWD<-function(ThetaEst,alpha=0.05, Design, sigmaq, tauq, time){
  
  #quantile of the standard normal distribution function
  Z<-qnorm(p=1-alpha/2)
  #number of cluster
  I<-dim(Design)[1]
  #number of timepoints
  TP<-dim(Design)[2]
  #number of all clusters and timepoints, which get the treatment
  U<-sum(Design)
  #W:sum over all timepoints: (sum of all clusters with treatment)^2
  W<-sum(colSums(Design)^2)
  #V:sum over all clusters: (sum of all timepoints with treatment)^2
  V<-sum(rowSums(Design)^2)
  
  if(time==TRUE){#with time trend
    #Variance of estimated Theta by lmm estimated with WLS and Designmatrixed encoded with nulls and zeros
    varTheta<-(I*sigmaq*(sigmaq+TP*tauq))/((I*U-W)*sigmaq+(U^2+I*TP*U-TP*W-I*V)*tauq) 
  }
  else{#no time trends
    #Variance of estimated TPheta by lmm estimated with WLS and Designmatrixed encoded with nulls and zeros
    varTheta<-(I*sigmaq*(sigmaq+TP*tauq))/((I*TP*U-U^2)*sigmaq+I*TP*(TP*U-V)*tauq) 
    
  }
  
  return(pnorm(q=ThetaEst/sqrt(varTheta)-Z))
  
  
}
