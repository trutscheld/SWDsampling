#' Design matrix for SWD model 
#'
#' @description Creates a design matrix for stepped wedge design (SWD) given the design parameters 
#' @param nC Number of clusters
#' @param nT Number of timepoint
#' @param nSw number of clusters switches from control to treatment at each timepoint
#' @return design matrix of size noCluster x noTimepoints
#' @examples
#' designMatrix.SWD(5,6,1)
#' designMatrix.SWD(10,6,2)
#' @export
designMatrix.SWD<-function(nC, nT, nSw){
  
  ma<-sapply(1:nT, function(i){
    
    noTr<-(i-1)*nSw
    noC<-nC-noTr
    return(c(rep(1,noTr), rep(0,noC)))
  })
  return(ma)
  
}


#' Design matrix for SWD model under an intervention effectiveness pattern 
#'
#' @description Creates a implementation matrix for a given stepped wedge design and intervention effectiveness pattern
#' @param nC Number of clusters
#' @param nT Number of timepoint
#' @param nSw number of clusters switches from control to treatment at each timepoint
#' @param pattern a vector for intervention effectiveness pattern gives the derivation from 100 percent effectiveness over time
#' @return intervention effectivness pattern matrix as design matrix of size noCluster x noTimepoints
#' @examples
#' implemMatrix.SWD(5,6,1, c(seq(0.4,1,0.2),1))
#' implemMatrix.SWD(10,6,2, c(seq(1,0.4,-0.2),1))
#' @export
implemMatrix.SWD<-function(nC, nT, nSw, pattern){
  
  if(length(pattern)+1==nT){
    ma<-sapply(1:nT, function(i){
      
      noTr<-(i-1)*nSw
      noC<-nC-noTr      
      vec<-c(if( i>1 ){rep(pattern[(i-1):1],each=noTr/(i-1))} , 
             rep(0,noC) )
      return(vec)
    })
    return(ma) 
  }else{
    stop("T he length of the pattern must be one less than the number of timepoints.")
  }
 
}
