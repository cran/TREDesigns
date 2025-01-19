#' Partially Balanced Ternary Residual Effect Designs for Number of Treatments, v ( >= 5)
#'
#' @param v Number of treatments (>= 5)
#'@description
#'This function generates a new class of partially balanced ternary residual effect designs (PBTREDs) for all treatments, v (>= 5) in periods, p = (v+3)/2 if v is odd, otherwise p = (v+2)/2, and the number of sequences, n = 2v.
#'
#' @return It returns a new class of PBTREDs along with its parameters, Information matrix (C), Average Variance Factor (AVF), and Canonical Efficiency Factor (CEF) for both treatment and residual effects.
#' @export
#'
#' @examples
#' library(TREDesigns)
#' PBtRED3(v = 5)
PBtRED3<-function(v){
  if(v<5){
    return(message("Please enter v (>= 5)."))
  }
  expand<-function(vec,times){
    final<-NULL
    for(i in 0:(times)){
      final<-rbind(final,vec+i)
    }
    final<-final%%(times+1)
    final[final==0]<-(times+1)
    return(final)
  }
  #######williams method
  if(v%%2==0){
    sequences<-matrix(NA,nrow=(v/2)+1,ncol=2)
    vv=ceiling(v/4)
  }else{
    sequences<-matrix(NA,nrow=((v+1)/2)+1,ncol=2)
    vv=ceiling((v+1)/4)
  }

  for(j in 1:(vv)){
    elements1<-c(1+j-1,v-j+1)
    elements2<-c(v-j+1,1+j-1)
    sequences[((2*j)-1):(2*j),1]<-elements1
    sequences[((2*j)-1):(2*j),2]<-elements2
  }

  sequences<-sequences[-nrow(sequences),]
  #################expand the initial blocks
  initial_block_elements<-sequences
  des<-NULL
  for(i in 1:ncol(initial_block_elements)){
    des<-cbind(des,(t(expand(initial_block_elements[,i],v-1))))
  }
  pbtcod4<-rbind(des[nrow(des),],des)
  if(v%%2==0){
    p=(v+2)/2
  }else{
    p=(v+3)/2
  }
  lm<-c(list("PBTRED"=pbtcod4,"v"=v,"p"=p,"n"=2*v,"r"=2*p,"t"=2),Study_RED(pbtcod4))
  return(lm)
}
