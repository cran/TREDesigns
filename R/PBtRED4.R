#' Partially Balanced Ternary Residual Effect Designs for Number of Treatments, v = 2m (m >= 4 & even)
#'
#' @param m Any even number (>= 4)
#'@description
#'This function generates a new class of partially balanced ternary residual effect designs (PBTREDs) for the number of treatments, v = 2m; where m >= 4 & even. The number of periods, p = (v+4)/2, while the number of sequences, n = v.
#' @return It returns a new class of PBTREDs along with its parameters, Information matrix (C), Average Variance Factor (AVF), and Canonical Efficiency Factor (CEF) for both treatment and residual effects.
#' @export
#'
#' @examples
#' library(TREDesigns)
#' PBtRED4(m = 4)
PBtRED4<-function(m){
  if(m%%2!=0||m<4){
    return(message("Please enter m >= 4 & even; where m = v/2."))
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
  v=2*m ##(m>=4 & even)
  base_array<-matrix(1:v,nrow=2,byrow = TRUE)
  #######williams method
  sequences<-matrix(NA,nrow=ncol(base_array),ncol=2)
  for(i in 1:nrow(base_array)){
    first<-base_array[i,1]
    last<-base_array[i,ncol(base_array)]
    for(j in 1:((ncol(base_array))/2)){
      elements<-c(first+j-1,last-j+1)
      sequences[((2*j)-1):(2*j),i]<-elements
    }
  }
  #################expand the initial blocks
  initial_block_elements<-sequences
  des<-list()
  for(i in 1:ncol(initial_block_elements)){
    des<-append(des,list(t(expand(initial_block_elements[,i],m-1))))
  }
  des[[2]]<-des[[2]]+m
  des<-do.call(cbind,des)
  first_row<-c(des[nrow(des),(m+1):(2*m)],des[nrow(des),1:m])
  pbtcod3<-rbind(first_row,des,first_row)
  row.names(pbtcod3)<-NULL
  lm<-c(list("PBTRED"=pbtcod3,"v"=v,"p"=m+2,"n"=v,"r"=(m+2)),Study_RED(pbtcod3))
  return(lm)
}
