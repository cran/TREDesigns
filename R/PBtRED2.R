#' Partially Balanced Ternary Residual Effect Designs for a Prime Number of Treatments, v (>=5); Series 1: v = 4m + 1, Series 2: v = 4m + 3 and Series 3: v = 6m + 1
#'
#' @param m Any number (>= 1) such that v (>=5) is prime.
#' @param series Choose series: 1, 2 or 3
#'@description
#'This function generates three series of partially balanced ternary residual effect designs (PBTREDs) for a prime number of treatments (v), where (v >= 5) with p periods, and n sequences.
#'
#'Parameters of Series 1: v = 4m+ 1, t = m, n = m(4m+ 1), p = 5, r = 5m
#'
#'Parameters of Series 2: v = 4m+ 3, t = 2, n = 2(4m+ 3), p = 2(m+1), r = 4(m+1)
#'
#'Parameters of Series 3: v = 6m+ 1, t = m, n = m(6m+1), p = 7, r =7m
#'
#' @return It returns the a new class of PBTREDs based on chosen m along with its parameters, Information matrix (C), Average Variance Factor (AVF), and Canonical Efficiency Factor (CEF) for both treatment and residual effects.
#' @export
#'
#' @examples
#' library(TREDesigns)
#' PBtRED2(m = 1, series=1)
PBtRED2<-function(m,series){
  #############
  if(series==1){
    v=4*m+1
  }
  if(series==2){
    v=4*m+3
  }
  if(series==3){
    v=6*m+1
  }
  ###############

  ########
  # Function to check if a number is prime
  is_prime <- function(n) {
    # Handle edge cases
    if (n <= 1) {
      return(FALSE) # 1 and numbers <= 0 are not prime
    }
    if (n == 2) {
      return(TRUE) # 2 is the only even prime number
    }
    if (n %% 2 == 0) {
      return(FALSE) # Exclude other even numbers
    }

    # Check divisors from 3 to square root of n
    max_divisor <- floor(sqrt(n))
    if (max_divisor >= 3) {  # Only run loop if max_divisor >= 3
      for (i in seq(3, max_divisor, by = 2)) {
        if (n %% i == 0) {
          return(FALSE)
        }
      }
    }

    return(TRUE) # If no divisors are found, it's prime
  }

  #######

  if(is_prime(v)==FALSE){
    return(message("Please enter m such that v is prime."))
  }
  ################
  # List of prime numbers up to 100
  primes <- c(2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97)

  # Define a vector to hold the primitive elements
  primitive_elements <- c()

  # Assign primitive elements for each prime based on given rules and extending
  for (s in primes) {
    if (s == 3 || s == 7 || s == 17 || s == 29 || s == 31 || s == 37) {
      pr <- 3
    } else if (s == 23) {
      pr <- 5
    } else {
      pr <- 2
    }
    primitive_elements <- c(primitive_elements, pr)
  }
  # Combine primes and their corresponding primitive elements into a data frame
  result <- data.frame(Prime = primes, PrimitiveElement = primitive_elements)
  pr=result[which(result[,1]==v),2]

  ###########series 1
  expand<-function(vec,times){
    final<-NULL
    for(i in 0:(times)){
      final<-rbind(final,vec+i)
    }
    final<-final%%(times+1)
    final[final==0]<-(times+1)
    return(final)
  }
  #expand(initial_blocks[,1],12)
  if(series==1){
    #v=4*m+1 ## v should be a prime number
    initial_blocks_of_powers<-matrix(c(0:(m-1)),1,m,byrow = TRUE)
    for(i in 2:4){
      initial_blocks_of_powers<-rbind(initial_blocks_of_powers,initial_blocks_of_powers[nrow(initial_blocks_of_powers),]+m)
    }
    #initial_blocks<-initial_blocks[1:4,]  #series 1 always in 4 rows
    initial_block_elements<-matrix(NA,nrow=4,ncol=m)
    for(i in 1:nrow(initial_block_elements)){
      for(j in 1:ncol(initial_block_elements)){
        initial_block_elements[i,j]<-pr^initial_blocks_of_powers[i,j]
      }
    }
    initial_block_elements<-initial_block_elements%%v
    initial_block_elements[initial_block_elements==0]<-v
    #################expand the initial blocks
    des<-NULL
    for(i in 1:ncol(initial_block_elements)){
      des<-cbind(des,t(expand(initial_block_elements[,i],v-1)))
    }
    des_series1<-rbind(des[nrow(des),],des)
    lm<-c(list("PBTRED"=des_series1,"v"=v,"p"=5,"n"=v*m,"r"=nrow(des_series1)*m,"t"=m),Study_tRED(des_series1))
    return(lm)
  }

  #####series 2
  if(series==2){
    #v=4*m+3  ## v should be a prime number
    seq1<-seq(0,(4*m),by=2)
    seq2<-seq(1,((4*m)+1),by=2)
    initial_blocks_of_powers<-(cbind(seq1,seq2))
    #########
    initial_block_elements<-matrix(NA,nrow=length(seq1),ncol=2)
    for(i in 1:nrow(initial_block_elements)){
      for(j in 1:ncol(initial_block_elements)){
        initial_block_elements[i,j]<-pr^initial_blocks_of_powers[i,j]
      }
    }
    initial_block_elements<-initial_block_elements%%v
    initial_block_elements[initial_block_elements==0]<-v
    #################expand the initial blocks
    des<-NULL
    for(i in 1:ncol(initial_block_elements)){
      des<-cbind(des,t(expand(initial_block_elements[,i],v-1)))
    }
    des_series2<-rbind(des[nrow(des),],des)
    lm<-c(list("PBTRED"=des_series2,"v"=v,"p"=2*(m+1),"n"=v*2,"r"=4*(m+1),"t"=2),Study_tRED(des_series2))
    return(lm)
  }

  ###series 3
  if(series==3){
    #v=(6*m)+1
    #v=4*m+1 ## v should be a prime number
    initial_blocks_of_powers<-matrix(c(0:(m-1)),1,m,byrow = TRUE)
    for(i in 2:6){
      initial_blocks_of_powers<-rbind(initial_blocks_of_powers,initial_blocks_of_powers[nrow(initial_blocks_of_powers),]+m)
    }
    initial_block_elements<-matrix(NA,nrow=6,ncol=m)
    for(i in 1:nrow(initial_block_elements)){
      for(j in 1:ncol(initial_block_elements)){
        initial_block_elements[i,j]<-pr^initial_blocks_of_powers[i,j]
      }
    }
    initial_block_elements<-initial_block_elements%%v
    initial_block_elements[initial_block_elements==0]<-v
    #################expand the initial blocks
    des<-NULL
    for(i in 1:ncol(initial_block_elements)){
      des<-cbind(des,t(expand(initial_block_elements[,i],v-1)))
    }
    des_series3<-rbind(des[nrow(des),],des)
    lm<-c(list("PBTRED"=des_series3,"v"=v,"p"=7,"n"=v*m,"r"=7*(m),"t"=m),Study_tRED(des_series3))
    return(lm)
  }
}
