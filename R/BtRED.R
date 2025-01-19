#' Balanced Ternary Residual Effect Designs for Prime Number of Treatments
#'
#' @param v Prime Number of Treatments, v ( >= 5)
#'@description
#'This function generates a class of balanced ternary residual effect designs (BTREDs) for a prime number of treatments (v), where (v >= 5) with p periods, determined as (v+3)/2, and n sequences, take the value v(v-1).
#'
#' @return It returns a new class of BTREDs along with its parameters, Information Matrix (C), Average Variance Factor (AVF), and Canonical Efficiency Factor (CEF) for both treatment and residual effects.
#' @export
#'
#' @examples
#' library(TREDesigns)
#' BtRED(v = 7)

BtRED<-function(v){
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
  ##############

  if(is_prime(v)==FALSE || v<5){
    return(message("Please enter a prime number, v(>=5)."))
  }
  ##########
  retain<-(v+1)/2
  ####generalized MOLS
  MOLS<-function(v){
    MOLS<-list()
    for(i in 1:(v-1)){
      mols<-NULL
      seq<-c(1:v)
      j=0
      repeat{
        mols<-rbind(mols,c(seq+j))
        j=j+i
        if(nrow(mols)==v){
          mols=mols%%v
          mols[mols==0]<-v
          MOLS<-append(MOLS,list(mols))
          break
        }
      }
    }
    return(MOLS)
  }
  ###########
  list_mols<-MOLS(v)
  #dim(array)
  #list_mols<-unlist(apply(array,3,function(i)list(i)),recursive = FALSE)
  retained_row_of_mols<-lapply(list_mols,function(mat)mat[1:retain,])
  ########now write side by side
  des<-do.call(cbind,retained_row_of_mols)
  btcod<-rbind(des[nrow(des),],des)
  #####
  t=v-1
  lm<-c(list("BTRED"=btcod,"v"=v,"p"=nrow(btcod),"n"=ncol(btcod),"r"=nrow(btcod)*t,"t"=t),Study_RED(btcod))
  return(lm)
}

