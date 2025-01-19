#' Partially Balanced Ternary Residual Effect Designs for Prime Number of Treatments
#'
#' @param v Prime Number of Treatments, v ( >= 5)
#'@description
#'This function generates a class of partially balanced ternary residual effect designs (PBTREDs) for a prime number of treatments (v), where (v >= 5) with p periods, take the value (v+3)/2, and n sequences, take the value v(v-1)/2.
#'
#' @return It returns a new class of PBTREDs along with its parameters, Information Matrix (C), Average Variance Factor (AVF), and Canonical Efficiency Factor (CEF) for both treatment and residual effects.
#' @export
#'
#' @examples
#' library(TREDesigns)
#' PBtRED1(v = 5)
PBtRED1<-function(v){
  s=v
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

  if(is_prime(s)==FALSE || v<5){
    return(message("Please enter a prime number, v(>=5)."))
  }
  #############

  retain<-(v+1)/2
  ############


  list_mols<-MOLS(s)
  retain_mols<-(v-1)/2
  retain_rows<-(v+1)/2
  t=(v-1)/2

  #########
  selected_mols<-list_mols[1:retain_mols]
  #####
  selected_rows_in_mols<-lapply(selected_mols,function(mat)mat[1:retain_rows,])
  #######arrange side by side
  des<-do.call(cbind,selected_rows_in_mols)
  ####
  pbtcod1<-rbind(des[nrow(des),],des)
  ###
  lm<-c(list("PBTRED"=pbtcod1,"v"=v,"p"=nrow(pbtcod1),"n"=v*t,"r"=nrow(pbtcod1)*t,"t"=t),Study_RED(pbtcod1))
  return(lm)

}
