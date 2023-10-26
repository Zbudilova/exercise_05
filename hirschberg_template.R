library('seqinr')
library("Biostrings")


#' Align two sequences globally using Hirschberg's algorithm
#' 
#' @param X DNAString object representing NT or AA sequence to align
#' @param Y DNAString object representing NT or AA sequence to align
#' @param align A list of DNAString objects with alignment of input sequences
#' @param match An integer value of a score for matching bases
#' @param mismatch An integer value of score for mismatching bases
#' @param gap An integer value of penalty for gap insertion
HirschbergTemplate <- function(X, Y, align, match, mismatch, gap){
    
    first_align_row <- align[[1]] # initialize the first row of alignment
    second_align_row <- align[[2]] # initialize the second row of alignment
  
  
    if (length(X)==0) # length of X is equal to zero
    {
        for (i in (1:length(Y)))# for each character in Y
        {
            first_align_row <- xscat(first_align_row, "-") # add gap
            second_align_row <- xscat(second_align_row, Y[i]) # add character from Y
        }
        align <- c(DNAStringSet(first_align_row), DNAStringSet(second_align_row))
        print(align)
    }
    else if  (length(Y) == 0)# length of Y is equal to zero
    {
        for (i in (1:length(X)))# for each character in X
        {
            first_align_row <- xscat(first_align_row, X[i])# add character from X
            second_align_row <- xscat(second_align_row, "-") # add gap
        }
        align <- c(DNAStringSet(first_align_row), DNAStringSet(second_align_row))
        print(align)
    }
    else if ((length(X)==1)||(length(Y)==1)) # length of X and Y is equal to 1
    {
        first_align_row <- xscat(first_align_row, X[1]) # add character from X
        second_align_row <- xscat(second_align_row, Y[1])# add character from Y
        align <- c(DNAStringSet(first_align_row), DNAStringSet(second_align_row))
        print(align)
    }
    else
    {
        xlen <- length(X) # length of X
        xmid <- xlen/2 # half of the length of X
        ylen <- length(Y) # length of Y
        
        left_score <- NWScore(X[1:xmid], Y, match, mismatch, gap)# NW score for the first half of X and the whole Y
        right_score <- NWScore(reverse(X[(xmid+1):xlen]), reverse(Y), match, mismatch, gap)# NW score for the second half of X and the whole Y (both are reversed)
        suma <- left_score+rev(right_score)
        ymid <- match(max(suma), suma) - 1

        # The first half
        if (ymid == 0)# index of division for Y is equal to 0
        {
            align <- HirschbergTemplate(X[1:(length(X)/2)],DNAString(""),align,match, mismatch, gap) # call Hirschberg function for the first half of X and for an empty DNAString object
        }
        else
        {
            align <- HirschbergTemplate(X[1:(length(X)/2)],Y[1:ymid],align,match, mismatch, gap) # call Hirschberg function for the first half of X and for the first part of Y
        }
        
        # The second half
        if ((xmid + 1) > xlen) # X cannot be further divided
        {
            align <- HirschbergTemplate(DNAString(""),Y[(ymid+1):length(Y)],align,match, mismatch, gap)# call Hirschberg function for an empty DNAString object and the second half of Y
        }
        else if ((ymid + 1) > ylen) # Y cannot be further divided
        {
            align <- HirschbergTemplate(X[((length(X)/2)+1):length(X)],DNAString(""),align,match, mismatch, gap)# call Hirschberg function for the second half of X and for an empty DNAString object
        }
        else 
        {
            align <- HirschbergTemplate(X[((length(X)/2)+1):length(X)],Y[(ymid+1):length(Y)],align,match, mismatch, gap)# call hirschberg function for the second half of X and the second part of Y
        }
    }
    return(align)
}


NWScore <- function(X,Y,match,mismatch,gap){
  m <- length(X)
  n <- length(Y)
  S <- (0:n)*gap
  for (i in 2:(m+1)){
    s <- S[1]
    c <- S[1]+gap
    S[1] <- c
    for (j in 2:(n+1)){
      if (X[i-1]==Y[j-1]){
        score <-match
      } else {
        score <- mismatch 
      }
      c <- max(S[j]+gap,c+gap,s+score)
      s <- S[j]
      S[j] <- c
    }
  }
  return(S)
}



X <- DNAString("AGTACGCA")
Y <- DNAString("TATGC")
align <- list(DNAString(""),DNAString(""))
match <- 2
mismatch <- -1
gap <- -2

vysledek <- HirschbergTemplate(X, Y, align, match, mismatch, gap)
vysledek


