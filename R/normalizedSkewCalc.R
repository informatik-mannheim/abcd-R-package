# Copyright 2018 by the authors.
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#     http://www.apache.org/licenses/LICENSE-2.0
#  
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License. 

#' @author Markus Gumbel

#' Subsequence from a sequence with a random start position.
#' 
#' This function creates a subsequence from a sequence with a 
#' random start position. The length of the subsequence depends
#' on the n-plet size size.
#'
#' @param seq Primary sequence of length N as a string.
#' @param n n-plet size (1 <= n <= nMax).
#' @param nMax maximal n-plet size.
#' @param multipleN If set to true, the start position is a multiple
#' of n, otherwise not.
#'
#' @return
#'
#' @examples
rndsubseq = function(seq, n, nMax, multipleN = TRUE) {
  if (n == nMax) {
    seq
  } else {
    N = nchar(seq)    # |X|; size of sequence
    xns = N / nMax    # |X_n*|; min. number of n-plets
    l = n * xns       # length of sub sequence
    i1 = if (multipleN) {
      h = (N - l) / n # Number of n-plets of size n minus the length.
      ((sample(h, size = 1) - 1) * n) + 1
    } else {
      sample(N - l, size = 1) # btw. 1 and N-l. (ell not 1)
    }
    i2 = i1 + l
    substr(seq, i1, i2)
  }
}


#' Sample-size normalized skews (To do).
#' 
#' Like 11.01.2019 in paper.
#'
#' @param seq DNA sequence
#' @param nSizes Vector of n-plet sizes (e.g. 1,2,3)
#' @param skewF Function to calculate the skews (default: Java version)
#'
#' @return A pair of AT and CG skews which are sample sized normalized.
#' Each entry may contain skews for different n-plet sizes (which
#' implies different positions).
#'
#' @export
#' 
#' @examples
nplet.sizenormskews = function(seq, nSizes = c(1, 2, 3), cf = 3,
                               skewF = nplet.skewj) {
  
  nMax = max(nSizes) # maximal n-plet size.
  normSkewAT = list() # List with normalized skews for different sizes.
  normSkewCG = list() # These lists will be filled below. 
  
  i = 1
  for (n in nSizes) { # for all n-plet sizes...

    subSeq = rndsubseq(seq, n, nMax, multipleN = TRUE) # Get a subsequence Y(n, nMax)
    skews = skewF(subSeq, n)

    # First iteration:
    rSkewsAT = list() # List for repeatitons (iterations)
    rSkewsCG = list()
    rSkewsAT[[1]] = skews$at
    rSkewsCG[[1]] = skews$cg
    
    iterations = as.integer(nMax / n * cf)
    repeats = iterations - 1 # Remaining repeats

    if (repeats > 0) { # to avoid wrong range 2:1 gives 2,1
      for (r in 2:iterations) {
        subSeq = rndsubseq(seq, n, nMax, multipleN = TRUE) # Y(n, nMax)
        skews = skewF(subSeq, n)
        rSkewsAT[[r]] = skews$at
        rSkewsCG[[r]] = skews$cg
      }
    }
    # Make matrix for calculation of mean values.
    # rows: positions (1 <= k <= n)
    # columns: iterations (1 <= r <= I)
    mskewsAT = matrix(unlist(rSkewsAT), nrow = n)
    mskewsCG = matrix(unlist(rSkewsCG), nrow = n)
    avgskewsAT = c() # List for average skews (for repeats); one entry for each k
    avgskewsCG = c()
    for (k in 1:n) { # all positions
      avgskewsAT = c(avgskewsAT, mean(mskewsAT[k,])) # S_A~T(X,n,k)
      avgskewsCG = c(avgskewsCG, mean(mskewsCG[k,])) # S_C~G(X,n,k)
    }
    normSkewAT[[i]] = avgskewsAT # 1 to n entries; for each position k one.
    normSkewCG[[i]] = avgskewsCG
    i = i + 1
  }
  names(normSkewAT) = nSizes # Assign n-plet sizes.
  names(normSkewCG) = nSizes  
  list(at = normSkewAT, cg = normSkewCG)  
}

#' Normalized skews based on partioning. (deprecated)
#'
#' @param seq DNA sequence
#' @param nSizes Vector of n-plet sizes (e.g. 1,2,3)
#'
#' @return
#' @export
#'
#' @examples
nplet.normpartskews_ = function(seq, nSizes = c(2, 4, 8), 
                               skewF = nplet.skewj) {

  nMax = max(nSizes)
  N = nchar(seq) # Size of the sequence
  nplets = N %/% nMax # Number of n-plets

  skewsAT = list() # List with nplets for differenz sizes.
  skewsCG = list()  
  
  print(paste("|X|: ", N, 
              "; n*: ", nMax,
              "; |X_n*|: ", nplets, sep = ""))
  i = 1
  for (n in nSizes) { # for all n-plet sizes...
    
    rskewsAT = list() # List for all partitions
    rskewsCG = list()

    l = (nplets * n) # subsequence length.
    noPart = N %/% l # Number of partitions.

    stepwide = max(1, noPart %/% 512) # Max 64 partitions for now
    #partitions = seq(1, noPart, by=stepwide)
    print("")
    print(paste("n: ", n, "; # part.: ", noPart, 
                "; length: ", l))
    
    j = 1
    for (p in 1:noPart) {
      i1 = (p - 1) * l + 1
      i2 = i1 + l - 1
      subSeq = substr(seq, i1, i2)

      nSkewsATCG = skewF(subSeq, n)
      rskewsAT[[j]] = nSkewsATCG$at
      rskewsCG[[j]] = nSkewsATCG$cg
      j = j + 1
    }
    if (N %% l != 0) {
      i1 = noPart * l + 1
      i2 = N
      size = i2 - i1
      percsize = size / N
      print(paste("remaining seq. length: ", size, " bases; ", 
                  sprintf("%1.1f", percsize * 100), " %", 
                  sep = ""))
      if (abs(noPart * l + size - N) > 1) {
        print(noPart * l + size)
        print(N)
        stop("Internal error!")
      }
    }
    
    # Make matrix for calculation of mean values.
    # rows: positions (1 <= k <= n)
    # columns: iterations (1 <= r <= I)
    mskewsAT = matrix(unlist(rskewsAT), nrow = n)
    mskewsCG = matrix(unlist(rskewsCG), nrow = n)
    avgskewsAT = c() # List for average skews
    avgskewsCG = c()
    for (k in 1:n) { # all positions
      avgskewsAT = c(avgskewsAT, mean(mskewsAT[k,]))
      avgskewsCG = c(avgskewsCG, mean(mskewsCG[k,]))
    }
    skewsAT[[i]] = avgskewsAT # C_{A~T}(n, k)
    skewsCG[[i]] = avgskewsCG # C_{C~G}(n, k)
    i = i + 1
  }
  names(skewsAT) = nSizes # Assign n-plet sizes.
  names(skewsCG) = nSizes  
  list(at = skewsAT, cg = skewsCG)  
}


#' Sample-size normalized skews. (deprecated)
#' 
#' This is the old version. Probably deprecated.
#'
#' @param seq DNA sequence
#' @param nSizes Vector of n-plet sizes (e.g. 1,2,3)
#' @param cf coverage factor. Default 1.
#' @param skewF Function to calculate the skews (default: Java version)
#'
#' @return
#'
#' @examples
nplet.normrandskews_ = function(seq, nSizes = c(1, 2, 3), cf = 1,
                                multipleN = FALSE,
                                skewF = nplet.skewj) {
  
  # Center skews such that they have a mean value of 0.
  centerskew = function(skews) {
    skews - mean(skews) 
  }
  
  nMax = max(nSizes) # maximal n-plet size.
  skewsAT = list()   # List with nplets for differenz sizes.
  skewsCG = list()  
  
  i = 1
  for (n in nSizes) { # for all n-plet sizes...
    
    nskewsAT = c()
    nskewsCG = c()
    
    subSeq = rndsubseq(seq, n, nMax, multipleN = multipleN) # Get a subsequence
    #subSeq = rndsubseq(seq, n, nMax) # Get a subsequence
    nSkewsATCG = skewF(subSeq, n)
    
    # First iteration:
    
    # Centered around all positions k:
    nskewsAT = c(nskewsAT, centerskew(nSkewsATCG$at))
    nskewsCG = c(nskewsCG, centerskew(nSkewsATCG$cg))
    
    #iterations = 1
    iterations = (nMax %/% n) * cf
    repeats = iterations - 1 # Remaining repeats
    
    if (repeats > 0) { # to avoid wrong range 2:1 gives 2,1
      for (r in 2:iterations) {
        subSeq = rndsubseq(seq, n, nMax) # 
        nSkewsATCG = skewF(subSeq, n)
        nskewsAT = c(nskewsAT, centerskew(nSkewsATCG$at))
        nskewsCG = c(nskewsCG, centerskew(nSkewsATCG$cg))
      }
    }
    # Add the values for the current n-plet size:
    skewsAT[[i]] = nskewsAT
    skewsCG[[i]] = nskewsCG
    i = i + 1
  }
  names(skewsAT) = nSizes # Assign n-plet sizes.
  names(skewsCG) = nSizes  
  list(at = skewsAT, cg = skewsCG)  
}
