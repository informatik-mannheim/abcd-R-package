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

#' Sample-size normalized skews.
#' This is the old version. Probably deprecated.
#'
#' @param seq DNA sequence
#' @param nSizes Vector of n-plet sizes (e.g. 1,2,3)
#' @param cf coverage factor. Default 1.
#' @param skewF Function to calculate the skews (default: Java version)
#'
#' @return
#' @export
#'
#' @examples
nplet.normrandskews_ = function(seq, nSizes = c(1, 2, 3), cf = 1,
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
    
    subSeq = rndsubseq(seq, n, nMax) # Get a subsequence
    nSkewsATCG = skewF(subSeq, n)
    
    # First iteration:
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
    skewsAT[[i]] = nskewsAT
    skewsCG[[i]] = nskewsCG
    i = i + 1
  }
  names(skewsAT) = nSizes # Assign n-plet sizes.
  names(skewsCG) = nSizes  
  list(at = skewsAT, cg = skewsCG)  
}


#' Sample-size normalized skews (To do).
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
      avgskewsAT = c(avgskewsAT, mean(mskewsAT[k,]))
      avgskewsCG = c(avgskewsCG, mean(mskewsCG[k,]))
    }
    normSkewAT[[i]] = avgskewsAT # 1 to n entries; for each position k one.
    normSkewCG[[i]] = avgskewsCG
    i = i + 1
  }
  names(normSkewAT) = nSizes # Assign n-plet sizes.
  names(normSkewCG) = nSizes  
  list(at = normSkewAT, cg = normSkewCG)  
}

#' Create normalized skews based on partioning.
#'
#' @param seq DNA sequence
#' @param nSizes Vector of n-plet sizes (e.g. 1,2,3)
#'
#' @return
#' @export
#'
#' @examples
nplet.normpartskews = function(seq, nSizes = c(2, 4, 8), 
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

#' Plot the standard deviation of n-plet skews.
#'
#' @param skews List of n-plet skews. The first entry is the 
#' AT skew (label: at), the second entry is the CG skew (label: cg)
#' @param seqId A sequence id shown in the plot. 
#' Default is "Unkn. sequence".
#' @param yMax Maximum range of the y-axis. Default is 0.01.
#'
#' @return
#' @export
#'
#' @examples
nplet.sdskewplot = function(skews, nSizes, seqId = "Unkn. sequence",
                            yMax = 0.01) {
  
  stdev = function(skews) {
    vat = lapply(skews$at, sd)
    vcg = lapply(skews$cg, sd)
    
    list(at=vat, cg=vcg)
  }
  
  v = stdev(skews)
  # plot(nSizes, v$at, type = "p",
  #      main = seqId, xlab = "n", ylab = "std. dev. of A~T skews",
  #      ylim = c(0, yMax))
  # plot(nSizes, v$cg, type = "p",
  #      main = seqId, xlab = "n", ylab = "std. dev. of C~G skews",
  #      ylim = c(0, yMax))

  par(mfrow=c(1, 2)) # Skews side by side
  barplot(unlist(v$at), names.arg = nSizes, space = .5,
          main = seqId, xlab = "n-plet size", 
          ylab = "std. dev. of A~T skews",
          ylim = c(0, yMax))
  
  barplot(unlist(v$cg), names.arg = nSizes, space = .5,
       main = seqId, xlab = "n-plet size", 
       ylab = "std. dev. of C~G skews",
       ylim = c(0, yMax))  
  par(mfrow=c(1, 1)) # Switch it off again?
}

#' Create the box and whisker plots for a set of
#' n-plet sizes for a sequence.
#'
#' @param seq DNA sequence
#' @param nSizes Vector of n-plet sizes (e.g. 1,2,3)
#' @param yMax +/- max. y range.
#'
#' @return
#' @export
#'
#' @examples
nplet.normskewplot = function(nSkews, seqId = "Unkn. sequence",
                              yMax = 0.25) {
  yRange = c(-yMax, yMax)

  nat = names(nSkews$at[1:length(nSkews$at)]) # Retrieve n-plet sizes
  ncg = names(nSkews$cg[1:length(nSkews$at)])
  
  meanAT = lapply(nSkews$at, mean) # Mean values
  meanCG = lapply(nSkews$cg, mean)
  
  par(mfrow=c(1, 2)) # Skews side by side
  
  title = seqId
  boxplot(nSkews$at, main = title, names = nat, 
          xlab = "n-plet size", ylab = "norm. A~T skew",
          ylim = yRange)
  points(1:length(nat), meanAT, pch = 23, # circle
         cex = .75, bg = "gray")  
  grid (0,NULL, lty = 6, col = "cornsilk2") 
  
  boxplot(nSkews$cg, main = title, names = ncg, 
          xlab = "n-plet size", ylab = "norm. C~G skew",
          ylim = yRange)
  points(1:length(ncg), meanCG, pch = 23, cex = .75, bg = "gray")    
  grid (0,NULL, lty = 6, col = "cornsilk2") 
  par(mfrow=c(1, 1)) # Switch it off again?
}

#' Box and whisker plot for normalized A~T and C~G skews.
#'
#' @param filename Path to a fasta file.
#' @param nSizes Vector of n-plet sizes (e.g. 1,2,3)
#' @param yMax +/- max. y range.
#' 
#' Additional parameters like in \code{\link{nplet.normpartskews}} 
#' function possible.
#'
#' @return
#' @export
#'
#' @examples
nplet.filenormskewplot = function(filename, yMax = 0.25, 
                                  seqId = "Unkn. sequence", ...) {
  seq = sequence.fastaread(filename)
  nskews = nplet.normpartskews(seq, ...)
  nplet.normskewplot(nskews, yMax = yMax, seqId = seqId)
}
