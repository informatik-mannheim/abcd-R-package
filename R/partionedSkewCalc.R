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

#' Create sub sequences.
#' @author Markus Gumbel
#' @param seq Primary sequence of length N.
#' @param pCount Number of partitions (1 <= pCount).
#' @param p i-th partition (1 <= i <= pCount).
#'
#' @return
#' @export
#'
#' @examples
subseqpart = function(seq, pCount = 1, p = 1) {
  if (pCount == 1) {
    seq
  } else {
    N = nchar(seq)
    pLength = N %/% pCount # int. division
    i1 = (p - 1) * pLength + 1
    i2 = i1 + pLength - 1
    substr(seq, i1, i2)
  }
}

#' A~T and C~G skews for partitions of a sequence.
#' 
#' @param seq DNA sequence
#' @param nSize n-plet size. Default is 1.
#' @param pCount Number of partitions. Default is 32.
#'
#' @return
#' @export
#'
#' @examples
nplet.partskew = function(seq, nSize = 1, pCount = 32,
                          skewF = nplet.skewj) {
  skewsAT = list()
  skewsCG = list()
  
  partitions = 1:pCount
  
  i = 1
  for (p in partitions) { # for all n-plet sizes...
    subSeq = subseqpart(seq, pCount, p) # 
    nSkewsATCG = skewF(subSeq, nSize)
    skewsAT[[i]] = nSkewsATCG$at
    skewsCG[[i]] = nSkewsATCG$cg
    i = i + 1
  }
  names(skewsAT) = partitions
  names(skewsCG) = partitions
  list(at = skewsAT, cg = skewsCG)
}

#' Box and whisker plots of A~T and C~G skews for a list of
#' partitions (subsequences).
#'
#' @param seq DNA sequence
#' @param nSizes Vector of n-plet sizes (e.g. 1,2,3)
#' @param pCount Number of partitions
#' @param yMax +/- max. y range.
#'
#' @return
#' @export
#'
#' @examples
nplet.partskewplot = function(seq, nSize = 1, 
                              pCount = 32, yMax = 1,
                              seqId = "Unknown sequence",
                              skewF = nplet.skewj) {
  skews = nplet.partskew(seq, nSize, pCount, skewF)
  yRange = c(-yMax, yMax) # Range on y-axis.

  meanAT = lapply(skews$at, mean) # Mean values
  meanCG = lapply(skews$cg, mean)
  
  dtick = pCount / 8 # Number of ticks
  xticks = c(1, seq(dtick, pCount, by = dtick))
  
  par(mfrow=c(1, 2)) # Skews side by side
  # Better ignore parameter during development: ylim = yRange
  title = paste(seqId, "\npart. skew; n = ", nSize)
  boxplot(skews$at, main=title, ylim = yRange, 
          xlab = "partition", ylab = "A~T skew",
          xaxt = "n")
  points(1:pCount, meanAT, pch = 23, cex = .15, bg = "gray")  
  grid (0,NULL, lty = 6, col = "cornsilk2") 
  axis(1, at = xticks, labels = xticks)
  
  boxplot(skews$cg, main=title, ylim = yRange, 
          xlab = "partition", ylab = "C~G skew",
          xaxt = "n")
  points(1:pCount, meanCG, pch = 23, cex = .15, bg = "gray")    
  grid (0,NULL, lty = 6, col = "cornsilk2") 
  axis(1, at = xticks, labels = xticks)
  
  par(mfrow=c(1, 1)) # Switch it off again?
}

#' Box and whisker plots for A~T and C~G skews for a list of
#' partitions read from a fasta file.
#'
#' @param filename Path to fasta filename.
#' @param nSizes Vector of n-plet sizes (e.g. 1,2,3)
#' @param pCount Number of partitions
#' @param yMax +/- max. y range.
#'
#' @return
#' @export
#'
#' @examples
nplet.filepartplot = function(filename, nSize = 1, 
                              pCount = 32, yMax = 0.1,
                              skewF = nplet.skewj) {
  seq = sequence.fastaread(filename)
  nplets.partskewplot(seq, nSize, pCount, yMax)
}
