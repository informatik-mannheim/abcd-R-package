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

# This section shows the slow R calculation for skews.

#' Calculates the base distribution of a sequence.
#'
#' @param seq The sequence as a vector chars or as regular
#' string. Considered symbols are in upper or lower case:
#' A, T, C, G. U is also considered and replaced by a T.
#'
#' @return Vector of distribution with named entries a, t, c, g.
#' @export
#'
#' @examples
sequence.basedist = function(seq) {
  if (length(seq) == 1) {
    seq = strsplit(seq, "")[[1]]
  }
  baseA = c(seq[seq=="a"], seq[seq=="A"])  # New sequence with A bases only.
  baseT = c(seq[seq=="t"], seq[seq=="T"], seq[seq=="u"], seq[seq=="U"])
  baseC = c(seq[seq=="c"], seq[seq=="C"])
  baseG = c(seq[seq=="g"], seq[seq=="G"])
  nA = length(baseA) # Number of bases.
  nT = length(baseT)
  nC = length(baseC)
  nG = length(baseG)
  
  c(a=nA, t=nT, c=nC, g=nG)
}

#' Calculates the A~T and C~G skews. 
#' 
#' This is the slow R-based version.
#'
#' @param seq The sequence as a vector chars or as regular
#' string. Considered symbols are in upper or lower case:
#' A, T, C, G. U is also considered and replaced by a T.
#'
#' @return Vector with two elements: (skew A-T, skew C-G)
#' @export
#'
#' @examples
skew = function(seq) {
  baseCount = sequence.basedist(seq)

  nA = baseCount[["a"]]; nT = baseCount[["t"]]
  nC = baseCount[["c"]]; nG = baseCount[["g"]]
  
  # The skew as defined.
  skewAT = (nA - nT) / (nA + nT)
  skewCG = (nC - nG) / (nC + nG)

  c(at = skewAT, cg = skewCG)
}

# Here follow the R only implementations.

#' Calculates the skews for an n-plet of size n.
#'
#' This is the slow R-based version.
#' 
#' @param seq The DNA sequence to analyze. This can be a vector
#' of characters or a regular string.
#' @param n n-plet size. Default is 1.
#'
#' @return List with two entries: First entry is the A-T skew,
#' second entry is the C-G skew. Eeach entry consists of a list
#' of skews for each position k in the n-plet (1 <= k <= n).
#' @export
#'
#' @examples
nplet.skew = function(seq, n = 1) {
  if (length(seq) == 1) {
    seq = strsplit(seq, split = "")[[1]]
  } 

  skewsAT = c() # Empty result vector.
  skewsCG = c()
  
  N = length(seq) # Bases in the sequence.
  # Create list of size N with 
  # pattern 1, 2, ..., n, 1, 2, ..., n, ...
  indizes = rep(1:n, N / n) 
  # Create n new sequences based on theses indizes:
  posSeq = split(seq, indizes) 
  # Iterate over the positions k (1 <= k <= n):
  for (k in 1:n) {
    #print(paste("k: ", k))
    skewATCG = skew(posSeq[[k]])
    skewsAT = c(skewsAT, skewATCG["at"])
    skewsCG = c(skewsCG, skewATCG["cg"])
  }
  names(skewsAT) = 1:n
  names(skewsCG) = 1:n
  list(at=skewsAT, cg=skewsCG)
}

#' Calculates the skews for a list of n-plet sizes.
#' 
#' This function iterates over the n-plet sizes and calculates  
#' the skews for each n-plet size.
#' 
#' @param seq The DNA sequence to analyze. This can be a vector
#' of characters or a regular string.
#' @param nSizes Vector of n-plet sizes. Default is (1,2,3).
#' @param skewF Function for the calcuation of the skew: either 
#' \code{\link{nplet.skew}} for the R version 
#' or \code{\link{nplet.skewJ}} for the Java version.
#' Default is \code{nplet.skewJ}.
#'
#' @return
#' @export
#'
#' @examples
nplet.skews = function(seq, nSizes = c(1,2,3), skewF = nplet.skewj) {
  skewsAT = list() # List with nplets for differenz sizes.
  skewsCG = list()
  
  i = 1 # for index
  for (n in nSizes) { # for all n-plet sizes...
    nSkewsATCG = skewF(seq, n)
    skewsAT[[i]] = nSkewsATCG$at
    skewsCG[[i]] = nSkewsATCG$cg
    i = i + 1
  }
  names(skewsAT) = nSizes # Assign n-plet sizes.
  names(skewsCG) = nSizes  
  list(at = skewsAT, cg = skewsCG)
}
