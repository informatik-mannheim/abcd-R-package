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
