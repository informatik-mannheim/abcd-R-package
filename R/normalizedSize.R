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

#' Calculate the n-plet skew for specific sequence length.
#' The subsequence starts at 5% of the entire sequence. This is
#' to avoid having many unclassified bases (N) which are usually
#' found at the beginnging of a sequence.
#'
#' @param seq Sequence as a string.
#' @param sizes Vector of sizes. The maximal size must be less than 95% 
#' of the size of the sequence.
#' @param n The n-plet size. Default is 10^3, 10^4, 10^5 and 10^6.
#' @param skewF Function for the calcuation of the skew: either nplet.skew
#' for the R version for nplet.skewJ for the Java version. 
#' Default is nplet.skewj.
#'
#' @return List with two entries: First entry contains the A-T skews for all
#' lenghts, second entry the C-G skews accordingly. 
#' @export
#'
#' @examples
nplet.nskewpersize = function(seq, sizes = c(1E3, 1E4, 1E5, 1E6), 
                              n = 3,
                              skewF = nplet.skewj) {
  
  N = nchar(seq) # Length of sequence.
  if (N < max(sizes)) {
    stop("Sequence is less than maximum size.")
  }
  
  nskewsAT = list()
  nskewsCG = list()
  i = 1
  for (size in sizes) {
    i1 = as.integer(0.05 * n)
    i2 = i1 + size - 1
    subseq = substr(seq, i1, i2)
    nskewsAT[[i]] = skewF(subseq, n)$at
    nskewsCG[[i]] = skewF(subseq, n)$cg
    i = i + 1
  }
  names(nskewsAT) = sizes
  names(nskewsCG) = sizes
  list(at=nskewsAT, cg=nskewsCG)
}

#' Calculates the standard deviation for n-plet skews.
#'
#' @param skews Data structure of type: <base-tuple>$<sizes>$<skews>, e.g.
#' `at`$1000$c(.01, 0.2, 0.1)
#'
#' @return
#' @export
#'
#' @examples
nplet.stddevpersize = function(skews) {
  list(at = lapply(skews$at, sd), cg = lapply(skews$cg, sd))
}

#' Plot the standard deviation of skews for a specific length as
#' a bar diagram.
#'
#' @param nskewssd Nested data strucuture as returned by nplet.stddevpersize.
#' @param seqid A brief description of the sequence.
#' @param yMax Maximum range of y-axis. Default is 0.05.
#'
#' @return
#' @export
#'
#' @examples
nplet.stddevplotpersize = function(nskewssd, 
                                   seqId = "Unkn. sequence",
                                   yMax = 0.05,
                                   xtickssize = 2) {
  nat = names(nskewssd$at[1:length(nskewssd$at)]) # Retrieve seq. sizes
  ncg = names(nskewssd$cg[1:length(nskewssd$cg)])
  
  dtick = length(nat) / xtickssize # Number of ticks - 1
  xticks = seq(dtick, length(nat), by = dtick)
  xlabels = c(nat[1], nat[xticks])
  xlabels = lapply(xlabels, function (l) sprintf("%1.1e", as.numeric(l)))
  space = 0.5
  pos = c(1, (xticks * (1 + space)) - space) # Position of x labels.

  par(mfrow=c(1, 2)) # Skews side by side
  barplot(unlist(nskewssd$at), names.arg = nat, space = .5,
          main = seqId, xlab = "seq. size", 
          ylab = "std. dev. of A~T skews",
          ylim = c(0, yMax), xaxt="n")  
  axis(1, labels = xlabels, at = pos, las = 1)
  grid (0,NULL, lty = 6, col = "cornsilk2") 
  
  barplot(unlist(nskewssd$cg), names.arg = ncg, space = .5,
          main = seqId, xlab = "seq. size", 
          ylab = "std. dev. of C~G skews",
          ylim = c(0, yMax), xaxt="n")  
  axis(1, labels = xlabels, at = pos, las = 1)  
  grid (0,NULL, lty = 6, col = "cornsilk2") 
  par(mfrow=c(1, 1)) # Switch it off again?
}
