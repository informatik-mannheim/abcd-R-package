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

#' Box and whisker plots for A~T and C~G n-plet-skews.
#'
#' @param nSkews Data structure for n-skews as 
#' returned by \code{\link{nplet.skews}}.
#' @param seqId A short description/label for the sequence.
#' @param ylabPrefix A prefix for the actual y label. E.g. "normalized ".
#'
#' @export
#'
#' @examples
nplet.skewplot = function(nSkews, yMax = 1, seqId = "Unknown sequence",
                          ylabPrefix = "") {
  
  yRange = c(-yMax, yMax)
  
  nat = names(nSkews$at[1:length(nSkews$at)]) # Retrieve n-plet sizes
  ncg = names(nSkews$cg[1:length(nSkews$at)])
  
  meanAT = lapply(nSkews$at, mean) # Mean values
  meanCG = lapply(nSkews$cg, mean)
  
  ylabTextAT = paste(ylabPrefix, "A~T skew")
  ylabTextCG = paste(ylabPrefix, "C~G skew")
  
  par(mfrow=c(1, 2)) # Skews side by side
  # Better ignore parameter during development: ylim = yRange
  boxplot(nSkews$at, main = seqId,
          names = nat, ylim = yRange,
          xlab = "n-plet size", ylab = ylabTextAT, xaxt="n")
  points(1:length(nat), meanAT, pch = 23, # circle
         cex = .75, bg = "gray") 
  grid (0,NULL, lty = 6, col = "cornsilk2")
  axis(1, labels = nat, at=1:length(nat), las=2)
  
  boxplot(nSkews$cg, main = seqId,
          names = ncg, ylim = yRange,
          xlab = "n-plet size", ylab = ylabTextCG, xaxt="n")
  points(1:length(ncg), meanCG, pch = 23, cex = .75, bg = "gray")  
  grid (0,NULL, lty = 6, col = "cornsilk2") 
  axis(1, labels = ncg, at=1:length(ncg), las=2)
  par(mfrow=c(1, 1)) # Switch it off again?
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
                            yMax = 0.01, ylabPrefix = "") {
  
  stdev = function(skews) {
    vat = lapply(skews$at, sd)
    vcg = lapply(skews$cg, sd)
    
    list(at=vat, cg=vcg)
  }
  
  v = stdev(skews)

  ylabTextAT = paste(ylabPrefix, "std. dev. of A~T skew")
  ylabTextCG = paste(ylabPrefix, "std. dev. of C~G skew")
  
  par(mfrow=c(1, 2)) # Skews side by side
  barplot(unlist(v$at), names.arg = nSizes, space = .5,
          main = seqId, xlab = "n-plet size", 
          ylab = ylabTextAT,
          ylim = c(0, yMax))
  
  barplot(unlist(v$cg), names.arg = nSizes, space = .5,
          main = seqId, xlab = "n-plet size", 
          ylab = ylabTextCG,
          ylim = c(0, yMax))  
  par(mfrow=c(1, 1)) # Switch it off again?
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
