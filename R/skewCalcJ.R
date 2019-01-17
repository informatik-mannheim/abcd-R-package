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

# skew related functions which are implemented in Java.

#' A~T and C~G skews for an n-plet. 
#' 
#' This is the Java-based version.
#'
#' @param seq DNA sequence as a string. Considered symbols are 
#' in upper or lower case: ' A, T, C, G. 
#' U is also considered and replaced by a T.
#' @param n n-plet size. Default is 1.
#'
#' @return list of skews: first entry: vector of A~T skews,
#' second entry: vector of C~G skews.
#' @export
#'
#' @examples
nplet.skewj = function(seq, n = 1) {
  # Create Java object:
  nc = new(J("bio/gcat/abcdtool/analysis/NPletCalc"), 
           seq, as.integer(n))
  skewsAT = nc$getSkewsR("A", "T")
  skewsCG = nc$getSkewsR("C", "G")
  names(skewsAT) = 1:n
  names(skewsCG) = 1:n
  res = list("at" = skewsAT, "cg" = skewsCG)
  gc()
  res
}

#' Skew for two given bases for an n-plet. 
#' 
#' This is the Java-based version.
#'
#' @param seq DNA sequence as a string. Considered symbols are 
#' in upper or lower case: ' A, T, C, G. 
#' U is also considered and replaced by a T.
#' @param base1 The first base as a string.
#' @param base2 The second base as a string.
#' @param n n-plet size. Default is 1.
#'
#' @return list of skews: first entry: vector of A~T skews,
#' second entry: vector of C~G skews.
#' @export
#'
#' @examples
nplet.skewbybase = function(seq, base1, base2, n = 1) {
  # Create Java object:
  skewb = new(J("bio/gcat/abcdtool/analysis/NPletCalc"),
             seq, as.integer(n))
  res = skewb$getSkewsR(base1, base2)
  gc()
  res
}