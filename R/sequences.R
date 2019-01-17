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

#' Conditional probabilities for bases in a sequence.
#'
#' @param seq 
#'
#' @return
#' @export
#'
#' @examples
sequence.condbaseprob = function(seq) {
  cp = new(J("bio/gcat/abcdtool/sequences/generator/ConditionalProbabilities"))
  ja = cp$createConditionalProbabilityMatrix(seq)
  rm = .jevalArray(ja, simplify = TRUE)
  basenames = c("A", "T", "C", "G", "N")
  rownames(rm) = basenames
  colnames(rm) = basenames
  rm
}

#' Get a random sequence with conditional probabilities like in 
#' Human chromosome 1.
#'
#' @param size Number of bases in the sequence.
#'
#' @return Sequence as a string. Bases are upper case.
#' @export
#'
#' @examples
sequence.rndbycondprob = function(size, condprob) {
  # Create Java object:
  cp = .jarray(condprob, dispatch = TRUE)
  rg = new(J("bio/gcat/abcdtool/sequences/generator/RandomSeqStringGenerator"))
  rg$rndCondSeqString(rJava::.jlong(size),  # size is cast into a long with .jlong()
                      cp) 
}

#' Reads a single fasta file entry.
#' 
#' This function reads the first entry, i.e. the first sequence of
#' a fasta file. If there are more than one entry, they are ignored.
#'
#' @param filename Path to fast file
#'
#' @return Sequence as a string.
#' @export
#'
#' @examples
sequence.fastaread = function(filename) {
  seqinr::read.fasta(filename, as.string = TRUE, 
                     seqonly = TRUE, forceDNAtolower = FALSE)[[1]]
}

#' Creates a fasta file from a sequence.
#'
#' @param filename 
#' @param seq 
#' @param headertext 
#'
#' @return
#' @export
#'
#' @examples
sequence.fastawrite = function(filename, seq, headertext) {
  ffw = new(J("bio/gcat/abcdtool/sequences/writer/FastaFileWriter"))
  ffw$writeSequence(seq, filename, headertext)
}
