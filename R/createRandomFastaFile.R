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

# Note: package seqinr converts bases to lower case.
NUCBASES = c("A", "T", "C", "G")

writeFastaFile = function(filename, header, size, probs) {

  idxRange = 1:4 # Index for the four bases.
  
  rndCondBase = function(base) {
    nextProbs = probs[, base] # Probs for next base.
    rndIdx = sample(idxRange, 1, prob=nextProbs) # Pick a random position.
    NUCBASES[rndIdx]
  }
  
  # Not used at the moment:
  rndBase = function() {
    rndIdx = sample(1:4, 1) # Random number btw. 1 and 4.
    bases[rndIdx]
  }
  
  # Inner function, great!
  rndBases = function(size) {
    b = rep(NA, size)
    base = "A" # We always start with an A as the prev. base.
    for (i in 1:size) {
      base = rndCondBase(base)
      b[i] = base
    }
    paste(b, collapse = "") # Make a single string.
  }
  colSize = 80 # Default number of bases per row.
  fullLines = size %/% colSize # Number of lines.
  restChars = size %% colSize  # Number of remaining bases.
  
  lines = rep(NA, fullLines + 1)
  lines[1] = paste(header, ", size: ", size)
  
  for (i in 2:(fullLines + 1)) {
    lines[i] = rndBases(colSize)
  }
  lines = c(lines, rndBases(restChars))
  
  # Now persist lines to file:
  fileConn <- file(filename)
  writeLines(lines, fileConn)
  close(fileConn)
}

defaultWriteFastaFile <- function(size) {
  t1 = Sys.time() # Start time.
  filename = "H_habilis_chr1_rnd.fasta"
  writeFastaFile(filename, "> random", size, chr1CondProbs)
  print(Sys.time() - t1)
}

chr1CondProbs = matrix(
  c(0.3265463, 0.2552814, 0.1729273, 0.24524500,
    0.2164471, 0.3278836, 0.2058872, 0.24978216, 
    0.3489394, 0.3422078, 0.2594270, 0.04942584, 
    0.2877819, 0.2417155, 0.2108700, 0.25963262), 
  nrow = 4, ncol = 4,
  dimnames = list(NUCBASES, NUCBASES)
)

chr1Probs = matrix(
  c(0.290964, 0.291815, 0.208496, 0.208724,
    0.290964, 0.291815, 0.208496, 0.208724, 
    0.290964, 0.291815, 0.208496, 0.208724,
    0.290964, 0.291815, 0.208496, 0.208724), 
  nrow = 4, ncol = 4,
  dimnames = list(NUCBASES, NUCBASES)
)