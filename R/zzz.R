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

# Access Java routines.
# requires rJava 0.9.9 or higher

# Some note:
# If a method cannot be found but it is defintely there, just restart R!
# Do not forget to wrap the primtive Java data types.
# int -> as.integer()
# char -> .jchar()
# etc.

#' Title
#'
#' @param libname 
#' @param pkgname 
#'
#' @return
#' @export
#'
#' @examples
.onLoad = function(libname, pkgname) {
  javaInit("")
  rJava::.jpackage(pkgname, jars='*', morePaths='', 
                   nativeLibrary=FALSE, lib.loc=libname)
}

javaInitDevelopment = function() {
  javaInit("c:\\Users\\Markus\\
           Local-Docs\\src\\jvm\\ABCD-Tool\\target\\classes\\")
}

javaInit = function(classPath) {
  # force.init = TRUE throws an error?
  res = rJava::.jinit(classpath=classPath, parameters="-Xmx4096m",
                      force.init = FALSE)
  if (res == 0) {
    packageStartupMessage("JVM is ready.")
  } else {
    packageStartupMessage((paste("JVM is not available. Error is ", res)))
  }
}