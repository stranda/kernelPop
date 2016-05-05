#
# a convenience function that provides the highest column number for demographic information in
# the individuals matrix
#
landscape.democol <- function()
  {
    as.integer(.Call("num_demo_cols",PACKAGE="kernelPop2"))
  }
