#probably needs a good bit of error checking written in
#
landscape.new.expression <-
function(rland,expmat,hsq)
{
  rland$expression <- list(expmat=expmat,hsq=hsq)
  rland
}

