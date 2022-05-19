#' Print the name and the value of a variable.
#' 
#' @param x The object to display.
#' @return The input variable.
#' @examples
#' \dontrun{
#' data.dir = '/path/to/dir'
#' Show(data.dir)
#' }
#' 
Show = function(x) {
  name = as.character(substitute(x))
  repr = capture.output(dput(x))
  cat(sprintf("%s = %s\n", name, repr))
  invisible(x)
}

#' Assert a length-1 vector.
#' 
#' Execute an error if the length of the input is 0 and 
#' execute a warning if the length of the input is bigger 
#' than 1. Return the first element of the input.
#' 
#' @param x The object for checking.
#' @return The first element of the input.
#' @examples
#' \dontrun{
#' data.dir = '.'
#' Show(data.dir)
#' assert_length_1(data.dir)
#' data.dir = c('.', '/path/to/dir')
#' Show(data.dir)
#' assert_length_1(data.dir)
#' }
#' 
assert_length_1 = function(x) {
  if (length(x) < 1) {
    name = as.character(substitute(x))
    stop(sprintf('No %s was found!', name))
  } else if (length(x) > 1) {
    name = as.character(substitute(x))
    warning(sprintf('More than one %s!', name))
  }
  return(x[[1]])
}

