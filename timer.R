timer <- function(input){
  time.start <- Sys.time()
  # input should be of the form
  # "name <- function(value)"
  # so that it actually is saved
  input
  time.end <- Sys.time()
  timediff <- time.end-time.start
  return(timediff)
}