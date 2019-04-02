#' #' @export
#' #' @importClassesFrom SingleCellExperiment 
#' setClass("SCEObjects", 
#'          slots = c(InputData = "list"),
#'          contains = "SingleCellExperiment")
#' setValidity("SCEObjects", function(object){
#'   if (!is(object,"list")) {
#'     stop("Error: input data has to be provided as a list of SingleCellExperiment objects \n")
#'   } 
#' })
#' 
#' #any class has to have a constructor and it is defined in the same page as the class
#' #extends= c("SingleCellExperiment")
#' #slots is where the data usually is
#' 
#' #S4 methods
#' #a method is a funciton that allows you to run different sets of code, based on different values of the argument
