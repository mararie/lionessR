#' LIONESS
#'
#' This function uses the LIONESS equation to estimate single-sample networks.
#' @param x Numeric matrix with samples in columns.
#' @param f Network reconstruction function. Defaults to Pearson correlation.
#' @keywords lioness
#' @export
#' @examples
#' exp <- matrix(sample(1000,1000)/1000, 100, 10)
#' colnames(exp) <- paste("sample", c(1:ncol(exp)), sep="_")
#' row.names(exp) <- paste("gene", c(1:nrow(exp)), sep="_")
#' lionessResults <- lioness(exp, netFun)

lioness <- function(x, f=netFun, ...){
 
  if(!is.function(f)){ stop("please use a function") }
  if(!is.matrix(x)) { print("please use a numeric matrix as input") }
    
  nrsamples <- ncol(x)
  samples <- colnames(x)

  # this applies netFun and extracts the aggregate network
  net <- f(x)
  agg <- c(net)

  # prepare the lioness output
  lionessOutput <- matrix(NA, nrow(net)*ncol(net), nrsamples+2)
  colnames(lionessOutput) <- c("reg", "tar", samples)
  lionessOutput[,1] <- rep(row.names(net), ncol(net))
  lionessOutput[,2] <- rep(colnames(net), each=nrow(net))
  lionessOutput <- as.data.frame(lionessOutput, stringsAsFactors=F)
  lionessOutput[,3:ncol(lionessOutput)] <- sapply(lionessOutput[,3:ncol(lionessOutput)], as.numeric)

  # run function f and the LIONESS equation
  for(i in 1:nrsamples){
    ss <- c(f(x[,-i])) # apply netFun on all samples minus one
    lionessOutput[,i+2] <- nrsamples*(agg-ss)+ss # apply LIONESS equation
  }
  return(lionessOutput)  
}
