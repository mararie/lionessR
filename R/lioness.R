#' LIONESS
#'
#' This function uses the LIONESS equation to estimate single-sample networks.
#' @param x Numeric matrix with samples in columns.
#' @param f Network reconstruction function. Defaults to Pearson correlation.
#' @param assay Specified the assay used in lionessR for SummerizedExperiment object. Default: counts. 
#' @keywords lioness
#' @export
#' @examples
#' exp <- matrix(sample(1000,1000)/1000, 100, 10)
#' colnames(exp) <- paste("sample", c(1:ncol(exp)), sep="_")
#' row.names(exp) <- paste("gene", c(1:nrow(exp)), sep="_")
#' lionessResults <- lioness(exp, netFun)

lioness <- function(x, f=netFun, assay='counts'){

    is.se <- inherits(x, "SummarizedExperiment")
    is.matrix <- is.matrix(x)
    
    if(!is.function(f)){ stop("please use a function") }
    if(is.matrix(x)) { print("take numeric matrix as input, ignore parameter for assay")}
    if(is.se) { x <- SummarizedExperiment::assays(x)[[assay]] }
    if(!is.matrix(x)) { print("please use a numeric matrix as input")}
        
    if(is.null(colnames(x))){ colnames(x) = seq_len(ncol(x)) }
    
    nrsamples <- ncol(x)
    samples <- colnames(x)

    # this applies netFun and extracts the aggregate network
    net <- netFun(x)
    agg <- c(net)

    # prepare the lioness output
    lionessOutput <- matrix(NA, nrow(net)*ncol(net), nrsamples+2)
    colnames(lionessOutput) <- c("reg", "tar", samples)
    lionessOutput[,1] <- rep(row.names(net), ncol(net))
    lionessOutput[,2] <- rep(colnames(net), each=nrow(net))
    lionessOutput <- as.data.frame(lionessOutput, stringsAsFactors=FALSE)
    lionessOutput[,3:ncol(lionessOutput)] <- vapply(lionessOutput[, 3:ncol(lionessOutput)], 
                                                    as.numeric, 
                                                    vector('numeric', nrow(lionessOutput)))

    # run function f and the LIONESS equation
    for(i in seq_len(nrsamples)){
        ss <- c(f(x[,-i])) # apply netFun on all samples minus one
        lionessOutput[,i+2] <- nrsamples*(agg-ss)+ss # apply LIONESS equation
    }
    return(lionessOutput)  
}
