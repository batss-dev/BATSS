#' @name print.batss
#' @title Print function for BATSS outputs
#' @description Print method function for objects of class 'batss' (i.e., output of the function [batss.glm]).
#' @param x An object of class 'batss'.
#' @param ... Additional arguments affecting the print produced.
#' @seealso [batss.glm()], the function generating S3 objects of class 'batss'. 
#' @returns Prints information for objects of class 'batss'.
#' @export
print.batss = function(x, ...){
    cat("\n Bayesian adaptive trial Monte Carlo simulation.\n")
    if(!is.null(x$par$RAR)){
        cat(" with response adapative randomisation (",
            x$call$RAR,
            ")\n",sep="")
    }
    cat(" (",length(x$par$seed)," Monte Carlo samples)\n",sep="")
    cat("\n Variables:\n")
    for(i in 2:length(x$call$var)){
        cat("  *",names(x$call$var)[i],":",as.character(x$call$var)[i],"\n")
    }
    if(!is.null(x$par$RAR)){
        cat("\nGroup randomisation:\n")
        cat("  *",x$call$RAR,"\n")
    }
    cat("\n Model:\n")
    cat("  *",format(x$call$model),if (x$type=="surv") "\n" else paste("(with",ifelse(is.null(x$call$link),"identity",x$call$link), "link)\n"))
    cat("\n Call:\n ")
    print(x$call)
}






