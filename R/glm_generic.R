#' @name print.bats
#' @title Print function for BATS outputs
#' @description Print method function for objects of class 'bats'.
#' @param x An object of class 'bats'.
#' @param ... Additional arguments affecting the print produced.
#' @export
print.bats = function(x, ...){
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
        cat("  *",x$par$RAR,"\n")        
    }
    cat("\n Model:\n")
    cat("  *",format(x$call$model)," (with",x$call$link, "link)\n")
    cat("\n Call:\n ")
    print(x$call)
}


 



