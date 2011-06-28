
select.mpt <- function(mpt.results, output = c("standard", "full"), round.digit = 6) {
	if(!is.list(mpt.results)) stop("mpt.results need to be a list.")
	if(length(mpt.results)< 2) stop("length(mpt.results) needs to be >= 2.") 
	n.models <- length(mpt.results)
	observed.data <- lapply(mpt.results, function(x) return(x[["data"]][["observed"]]))
	equal.data <- sapply(observed.data, function(x) identical(observed.data[[1]],x))
	if(!all(equal.data)) stop("observed data for the models differ.")
	n.data <- dim(observed.data[[1]])[1]
	if (!is.null(names(mpt.results))) m.names <- names(mpt.results)
	else m.names <- 1:n.models
	
	if (n.data == 1) {
		c.fia <- sapply(mpt.results, function(x) any(grepl("^FIA$", colnames(x[["information.criteria"]]))))
		if (any(c.fia)) if (all(!c.fia)) warning(paste("FIA not available for model(s):", paste(m.names[which(c.fia == FALSE)], collapse = ", ")))
		n.parameters <- sapply(mpt.results, function(x) x[["model.info"]][["n.parameters"]])
		if (any(c.fia)) {
			FIA <- vapply(mpt.results, function(x) tryCatch(x[["information.criteria"]][["FIA"]], error = function(e) NA), 0)
			delta.FIA <- FIA - min(FIA, na.rm = TRUE)
		}
		AIC <- sapply(mpt.results, function(x) x[["information.criteria"]][["AIC"]])
		BIC <- sapply(mpt.results, function(x) x[["information.criteria"]][["BIC"]])
		delta.AIC <- AIC - min(AIC)
		denom.wAIC <- sum(exp(-0.5*(delta.AIC)))
		wAIC <- sapply(delta.AIC, function(x) exp(-0.5*(x))/denom.wAIC)
		delta.BIC <- BIC - min(BIC)
		denom.wBIC <- sum(exp(-0.5*(delta.BIC)))
		wBIC <- sapply(delta.BIC, function(x) exp(-0.5*(x))/denom.wBIC)
		df.out <- data.frame(model = m.names, n.parameters)
		if (any(c.fia)) {
			df.out <- cbind(df.out, delta.FIA)
			if (output[1] != "standard") df.out <- cbind(df.out, FIA)
		}
		if (output[1] == "standard") {
			df.out <- cbind(df.out, delta.AIC, wAIC, delta.BIC, wBIC)   
		} else {
			df.out <- cbind(df.out, delta.AIC, wAIC, AIC, delta.BIC, wBIC, BIC)
		}
	} else {
		c.fia <- sapply(mpt.results, function(x) any(grepl("^FIA$", colnames(x[["information.criteria"]][["individual"]]))))
		if (any(c.fia)) if (all(!c.fia)) warning(paste("FIA not available for model(s):", paste(m.names[which(c.fia == FALSE)], collapse = ", ")))
		n.parameters <- sapply(mpt.results, function(x) x[["model.info"]][["aggregated"]][["n.parameters"]])
		if (any(c.fia)) {
			FIA.sum <- vapply(mpt.results, function(x) {if (any(grepl("^FIA$", colnames(x[["information.criteria"]][["sum"]])))) x[["information.criteria"]][["sum"]][["FIA"]] else NA}, 0)
			delta.FIA.sum <- FIA.sum - min(FIA.sum, na.rm = TRUE)
			FIA.aggregated <- vapply(mpt.results, function(x) {if (any(grepl("^FIA$", colnames(x[["information.criteria"]][["aggregated"]])))) x[["information.criteria"]][["aggregated"]][["FIA"]] else NA}, 0)
			delta.FIA.aggregated <- FIA.aggregated - min(FIA.aggregated, na.rm = TRUE)
			FIAs <- vapply(mpt.results, function(x) {if (any(grepl("^FIA$", colnames(x[["information.criteria"]][["individual"]])))) x[["information.criteria"]][["individual"]][["FIA"]] else rep(NA, n.data)}, rep(0, n.data))
			FIA.best <- rowSums(apply(FIAs, 1, function(x) x == min(x, na.rm = TRUE)))
		}
		AIC.sum <- sapply(mpt.results, function(x) x[["information.criteria"]][["sum"]][["AIC"]])
		BIC.sum <- sapply(mpt.results, function(x) x[["information.criteria"]][["sum"]][["BIC"]])
		AIC.aggregated <- sapply(mpt.results, function(x) x[["information.criteria"]][["aggregated"]][["AIC"]])
		BIC.aggregated <- sapply(mpt.results, function(x) x[["information.criteria"]][["aggregated"]][["BIC"]])
		delta.AIC.sum <- AIC.sum - min(AIC.sum)
		denom.wAIC.sum <- sum(exp(-0.5*(delta.AIC.sum)))
		wAIC.sum <- sapply(delta.AIC.sum, function(x) exp(-0.5*(x))/denom.wAIC.sum)
		delta.BIC.sum <- BIC.sum - min(BIC.sum)
		denom.wBIC.sum <- sum(exp(-0.5*(delta.BIC.sum)))
		wBIC.sum <- sapply(delta.BIC.sum, function(x) exp(-0.5*(x))/denom.wBIC.sum)
		delta.AIC.aggregated <- AIC.aggregated - min(AIC.aggregated)
		denom.wAIC.aggregated <- sum(exp(-0.5*(delta.AIC.aggregated)))
		wAIC.aggregated <- sapply(delta.AIC.aggregated, function(x) exp(-0.5*(x))/denom.wAIC.aggregated)
		delta.BIC.aggregated <- BIC.aggregated - min(BIC.aggregated)
		denom.wBIC.aggregated <- sum(exp(-0.5*(delta.BIC.aggregated)))
		wBIC.aggregated <- sapply(delta.BIC.aggregated, function(x) exp(-0.5*(x))/denom.wBIC.aggregated)
		AICs <- sapply(mpt.results, function(x) x[["information.criteria"]][["individual"]][["AIC"]])
		AIC.best <- rowSums(apply(AICs, 1, function(x) x == min(x)))
		BICs <- sapply(mpt.results, function(x) x[["information.criteria"]][["individual"]][["BIC"]])
		BIC.best <- rowSums(apply(BICs, 1, function(x) x == min(x)))
		df.out <- data.frame(model = m.names, n.parameters)
		if (any(c.fia)) {
			df.out <- cbind(df.out, delta.FIA.sum, FIA.best)
			if (output[1] != "standard") df.out <- cbind(df.out, FIA.sum)
		}
		if (output[1] == "standard") {
			df.out <- cbind(df.out, delta.AIC.sum, wAIC.sum, AIC.best, delta.BIC.sum, wBIC.sum, BIC.best)   
		} else {
			df.out <- cbind(df.out, delta.AIC.sum, wAIC.sum, AIC.best, AIC.sum, delta.AIC.aggregated, wAIC.aggregated, AIC.aggregated, delta.BIC.sum, wBIC.sum, BIC.best, BIC.sum, delta.BIC.aggregated, wBIC.aggregated, BIC.aggregated)
		}
	}
	
	rownames(df.out) <- NULL
	
	df.out[,-1] <- round(df.out[,-1], round.digit)
	df.out
	
}

