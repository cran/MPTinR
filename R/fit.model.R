   
fit.model <- function(data, model.filename, restrictions.filename = NULL, n.optim = 5, fia = NULL, ci = 95, starting.values = NULL, lower.bound = 0, upper.bound = 1, output = c("standard", "fia", "full"), reparam.ineq = TRUE, fit.aggregated = TRUE, sort.param = TRUE, show.messages = TRUE, model.type = c("easy", "eqn", "eqn2"),  multicore = c("none", "individual", "n.optim"), sfInit = FALSE, nCPU = 2, control = list(), use.gradient = TRUE, use.hessian = FALSE){
	
	llk.model <- function(Q, unlist.model, data, param.names, n.params, lower.bound, upper.bound, llk.gradient, llk.hessian, tmp.env){
		if (length(upper.bound) == 1) {
			Q[Q > upper.bound] <- upper.bound
		} else {
			Q[Q > upper.bound] <- upper.bound[Q > upper.bound]
		}
		if (length(lower.bound) == 1) {
			Q[Q < lower.bound] <- lower.bound
		} else {
			Q[Q < lower.bound] <- lower.bound[Q < lower.bound]
		}
		#tmpllk.env <- new.env()
		for (i in seq_len(n.params))  assign(param.names[i],Q[i], envir = tmp.env)
		
		model.eval <- vapply(unlist.model, eval, envir = tmp.env, 0)
		if (any(model.eval < 0)) stop(paste("Model not constructed well. Line ", which(model.eval < 0), " produces probabilities < 0!", sep = ""))
		llk <- data * log(model.eval)
		llk[data == 0] <- 0
		llk <- sum(llk)
		if (is.na(llk)) llk <- -1e10
		if (llk == -Inf) llk <- -1e10
		return(-llk)
	}
	
	model.predictions <- function(Q, unlist.model, param.names, n.params, lower.bound, upper.bound, llk.gradient, llk.hessian, tmp.env){
		#tmpllk.env <- new.env()
		for (i in seq_len(n.params))  assign(param.names[i],Q[i], envir = tmp.env)
		vapply(unlist.model, eval, envir = tmp.env, 0)
	}
	
	llk.gradient.funct <- function(Q, unlist.model, data, param.names, n.params, lower.bound, upper.bound, llk.gradient, llk.hessian, tmp.env){
		if (length(upper.bound) == 1) {
			Q[Q > upper.bound] <- upper.bound
		} else {
			Q[Q > upper.bound] <- upper.bound[Q > upper.bound]
		}
		if (length(lower.bound) == 1) {
			Q[Q < lower.bound] <- lower.bound
		} else {
			Q[Q < lower.bound] <- lower.bound[Q < lower.bound]
		}
		#tmpllk.env <- new.env()
		for (i in 1:n.params)  assign(param.names[i],Q[i], envir = tmp.env)
		
		model.eval <- vapply(llk.gradient, eval, 0, envir = tmp.env)
		model.eval[is.na(model.eval)] <- -1e10
		model.eval[model.eval == -Inf] <- -1e10
		return(-model.eval)
	}

	llk.hessian.funct <- function(Q, unlist.model, data, param.names, n.params, lower.bound, upper.bound, llk.gradient, llk.hessian, tmp.env){
		if (length(upper.bound) == 1) {
			Q[Q > upper.bound] <- upper.bound
		} else {
			Q[Q > upper.bound] <- upper.bound[Q > upper.bound]
		}
		if (length(lower.bound) == 1) {
			Q[Q < lower.bound] <- lower.bound
		} else {
			Q[Q < lower.bound] <- lower.bound[Q < lower.bound]
		}
		#tmpllk.env <- new.env()
		for (i in 1:n.params)  assign(param.names[i],Q[i], envir = tmp.env)
		
		model.eval <- apply(llk.hessian, c(1,2), function(x) eval(x[[1]], envir = tmp.env))
		#model.eval <- apply(llk.hessian, c(1,2), eval, envir = tmpllk.env)
		model.eval[is.na(model.eval)] <- -1e10
		model.eval[model.eval == -Inf] <- -1e10
		return(-model.eval)
	}
	
	###########################################################
	## objective, gradient, hessian above, preparation below ##
	###########################################################
	
	# check if model or restrictions are connections and if it is needed later again:
	class.model <- class(model.filename)
	if (!is.null(fia) & ("connection" %in% class.model)) {
		tmp.model <- readLines(model.filename)
		model.filename <- textConnection(tmp.model)
	}
	class.restr <- class(restrictions.filename)
	if (!is.null(restrictions.filename) & !is.null(fia) & ("connection" %in% class.restr)) {
		tmp.restr <- readLines(restrictions.filename)
		restrictions.filename <- textConnection(tmp.restr)
	}
	
	
	tree <- .get.mpt.model(model.filename, model.type)
	if(is.null(data)) stop("Model seems to be constructed well (i.e., all probabilities sum to 1), but data is NULL.")
	
	if(is.vector(data)) {
		data <- array(data, dim = c(1, length(data)))
		multiFit <- FALSE
	} else
		if(dim(data)[1] == 1) {
			if (is.data.frame(data)) data <- as.matrix(data)
			data <- array(data, dim = c(1,length(data)))
			multiFit <- FALSE
		} else 
			if(is.matrix(data) | is.data.frame(data)) {
				if (is.data.frame(data)) data <- as.matrix(data)
				multiFit <- TRUE
			} else stop("data is neither vector, nor matrix, nor data.frame!")
	
	if (sum(sapply(tree, length)) != length(data[1,])) stop(paste("Size of data does not correspond to size of model (i.e., model needs ", sum(sapply(tree, length)), " datapoints, data gives ", length(data[1,]), " datapoints).", sep = ""))
	
	orig.params <- NULL
	use.restrictions <- FALSE
		
	if (!is.null(restrictions.filename)) {
		orig.params <- .find.MPT.params(tree)
		new.restrictions <- .check.restrictions(restrictions.filename, tree)
		if (length(new.restrictions) > 0) use.restrictions <- TRUE
		if (!reparam.ineq) {
			res.no.ineq <- new.restrictions
			for (res in 1:length(new.restrictions)) if (new.restrictions[[res]][3] == "<") res.no.ineq[[1]] <- NULL
			if (length(res.no.ineq) == 0) use.restrictions <- FALSE
			else new.restrictions <- res.no.ineq
			}
		if (use.restrictions) tree <- .apply.MPT.restrictions(tree, new.restrictions)
		restrictions <- new.restrictions
	}
	
	# check if restrictions is connection and needed again then construct anew here:
	if (!is.null(fia) & ("connection" %in% class.model)) {
		model.filename <- textConnection(tmp.model)
	}
	if (!is.null(restrictions.filename) & !is.null(fia) & ("connection" %in% class.restr)) {
		restrictions.filename <- textConnection(tmp.restr)
	}
	
	# make arguments:
	
	param.names <- .find.MPT.params(tree)
	length.param.names <- length(param.names)
	categories.per.type <- vapply(tree, length, 0)
	
	# gradient and hessian:
	
	llk.function <- tryCatch(.make.llk.function(tree, param.names, length.param.names), error = function(e) {warning("likelihood function cannot be build"); NULL})
	llk.gradient <- tryCatch(.make.llk.gradient(llk.function, param.names, length.param.names), error = function(e) {warning("gradient function cannot be build (probably derivation failure, see ?D"); NULL})
	llk.hessian <- tryCatch(.make.llk.hessian(llk.function, param.names, length.param.names), error = function(e) {warning("hessian function cannot be build (probably derivation failure, see ?D"); NULL})
	
	if (!is.null(fia)) {
		if (multiFit) {
			data.new <- rbind(data, apply(data,2,sum))
			fia.tmp <- get.mpt.fia(data.new, model.filename, restrictions.filename, fia, model.type)
			fia.df <- fia.tmp[-dim(fia.tmp)[1],]
			fia.agg.tmp <- fia.tmp[dim(fia.tmp)[1],]
			fia.df <- list(fia.df, fia.agg.tmp)
		} else {
			fia.df <- get.mpt.fia(data, model.filename, restrictions.filename, fia, model.type)
		}
	}
	
	# call the workhorse:	
	fit.mptinr(data = data, objective = llk.model, param.names = param.names, categories.per.type = categories.per.type, gradient = llk.gradient.funct, use.gradient = use.gradient, hessian = llk.hessian.funct, use.hessian = use.hessian, prediction = model.predictions, n.optim = n.optim, fia.df = if(!is.null(fia)) fia.df, ci = ci, starting.values = starting.values, lower.bound = lower.bound, upper.bound = upper.bound, output = output, orig.params = orig.params, fit.aggregated = fit.aggregated, sort.param = sort.param, show.messages = show.messages, use.restrictions = use.restrictions, restrictions = restrictions, multicore = multicore, sfInit = sfInit, nCPU = nCPU, control = control, unlist.model = unlist(tree), llk.gradient = llk.gradient, llk.hessian = llk.hessian)
	
}

