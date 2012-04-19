
make.eqn <- function(model.filename, eqn.filename) {
	model <- .read.MPT.model(model.filename)
	model.df <- .make.model.df(model)
	write.table(dim(model.df)[1], eqn.filename, row.names = FALSE, col.names = FALSE)
	write.table(model.df, eqn.filename, append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
}

make.mdt <- function(data, mdt.filename, index, prefix = "dataset") {
	if (!is.vector(data)) stop("data needs to be a vector")
	df <- data.frame(seq_len(length(data)), data)
	colnames(df) <- c(prefix, index)
	write.table(df, file = mdt.filename, row.names = FALSE, quote = FALSE)
}
