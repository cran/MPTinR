
make.eqn <- function(model.filename, eqn.filename) {
	model <- .read.MPT.model(model.filename)
	model.df <- .make.model.df(model)
	write.table(dim(model.df)[1], eqn.filename, row.names = FALSE, col.names = FALSE)
	write.table(model.df, eqn.filename, append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
}