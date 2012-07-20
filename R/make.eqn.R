
make.eqn <- function(model.filename, eqn.filename) {
	model <- .read.MPT.model(model.filename)
	model.df <- .make.model.df(model)
	write.table(dim(model.df)[1], eqn.filename, row.names = FALSE, col.names = FALSE)
	write.table(model.df, eqn.filename, append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
}

make.mdt <- function(data, mdt.filename, index, prefix = "dataset") {
	my.con <- file(mdt.filename, open = "w")
	if (is.vector(data)) { 
		df <- data.frame(seq_len(length(data)), data)
		colnames(df) <- c(prefix, index)
		write.table(df, file = my.con, row.names = FALSE, quote = FALSE)
		writeLines("===", con = my.con)
	}
	if (is.matrix(data) | is.data.frame(data)) {
		for (c in seq_len(nrow(data))) {
			df <- data.frame(seq_len(length(data[c,])), data[c,])
			colnames(df) <- c(prefix, c)
			if (c == 1) {
				suppressWarnings(write.table(df, file = my.con, row.names = FALSE, quote = FALSE))
				writeLines("===", con = my.con)
			} else {
				suppressWarnings(write.table(df, file = my.con, row.names = FALSE, quote = FALSE, append = TRUE))
				writeLines("===", con = my.con)
			}
		}
	}
	close(my.con)
}
