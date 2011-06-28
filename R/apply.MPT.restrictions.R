
.apply.MPT.restrictions <- function(tree, restrictions) {
	safeDeparse <- function(expr){
		ret <- paste(deparse(expr), collapse="")
		#rm whitespace
		ret <- gsub("[[:space:]][[:space:]]+", " ", ret)
		gsub("^expression\\(", "", gsub("[[:space:]][[:space:]]", " ", gsub("\\)$", "", ret)))
	}
	for (c1 in 1:length(tree)){
		for (c2 in 1:length(tree[[c1]])) {
			for (c3 in 1:length(restrictions)) {
				tree[[c1]][c2] <- parse(text = gsub(paste("\\<",restrictions[[c3]][1], "\\>", sep =""),restrictions[[c3]][2],safeDeparse(tree[[c1]][c2])))[1]
			}
		}
	}
	return(tree)
}
