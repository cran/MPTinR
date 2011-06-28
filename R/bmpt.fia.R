
bmpt.fia <- function(s, parameters, category, N, ineq0 = NULL, Sample = 2e+05) {

	t0 <- Sys.time()
	print(paste("Computing FIA: Iteration begins at ", t0, sep = ""))
	flush.console()

	s <- strsplit(s, "")

	s <- s[[1]]

	if (any(!(s %in% c("p", "C")))) stop("Characters in the string must be either p or C.")

	type <- ifelse(s == "p", 1, 0)

	if (sum(type) != length(parameters)) stop("number of parameters not consistent in input arguments.")

	fixid <- (parameters <= 0)
	if (any(fixid)) parameter0 <- -parameters[fixid]
	else parameter0 <- 0
	parameters <- parameters[!fixid]

	if(any(floor(parameters) != parameters)) stop("use positive integers for free parameter assignment.")
	if(any(parameter0 < 0 | parameter0 > 1)) stop ("use negative values between 0 and 1 to fix a parameter.")
	if(sum(!type)!=length(category)) stop("number of categories not consistent in input arguments.")
	if (!is.null(ineq0)) if (dim(ineq0)[2] != 2) stop("ineq0 should have two columns.")

	######################

	L <- length (s)
	code <- matrix(0, L, L)

	if (type[1] == 0 & L == 1) return()
	if (type[1] == 0 & L != 1) stop("Not a BMPT model.")

	p <- 1
	u <- 1
	for (i in 2:L) {
		code[i,] <- code[p,]
		code[i,p] <- u
		if (type[i] == 1) {
			u <- 1
			p <- i
		} else {
			u <- -1
			ind <- i-1
			while (ind > 0) {
				if (ind <= 0 & i < L) stop( "Not a BMPT model.")
				
				if (type[ind] == 1) {
					p <- ind
					break
				} else {
					if (type[ind] == 0) {
						if (type[ind-1] !=1) stop("Not a BMPT model.")
						type[c(ind-1,ind, ind+1)] <- -1
						ind <- ind-2
					} else {
						if (type[ind] == -1) {
							type[ind+1] <- -1
							while (type[ind] == -1) ind <- ind-1
							if (type[ind] != 1) stop ("Not a BMPT model.")
							type[ind] <- -1
							ind <- ind-1
						}
					}
				}
			}
		}
	}

	if (ind > 0) stop ("Not a BMPT model.")

	code <- code[s == "C", s == "p"]

	##################

	code1 <- code[,!fixid, drop = FALSE]
	code0 <- code[,fixid, drop = FALSE]

	###################

	P1 <- length(parameters)
	assignment <- matrix(0, P1, P1)
	id = 1:P1
	pos <- matrix(0, 1, P1)

	i <- 1
	while (length(parameters) > 0) {
		a <- min(parameters)
		ind <- (parameters==a)
		assignment[id[ind],i] <- 1
		parameters <- parameters[-which(ind==TRUE)]
		id <- id[-which(ind==TRUE)]
		pos[i] <- a
		i <- i+1
	}
	assignment <- assignment[,1:(i-1)]
	pos <- pos[1:(i-1)]

	A <- (code1 == 1) %*% assignment 
	B <- (code1 == -1) %*% assignment 

	##################
	# fixed parameters
	# (this part was slightly changed for dealing with cases without fixed parameters, v 0.6.3)
	if (any(fixid)) {
		As <- matrix(1, dim(code0)[1], 1) %*% parameter0
		c <- apply((As^(code0==1)),1,prod)*apply((1-As)^(code0==-1),1,prod)
	} else {
		c <- rep(1, dim(code)[1])
	}

	#################
	if (!is.null(ineq0)) {
		nineq <- dim(ineq0)[1]
		Ineq <- matrix(0, nineq, dim(assignment)[2])
		for (i in 1:nineq) {
			Ineq[i, pos==ineq0[i,1]] <- -1  ####### erased transpose operator from originial code
			Ineq[i, pos==ineq0[i,2]] <- 1  ####### erased transpose operator from original code
		}
	}

	###

	t1 <- dim(A)
	M <- t1[1]
	S <- t1[2]
	C <- max(category)
	pattern<-matrix(0, C, M)
	for (i in 1:C) {
		if (!any(category==i)) stop("argument category should involve consecutive numbers from 1 to the number of categories")
		pattern[i, category==i] <- 1
	}

	sample <- 1
	integral <- 0
	vr <- 0
	count <- 0
	while (sample <= Sample) {
		theta <- rbeta(S, 0.5, 0.5)
		if (is.vector(theta) & length(theta) == 1) theta <- as.matrix(theta)
		count <- count +1
		if (!is.null(ineq0)) ineqeff <- (theta%*%(t(Ineq))<0)
		else ineqeff <- FALSE
		if (any(ineqeff) | any(theta == 0 | theta == 1)) next
		############
		# calculate the integrant, which is part of the Fisher information
		Theta <- matrix(1, M, 1) %*% theta
		p <- apply(Theta^A,1,prod)*apply((1-Theta)^B,1,prod)*c
		#if (dim(p)[2] == 1) V <- pattern%*%diag(p[,1], dim(p)[1], dim(p)[1])
		#else V <- pattern%*%diag(p)
		V <- pattern%*%diag(p)
		pc <- rowSums(V)
		delta0 <- V %*% (A-(A+B) %*% diag(theta))
		D = 1/pc
		D[is.infinite(D)] <- 0
		I <- t(delta0) %*% diag(D) %*% delta0 %*% solve(diag(theta*(1-theta))) *pi*pi
		
		detI <- abs(det(I))
		integral <- integral + sqrt(detI)
		vr <- vr + detI
		sample <- sample + 1
		if (floor(sample/10000) == sample/10000) {
			message(paste("Samples:", sample), "\r", appendLF=FALSE)
			flush.console()
		}
	}

	integral <- integral/Sample
	vr <- vr/Sample
	vr1 <- abs(vr-integral^2)/Sample
	lnInt <- log(integral)
	d1 <- sqrt(vr1)/integral*qnorm(.975)
	CI1 <- c(lnInt-d1, lnInt+d1)

	const <- Sample/count
	lnconst <- log(const)
	d2 <- sqrt((1-const)/Sample)* qnorm(.975)
	CI2 <- c(lnconst-d2, lnconst+d2)
	CFIA <- lnInt + lnconst+S/2*log(N/2/pi)
	d <- sqrt(d2^2 + d1^2)
	CI = c(CFIA-d, CFIA+d)

	out <- c(CFIA, CI, lnInt, CI1, lnconst, CI2)
	names(out) <- c("CFIA", "CI.l", "CI.u", "lnInt", "CI.lnint.l", "CI.lnint.u", "lnconst", "CI.lnconst.l", "CI.lnconst.u")
	
	t1 <- Sys.time()
	print(paste("Computing FIA: Iteration stopped at ", t1, sep = ""))
	print(t1-t0)
	
	out

}

