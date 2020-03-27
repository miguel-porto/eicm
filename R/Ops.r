"+.eicm.options" <- function(a, b) {
	n <- names(b)
	if("mask" %in% n) {
		# The + operator is equivalent to combine masks with AND
		if("env" %in% names(b$mask)) {
			if(is.scalar(a$mask$env)) {
				if(a$mask$env == 1L)
					a$mask$env <- b$mask$env
				# else remain all zero
			} else
				a$mask$env <- a$mask$env & b$mask$env
		}

		if("sp" %in% names(b$mask)) {
			if(is.scalar(a$mask$sp)) {
				if(a$mask$sp == 1L)
					a$mask$sp <- b$mask$sp
				# else remain all zero
			} else
				a$mask$sp <- a$mask$sp & b$mask$sp
		}
			
		excess <- names(b$mask)[is.na(pmatch(names(b$mask), c("env", "sp")))]
		if(length(excess) > 0)
			warning("These mask elements were not recognised and ignored: ", excess)

	}
	if("offset" %in% n) {	# we just silently replace offsets
		if("env" %in% names(b$offset))
			a$offset$env <- b$offset$env
		if("sp" %in% names(b$offset))
			a$offset$sp <- b$offset$sp
	}
	return(eicm.options(a))
}
