plot.eicm.matrix <- function(x, true.model, type=ifelse(is.null(true.model), "network", "coefficients"),
	labels=TRUE, exclude.orphans=TRUE, ...) {
	switch(pmatch(type, c("coefficients", "network")), {
		coefficientComparisonPlot(x, true.model=true.model, ...)
	}, {
		if (!requireNamespace("igraph", quietly = TRUE)) {
			message("For a network plot, install 'igraph'")
			return(invisible(NULL))
		}
		
		if(missing(true.model)) {
			distMatrix <- x$sp
			if(exclude.orphans) {
				exclude <- apply(distMatrix != 0, 1, sum) == 0 & apply(distMatrix != 0, 2, sum) == 0
				distMatrix <- distMatrix[!exclude, !exclude]
			}
			distMatrix <- t(distMatrix)
			if(ncol(distMatrix) < 2) return(invisible(NULL))
			net <- igraph::graph_from_adjacency_matrix(as.matrix(distMatrix), mode="directed", weighted="coef", diag=FALSE)

			oldpar <- graphics::par(no.readonly = TRUE)
			on.exit(graphics::par(oldpar))

			graphics::par(mar=rep(0, 4))
			if(is.null(igraph::E(net)$coef)) {
				igraph::plot.igraph(net, edge.arrow.size=0.8, vertex.shape=ifelse(labels, "none", "circle"),
					edge.label.family="Sans", vertex.size=5,
					layout=igraph::layout_with_fr,
					vertex.label=if(labels) names(igraph::V(net)) else NA, vertex.label.color="black")
			} else {
				igraph::E(net)$width <- (abs(igraph::E(net)$coef) + 0.4) * 2
				igraph::plot.igraph(net, edge.arrow.size=0.8, vertex.shape=ifelse(labels, "none", "circle"),
					edge.color=ifelse(igraph::E(net)$coef > 0, "#5555ffff", "#ff5555ff"),
					edge.label.family="Sans", vertex.size=5,
					layout=igraph::layout_with_fr,
					vertex.label=if(labels) names(igraph::V(net)) else NA, vertex.label.color="black")
			}
		} else {
			if(!inherits(true.model, "eicm")) stop("true.model must be of class 'eicm'")
			plotBiNetworkFromMatrices(x$sp, true.model$model$sp, labels=labels, exclude.orphans=exclude.orphans, ...)
		}
	})
}

plotBiNetworkFromMatrices <- function(fitted.distMatrix, true.distMatrix, labels=TRUE, exclude.orphans=TRUE,
	severe.threshold=0.5, layout=TRUE) {
	if(exclude.orphans) {
		merge <- true.distMatrix + fitted.distMatrix
		exclude <- apply(merge != 0, 1, sum) == 0 & apply(merge != 0, 2, sum) == 0
		merge <- merge[!exclude, !exclude]
		true.distMatrix <- true.distMatrix[!exclude, !exclude]
		fitted.distMatrix <- fitted.distMatrix[!exclude, !exclude]
	}
	
	truevals <- foldMatrix(true.distMatrix)
	fm1 <- truevals != 0
	fm2 <- foldMatrix(fitted.distMatrix) != 0

	merge <- fm1
	merge[fm1 & fm2] <- 3
	merge[fm1 & !fm2] <- 2
	merge[!fm1 & fm2] <- 1
	
	merge[fm1 & !fm2 & abs(truevals) < severe.threshold] <- 4
	
	net <- igraph::graph_from_adjacency_matrix(t(merge), mode="upper", weighted="coef", diag=FALSE)

	if(layout) {	
		oldpar <- graphics::par(no.readonly = TRUE)
		on.exit(graphics::par(oldpar))
		graphics::par(mar=c(0, 0, 2, 0))
	}
	igraph::plot.igraph(net, edge.arrow.size=0.5, vertex.shape=ifelse(labels, "none", "circle"),
		vertex.color="black", vertex.frame.color="white",
		edge.color=c("red", "gray", "#00aa00", "gray")[igraph::E(net)$coef],
		edge.width=c(2, 2, 2, 1)[igraph::E(net)$coef],
		edge.lty=c(1, 1, 1, 3)[igraph::E(net)$coef],
		edge.label.family="Sans", vertex.size=4, 
		layout=igraph::layout_with_fr,
		vertex.label=if(labels) gsub("sp", "", names(igraph::V(net))) else NA, vertex.label.color="black")
	return(net)
}

