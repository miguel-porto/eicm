plot.eicm.matrix <- function(x, true.model=NULL, type=ifelse(is.null(true.model), "network", "coefficients"),
	labels=TRUE, exclude.orphans=TRUE, ...) {
	switch(pmatch(type, c("coefficients", "network")), {
		coefficientComparisonPlot(x, true.model=true.model, ...)
	}, {
		if(!is.null(true.model) && !inherits(true.model, "eicm")) stop("true.model must be of class 'eicm'")
		
		plotNetworkFromMatrix(x$sp, true.adjacency=if(is.null(true.model)) NULL else true.model$model$sp, exclude.orphans=exclude.orphans, labels=labels, ...)
	})
}


#' Plot graphs from adjacency matrices
#'
#' Plots a graph from a weighted adjacency matrix, using \code{igraph}'s plotting functions,
#' optionally comparing it with another "true" adjacency matrix.
#'
#' When comparing two adjacency matrices
#'
#' @param adjacency the square adjacency matrix.
#' @param true.adjacency optional. The reference "true" square adjacency matrix to which to compare the first one.
#' @param exclude.orphans logical. Hide nodes without links?
#' @param labels logical. Draw default labels?
#' @param lwd a numerical value giving the amount by which the arrows' line width should be magnified
#'        relative to the default, when plotting the weighted graph (only used when
#'        \code{true.adjacency} is not provided).
#' @param edge.arrow.size size of the arrow heads. See \code{\link{igraph::plot.igraph}}
#' @param severe.threshold the absolute threshold above which the interaction weights are highlighted in the graph.
#'
#' @note The arrow direction depicts the direction of the interaction. Species in columns affect species in rows.
#'
#' The matrices should include row and column labels, otherwise the node labels may not correspond to the species index
#' (when \code{exclude.orphans = TRUE})
#'
#' @return The corresponding igraph network, invisibly.
#'
#' @examples
#' # generate two adjacency matrices with 15 species and 10 interactions
#' A <- matrix(0, ncol=15, nrow=15)
#' A[sample(length(A), 10)] <- runif(10)
#' 
#' B <- matrix(0, ncol=15, nrow=15)
#' B[sample(length(B), 10)] <- runif(10)
#' 
#' # set the species names
#' rownames(A) <- rownames(B) <-
#'   colnames(A) <- colnames(B) <- paste0("S", 1:15)
#' 
#' plotNetworkFromMatrix(A, B)
#' 
#' @export
plotNetworkFromMatrix <- function(adjacency, true.adjacency=NULL, labels=TRUE, exclude.orphans=TRUE,
	lwd=1, edge.arrow.size=0.8, severe.threshold=0.5) {
	if (!requireNamespace("igraph", quietly = TRUE)) {
		message("For a network plot, install 'igraph'")
		return(invisible(NULL))
	}
	diag(adjacency) <- 0
	
	if(is.null(true.adjacency)) {
		if(exclude.orphans) {
			exclude <- apply(adjacency != 0, 1, sum) == 0 & apply(adjacency != 0, 2, sum) == 0
			adjacency <- adjacency[!exclude, !exclude]
		}
		adjacency <- t(adjacency)
		if(ncol(adjacency) < 2) return(invisible(NULL))
		net <- igraph::graph_from_adjacency_matrix(as.matrix(adjacency), mode="directed", weighted="coef", diag=FALSE)

		oldpar <- graphics::par(mar=rep(0, 4))	# we don't want to restart a new layout each call, so we only save what we change
		on.exit(graphics::par(oldpar))
		if(is.null(igraph::E(net)$coef)) {
			igraph::plot.igraph(net, edge.arrow.size=edge.arrow.size, vertex.shape=ifelse(labels, "none", "circle"),
				edge.label.family="Sans", vertex.size=5,
				layout=igraph::layout_with_fr,
				vertex.label=if(labels) names(igraph::V(net)) else NA, vertex.label.color="black")
		} else {
			igraph::E(net)$width <- (abs(igraph::E(net)$coef) + 0.4) * lwd
			igraph::plot.igraph(net, edge.arrow.size=edge.arrow.size, vertex.shape=ifelse(labels, "none", "circle"),
				edge.color=ifelse(igraph::E(net)$coef > 0, "#5555ffff", "#ff5555ff"),
				edge.label.family="Sans", vertex.size=5,
				layout=igraph::layout_with_fr,
				vertex.label=if(labels) names(igraph::V(net)) else NA, vertex.label.color="black")
		}
		return(invisible(net))
	} else {
		diag(true.adjacency) <- 0
		if(exclude.orphans) {
			merge <- true.adjacency + adjacency
			exclude <- apply(merge != 0, 1, sum) == 0 & apply(merge != 0, 2, sum) == 0
			merge <- merge[!exclude, !exclude]
			true.adjacency <- true.adjacency[!exclude, !exclude]
			adjacency <- adjacency[!exclude, !exclude]
		}
		
		truevals <- foldMatrix(true.adjacency)
		fm1 <- truevals != 0
		fm2 <- foldMatrix(adjacency) != 0

		merge <- fm1
		merge[fm1 & fm2] <- 3
		merge[fm1 & !fm2] <- 2
		merge[!fm1 & fm2] <- 1
		
		merge[fm1 & !fm2 & abs(truevals) < severe.threshold] <- 4
		
		net <- igraph::graph_from_adjacency_matrix(t(merge), mode="upper", weighted="coef", diag=FALSE)

		oldpar <- graphics::par(mar=c(0, 0, 2, 0))
		on.exit(graphics::par(oldpar))
		igraph::plot.igraph(net, edge.arrow.size=edge.arrow.size, vertex.shape=ifelse(labels, "none", "circle"),
			vertex.color="black", vertex.frame.color="white",
			edge.color=c("red", "gray", "#00aa00", "gray")[igraph::E(net)$coef],
			edge.width=c(2, 2, 2, 1)[igraph::E(net)$coef],
			edge.lty=c(1, 1, 1, 3)[igraph::E(net)$coef],
			edge.label.family="Sans", vertex.size=4, 
			layout=igraph::layout_with_fr,
			vertex.label=if(labels) names(igraph::V(net)) else NA, vertex.label.color="black")

		if(labels) {			
			graphics::legend("topleft", inset=c(0, 0), horiz=FALSE, legend=c(
				"True positives",
				sprintf("False negatives (|x|>=%.1f)", severe.threshold),
				sprintf("False negatives (|x|<%.1f)", severe.threshold),
				"False positives"),
				col=c("#009900", "gray", "gray", "red"), lwd=c(2, 2, 1, 2), lty=c(1, 1, 3, 1), bty="n")
		}
		return(invisible(net))
	}
}
